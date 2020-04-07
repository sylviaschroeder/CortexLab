function [cellClasses, cellValues, surrMasks, surrPrctiles] = ...
    getCellClasses(classImage, ROIs, ops, stat)

cellClasses = NaN(length(ROIs),1);
cellValues = NaN(length(ROIs),1);
surrPrctiles = NaN(length(ROIs),100);

if ~isfield(ops, 'zoom')
    ops.zoom = 2;
end
if ~isfield(ops, 'scope')
    ops.scope = 'b';
end
if ~isfield(ops, 'distance')
    ops.distance = 3;
end
if ~isfield(ops, 'width')
    ops.width = 50;
end
if ~isfield(ops, 'classThresholds')
    ops.classThresholds = [30 70]; % thresholds in percentiles to classify
                                   % cells as not marked, unknown, and
                                   % marked
end
if ~isfield(ops, 'bloodThreshold')
    ops.bloodThreshold = 0; % 0-100% maximum pixel value to count as blood vessel
end
if ~isfield(ops, 'bloodSize')
    ops.bloodSize = 0; % 0-100%; minimum size to count as blood vessel
end

cellMasks = zeros([length(ROIs), size(classImage)]);
for iCell = 1:length(ROIs)
    tmp = zeros(size(classImage));
    tmp(ROIs{iCell}) = 1;
    cellMasks(iCell,:,:) = tmp;
end
allMasks = squeeze(sum(cellMasks,1));

microns2Pix = infoPixUm(size(classImage,1), ops.zoom, ops.scope);
opt.innerNeuropil = round(ops.distance * microns2Pix.rPU);
opt.outerNeuropil = round(ops.width * microns2Pix.rPU);
[stat, cellPix] = createCellMasks(stat, numel(ops.yrange), numel(ops.xrange));
[~,surr] = createNeuropilMasks(opt, stat, cellPix);
surrMasks = zeros(size(surr,1), size(classImage,1), size(classImage,2));
surrMasks(:, ops.yrange, ops.xrange) = surr > 0;

bloodThresh = prctile(classImage(:), ops.bloodThreshold);
bloodInd = classImage <= bloodThresh;
bloodInd = bwareafilt(bloodInd, round([numel(classImage) * ops.bloodSize / 100, ...
    numel(classImage)]));
bloodMask = ~bwmorph(bloodInd, 'close');

% xVals = linspace(min(classImage(:)), max(classImage(:)), 100);
for iCell = 1:length(ROIs)
    cellVal = median(classImage(squeeze(cellMasks(iCell,:,:) == 1)));
    cellValues(iCell) = cellVal;
    surrVals = classImage(squeeze(surrMasks(iCell,:,:) == 1) & ...
        bloodMask);
    threshs = prctile(surrVals, ops.classThresholds);
    surrPrctiles(iCell,:) = prctile(surrVals, 1:100);
    if cellVal <= threshs(1)
        cellClasses(iCell) = -1;
    elseif cellVal >= threshs(2)
        cellClasses(iCell) = 1;
    else
        cellClasses(iCell) = 0;
    end
end