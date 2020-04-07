
close all;
clear;
clc;

animal = 'MK012';
expDate = '2014-08-10';
iExp = 1812;
Plane = 2;

folder = fullfile(ppbox.getTmpFolder, animal, expDate, num2str(iExp));

[fileName, pathName] = uigetfile(sprintf('*plane%03.0f*ROI*.mat', Plane), [], folder, 'MultiSelect', 'on');
disp('Let''s analyze the following files:');
disp(pathName)
disp(fileName')
if ~iscell(fileName)
    fileName = {fileName};
end

nFiles = length(fileName);

for iFile = 1:nFiles
    data = load(fullfile(pathName, fileName{iFile}));
    info(iFile) = data.meta;
end

%% defining the sequence of colors
colors = 'rgbcmyw';
colorsMatrix = [1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1; 1 1 0; 1 1 1];

%% Lets see what area the sub-datasets cover of the overall targetFrame

figure
% subplot(1, 2, 1)
imshowpair(info(1).targetFrame, info(1).collage);
axis equal tight off
colormap gray(256);
hold on;

for iFile = 1:nFiles
    pos = info(iFile).cropPosDest;
    xx{iFile} = [pos.xmin, pos.xmin, pos.xmax, pos.xmax, pos.xmin];
    yy{iFile} = [pos.ymin, pos.ymax, pos.ymax, pos.ymin, pos.ymin];
    iColor = mod(iFile-1, length(colors))+1;
    plot(xx{iFile}, yy{iFile}, colors(iColor));
end

%% Showing all the datasets separately

figure
nRows = floor(sqrt(nFiles));
nColumns = ceil(nFiles/nRows);
for iFile = 1:nFiles
    pic = info(iFile).targetFrame;
    pic = pic - min(pic(:));
    pic = pic/max(pic(:));
    %     pic = repmat(pic, 1, 1, 3);
    nROIs = length(info(iFile).targetFrameROI);
    rois = pic;
    for iROI = 1:nROIs
        rois(info(iFile).targetFrameROI{iROI}) = 1;
    end
    
    finalPic{iFile} = zeros(size(pic, 1), size(pic, 2),  3);
    iColor = mod(iFile-1, length(colors))+1;
    for iChannel = 1:3
        if colorsMatrix(iColor, iChannel);
            finalPic{iFile}(:,:,iChannel) = rois;
        else
            finalPic{iFile}(:,:,iChannel) = pic;
        end
    end
    subplot(nRows, nColumns, iFile)
    image(finalPic{iFile})
    %     colormap gray(256)
    axis equal tight off
    hold on;
    for iF = 1:nFiles
        iColor = mod(iF-1, length(colors))+1;
        plot(xx{iF}, yy{iF}, colors(iColor));
    end
end

%% One more figure

figure
mergedPic = zeros(size(finalPic{1}));
for iFile = 1:nFiles
    mergedPic = mergedPic+finalPic{iFile};
end
mergedPic = mergedPic/nFiles;
image(mergedPic);
axis equal tight off
hold on;
for iF = 1:nFiles
    iColor = mod(iF-1, length(colors))+1;
    plot(xx{iF}, yy{iF}, colors(iColor));
end

%% Making a collage, based on registered Mean Intensity Projections
% assuming we have these files, if not - crash
% these files are automatically produced by layBin.m function

regCollage = info(1).targetFrame;
for iFile = 1:nFiles
    filename = fullfile(info(iFile).folderProcessed, [info(iFile).basenameRegistered, '_AVG.tiff']);
    AVG{iFile} = double(img.loadFrames(filename, 1, 1, 1));
    regCollage(info(iFile).validY, info(iFile).validX) = AVG{iFile};
    pp = prctile(AVG{iFile}(:), [0.1, 99.9]);
    AVG{iFile}(AVG{iFile}<pp(1)) = pp(1);
    AVG{iFile}(AVG{iFile}>pp(2)) = pp(2);
    AVG{iFile} = (AVG{iFile}-pp(1))/(pp(2)-pp(1));
    %     pos = info(iFile).cropPosDest;
end

pp = prctile(regCollage(:), [0.1 99.9]);
regCollage(regCollage(:)<pp(1)) = pp(1);
regCollage(regCollage(:)>pp(2)) = pp(2);
figure
imagesc(regCollage);
axis equal tight off;
colormap gray;
hold on;

for iFile = 1:nFiles
    xmin = info(iFile).validX(1);
    xmax = info(iFile).validX(end);
    ymin = info(iFile).validY(1);
    ymax = info(iFile).validY(end);
    iColor = mod(iFile-1, length(colors))+1;
    plot([xmin, xmin, xmax, xmax, xmin], [ymin, ymax, ymax, ymin, ymin], colors(iColor));
end

%% Collecting together all the info from multiple files

colors = 'rgb';
rectIdx = ones(1, length(info(1).targetFrameROI));
cellIdx = 1:length(info(1).targetFrameROI);
for iFile = 2:nFiles
    rectIdx = [rectIdx, iFile*ones(1, length(info(iFile).targetFrameROI))];
    cellIdx = [cellIdx, 1:length(info(iFile).targetFrameROI)];
end
ROIs = [info.targetFrameROI];
F = [info.F];
localROIs = cell(0);
roiXY = cell(0);
for iFile = 1:nFiles
    localROIs = [localROIs, info(iFile).ROI.CellMaps];
    roiXY = [roiXY, info(iFile).ROI.CellXYPixels];
end
sz = size(info(1).targetFrame);

%% Interactive part of deleting the unnecessary traces

nROIs = length(ROIs);
iROI = 1;
nRows = 2;
nColumns = 4;
sqsz = 50;

figure
while iROI<nROIs
    idx = iROI;
    clf;
    for jROI = iROI+1:nROIs
        if ~isempty(intersect(ROIs{iROI}, ROIs{jROI}))
            idx = [idx; jROI];
        end
    end
    if length(idx)>1
        rois = zeros(prod(sz), 3);
        ax = [0];
        for iIdx = 1:min(3, length(idx))
            rois(ROIs{idx(iIdx)}, iIdx) = 1;
            ax(end+1) = subplot(nRows, nColumns, iIdx+1);
            mask = AVG{rectIdx(idx(iIdx))}*0.9;
            frame = repmat(mask, 1, 1, 3);
            mask(localROIs{idx(iIdx)}) = 1;
            frame(:,:,iIdx) = mask;
            xAxis = info(rectIdx(idx(iIdx))).validX;
            yAxis = info(rectIdx(idx(iIdx))).validY;
            image(xAxis, yAxis, frame);
            axis equal tight off
            colormap gray
            title({sprintf('nPixels = %d, iRect = %d', length(localROIs{idx(iIdx)}), rectIdx(idx(iIdx))); ...
                sprintf('Press ''%d'' to delete me', iIdx)});
            subplot(nRows, nColumns, [nColumns+1:2*nColumns])
            plot(F(:, idx(iIdx)), colors(iIdx))
            hold on;
        end
        tmp = corrcoef(F(:, idx(1:min(3, length(idx)))));
        try
            title(sprintf('\\rho_{12}=%5.3f, \\rho_{13}=%5.3f, \\rho_{23}=%5.3f,', tmp(1, 2), tmp(1,3), tmp(2,3)));
        catch
            title(sprintf('\\rho_{12}=%5.3f', tmp(1, 2)));
        end
        rois = reshape(rois, sz(1), sz(2), 3);
        ax(1) = subplot(2, 4, 1);
        imagesc(rois);
        axis equal tight off
        linkaxes(ax, 'xy');
        xc = sum(sum(sum(rois, 3), 1).*(1:sz(2)))/sum(rois(:));
        yc = sum(sum(sum(rois, 3), 2).*(1:sz(1))')/sum(rois(:));
        zoom reset % this allows to zoom out beyond the xlim/ylim set in the following lines
        xlim([xc-sqsz/2, xc+sqsz/2]);
        ylim([yc-sqsz/2, yc+sqsz/2]);
        title(sprintf('Press ''0'' to move on w/o deleting'));
        
        FlushEvents('keyDown');
        key = GetChar();
        indices = 1:nROIs;
        switch key
            case '0'
                % dont delete anything, and just move on
                iROI = iROI + 1;
            case '1'
                indices = setdiff(1:nROIs, idx(1));
            case '2'
                indices = setdiff(1:nROIs, idx(2));
            case '3'
                if length(idx)>=3
                    indices = setdiff(1:nROIs, idx(3));
                else
                    % do nothing2
                end
            otherwise
                % do nothing
        end
        ROIs = ROIs(indices);
        localROIs = localROIs(indices);
        roiXY = roiXY(indices);
        F = F(:, indices);
        rectIdx = rectIdx(indices);
        cellIdx = cellIdx(indices);
        nROIs = length(ROIs);
    else
        iROI = iROI + 1;
    end
end

% return;

%% Building the new merged info structure

infoMerged = info(1);
infoMerged = rmfield(infoMerged, 'basenameRegistered');
infoMerged = rmfield(infoMerged, 'basenameRect');
infoMerged = rmfield(infoMerged, 'cropPosition');
infoMerged = rmfield(infoMerged, 'iCrop');
infoMerged.nCrops = nFiles;
infoMerged.basenameCrop = cell(nFiles, 1);
infoMerged.targetFrameCrop = cell(nFiles, 1);
infoMerged.cropPosSource = struct('ymin', [], 'ymax', [], 'xmin', [], 'xmax', []);
infoMerged.cropPosDest = struct('ymin', [], 'ymax', [], 'xmin', [], 'xmax', []);
infoMerged.validX = cell(nFiles, 1);
infoMerged.validY = cell(nFiles, 1);
infoMerged.dx = cell(1, nFiles);
infoMerged.dy = cell(1, nFiles);
for iFile = 1:nFiles
    infoMerged.basenameCrop{iFile} = info(iFile).basenameRect;
    infoMerged.targetFrameCrop{iFile} = info(iFile).targetFrameCrop;
    infoMerged.cropPosSource(iFile) = info(iFile).cropPosSource;
    infoMerged.cropPosDest(iFile) = info(iFile).cropPosDest;
    infoMerged.validX{iFile} = info(iFile).validX;
    infoMerged.validY{iFile} = info(iFile).validY;
    infoMerged.dx{iFile} = info(iFile).dx;
    infoMerged.dy{iFile} = info(iFile).dy;
end
infoMerged.F = F;
infoMerged.targetFrameROI = ROIs;

fNames = fieldnames(infoMerged.ROI);
for iFile = 1:nFiles
    idx = cellIdx(find(rectIdx == iFile));
    for iField = 1:length(fNames)
        if iFile == 1
            infoMerged.ROI.(fNames{iField}) = info(iFile).ROI.(fNames{iField})(idx);
        else
            infoMerged.ROI.(fNames{iField}) = cat(2, ...
                infoMerged.ROI.(fNames{iField}), info(iFile).ROI.(fNames{iField})(idx));
        end
    end
end

for iFile = 1:nFiles
    idx = cellIdx(find(rectIdx == iFile));
    for iChannel = 1:infoMerged.nChannels
        if isfield(infoMerged.chData(iChannel), 'F')
            
            if iFile == 1
                infoMerged.chData(iChannel).F = info(iFile).chData(iChannel).F(:, idx);
            else
                infoMerged.chData(iChannel).F = cat(2, ...
                    infoMerged.chData(iChannel).F, info(iFile).chData(iChannel).F(:, idx));
            end
        end
    end
end

infoMerged.registeredMIP = cell(nFiles, 1);
for iFile = 1:nFiles
    filename = fullfile(info(iFile).folderProcessed, [info(iFile).basenameRegistered, '_AVG.tiff']);
    infoMerged.registeredMIP{iFile} = img.loadFrames(filename, 1, 1, 1);
end
infoMerged.collageMIP = regCollage;
infoMerged.mergeInfo.iCrop = rectIdx;
infoMerged.mergeInfo.iCell = cellIdx;

%% saving the new structure (under the name 'meta'

filename = fullfile(infoMerged.folderProcessed, [infoMerged.basenamePlane, '_ROI.mat']);
meta = infoMerged;
save(filename, 'meta');
