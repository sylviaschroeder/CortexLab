function [ci, bloodMask, bloodHandle, lines] = ...
    plotCellClasses(classImage, ROIs, classes, ops, ax)

if nargin < 5
    figure
    ax = gca;
end

if ~isfield(ops, 'bloodThreshold')
    ops.bloodThreshold = 0; % 0-100% maximum pixel value to count as blood vessel
end
if ~isfield(ops, 'bloodSize')
    ops.bloodSize = 0; % 0-100%; minimum size to count as blood vessel
end
if ~isfield(ops, 'colors')
    ops.colors = [0 .7 0; 0 0 1; 0.75 0 0.75];
end

classImage = double(classImage);
ci = classImage - min(classImage(:));
ci = round(ci ./ max(ci(:)) .* 255 + 1);
ci = ind2rgb(ci, hot(256));

bloodThresh = prctile(classImage(:), ops.bloodThreshold);
bloodMask = classImage <= bloodThresh;
bloodMask = bwareafilt(bloodMask, round([numel(classImage) * ops.bloodSize / 100, ...
    numel(classImage)]));
bloodMask = double(bwmorph(bloodMask, 'close'));

axes(ax);
imshow(ci)
hold on
blood = ones(size(ci(:,:,1)));
bloodHandle = imshow(blood, ones(1,3)*0.5);
alpha(bloodHandle, bloodMask);

lines = NaN(1, length(ROIs));
if ~all(isnan(classes))
    colors = repmat(ops.colors(3,:), length(ROIs),1);
    colors(classes == 1,:) = repmat(ops.colors(1,:), sum(classes == 1), 1);
    colors(classes == -1,:) = repmat(ops.colors(2,:), sum(classes == -1), 1);
    for iCell = 1:length(ROIs)
        tmp = zeros(size(classImage));
        tmp(ROIs{iCell}) = 1;
        perim = bwperim(tmp, 8);
        [yPerim, xPerim] = find(perim);
        lines(iCell) = plot(xPerim,yPerim,'.', ...
            'Color', colors(iCell,:), 'MarkerSize', 3);
    end
end