function [maskImage, colors] = plotCellMasks(meta, doSave, plotIDs)

maskImage = zeros(length(meta.validY), length(meta.validX), 'uint8');
% maskImage = zeros(512, 512, 'uint16');
coloredNeurons = find(strcmp(meta.ROI.CellClasses, 's'));
grayNeurons = setdiff(1:length(meta.ROI.CellMaps), coloredNeurons);
centroids = zeros(length(meta.ROI.CellClasses), 2);
textCols = zeros(length(meta.ROI.CellClasses), 1);
for k = 1:length(coloredNeurons)
    if isempty(meta.ROI.CellMaps{coloredNeurons(k)})
        continue
    end
    maskImage(meta.ROI.CellMaps{coloredNeurons(k)}) = coloredNeurons(k);
%     maskImage(meta.targetFrameROI{coloredNeurons(k)}) = coloredNeurons(k);
    im = zeros(length(meta.validY), length(meta.validX), 'uint8');
    im(meta.ROI.CellMaps{coloredNeurons(k)}) = 1;
    stats = regionprops(im, 'Centroid');
    centroids(coloredNeurons(k),:) = stats(1).Centroid;
    textCols(coloredNeurons(k)) = 0.3;
    if k < 0.28 * length(coloredNeurons) || k > 0.8 * length(coloredNeurons)
        textCols(coloredNeurons(k)) = 1;
    end
end
for k = 1:length(grayNeurons)
    if isempty(meta.ROI.CellMaps{grayNeurons(k)})
        continue
    end
    maskImage(meta.ROI.CellMaps{grayNeurons(k)}) = length(meta.ROI.CellClasses) + 1;
    im = zeros(length(meta.validY), length(meta.validX), 'uint8');
    im(meta.ROI.CellMaps{grayNeurons(k)}) = 1;
    stats = regionprops(im, 'Centroid');
    centroids(grayNeurons(k),:) = stats(1).Centroid;
    textCols(grayNeurons(k)) = 1;
end
figure
imagesc(maskImage)
colors = colormap(jet(min(254,length(meta.ROI.CellClasses))));
% colors = ones(length(coloredNeurons),3);
colors = [0 0 0; colors];
if ~isempty(grayNeurons)
    colors = [colors; 0.5 0.5 0.5];
%     colors = [colors; 0 0 0];
end
colormap(colors)
axis equal tight off
if plotIDs == 1
    numCols = unique(textCols);
    for c = 1:length(numCols)
        ind = find(textCols == numCols(c));
        text(centroids(ind,1), centroids(ind,2), num2str(ind), ...
            'Color', [1 1 1] * numCols(c), ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold')
%         text(centroids(ind,1), centroids(ind,2), num2str(ind), ...
%             'Color', [1 1 1] * 0, 'BackgroundColor', 'w', 'Margin', 1, ...
%             'HorizontalAlignment', 'center', 'FontWeight', 'bold')
    end
    screenSize = get(0,'ScreenSize');
    set(gcf, 'Position', [10 50 1120 screenSize(4)-140])
end
title(['ROI masks of plane ' num2str(meta.iPlane)])

if doSave == 1
    imwrite(maskImage, colors, fullfile(meta.folderProcessed, ...
        strrep(meta.basenameRegistered, 'registered', 'cellMasks.tiff')))
end