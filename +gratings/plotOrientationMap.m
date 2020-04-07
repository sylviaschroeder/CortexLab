function figureHandles = plotOrientationMap(neuronOrientations, type, ...
    ROImaps, mapSize, minWeight)

% neuronOrientations    [neurons x 2]; 1st col: pref. orientation or
%                       direction, 2nd col.: goodness of tuning, e.g. adj.
%                       R-Squared
% type                  'ori' or 'dir'
% ROImaps               {neurons x 1}; each entry contains pixels of ROI
% mapsSize              [1 x 2]; size of frame ([validY validX])
% minWeight             minimum goodness of tuning to show orientation,
%                       otherwise ROI will be gray

if nargin < 5
    minWeight = 0.3;
end

map = zeros([mapSize,3]);
figureHandles = zeros(1,1);
colors = hsv(360);
if strcmp(type, 'ori')
    neuronOrientations(:,1) = round(mod(neuronOrientations(:,1),180)*2);
    str = 'Orientation';
elseif strcmp(type, 'dir')
    neuronOrientations(:,1) = round(mod(neuronOrientations(:,1),360));
    str = 'Direction';
else
    display('type must be ''ori'' or ''dir''!')
    return
end
for roi = 1:size(neuronOrientations,1)
    ori = neuronOrientations(roi,1);
    if ori == 0
        ori = 360;
    end
    col = colors(ori,:);
%     col = col .* neuronOrientations(roi,2) + [1 1 1].*.5 .* (1-neuronOrientations(roi,2));
%     col = col .* neuronOrientations(roi,2);
    if neuronOrientations(roi,2) < minWeight
        tmp = zeros(mapSize);
        tmp(ROImaps{roi}) = 1;
        ROImaps{roi} = find(bwperim(tmp));
        col = ones(1,3).*.5;
    end
    map(ROImaps{roi}) = col(1);
    map(prod(mapSize)+ROImaps{roi}) = col(2);
    map(2*prod(mapSize)+ROImaps{roi}) = col(3);
end

figureHandles(1) = figure('Position', [720 275 945 780]);
ax(1) = subplot(1,2,1);
image(map)
title([str ' map'])
axis image off
ax(2) = subplot(1,2,2);
image(reshape(colors, [size(colors,1),1,3]))
ylabel(str)
set(ax(2), 'Position', [0.85 0.11 0.05 0.82], 'YAxisLocation', 'right', ...
    'YTick', 1:30:360, 'YTickLabel', 0:15:179, 'XTick', [], 'box', 'off')
set(ax(1), 'Position', [0.13 0.11 0.7 0.82])