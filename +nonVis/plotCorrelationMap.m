function handle = plotCorrelationMap(rhos, ROImaps, mapSize)

% rhos          [neurons x 1]
% ROImaps       {neurons x 1}; each entry contains pixels of ROI
% mapsSize      [1 x 2]; size of frame ([validY validX])

stepSize = .05;
maxRho = ceil(max(abs(rhos)) / stepSize) * stepSize;
rhoValues = -maxRho:stepSize:maxRho;
labeled = [flip((length(rhoValues)+1)/2:-4:1), ...
    (length(rhoValues)+1)/2+4:4:length(rhoValues)];

map = zeros([mapSize,3]);
colors = jet(length(rhoValues));
for roi = 1:size(rhos,1)
    [~,rhoValue] = min(abs(rhoValues - rhos(roi)));
    col = colors(rhoValue,:);
    map(ROImaps{roi}) = col(1);
    map(prod(mapSize)+ROImaps{roi}) = col(2);
    map(2*prod(mapSize)+ROImaps{roi}) = col(3);
end

figure('Position', [720 275 945 780]);
ax(1) = subplot(1,2,1);
image(map)
axis image off
ax(2) = subplot(1,2,2);
image(reshape(colors, [size(colors,1),1,3]))
ylabel('Rho')
set(ax(2), 'Position', [0.85 0.11 0.05 0.82], 'YAxisLocation', 'right', ...
    'YTick', labeled, 'YTickLabel', rhoValues(labeled), ...
    'XTick', [], 'box', 'off', 'YDir', 'normal')
set(ax(1), 'Position', [0.13 0.11 0.7 0.82])
handle = ax(1);