function plotRFsOfBrainRegions(RFs, stimPosition)

colors = zeros(size(RFs,1), size(RFs,2), 3);
colors(:,:,1) = repmat(fliplr(0:size(RFs,1)-1)' ./ (size(RFs,1)-1), 1, size(RFs,2)); % red
colors(:,:,2) = repmat((0:size(RFs,1)-1)' ./ (size(RFs,1)-1), 1, size(RFs,2)); % green
colors(:,:,3) = repmat((0:size(RFs,2)-1) ./ (size(RFs,2)-1), size(RFs,1), 1); % blue

%% Plot contour of RF for each patch
RFmatrix = cat(3, RFs{:});
contourLevel = max(RFmatrix(:)) / 2;
xSpace = diff(stimPosition(1:2)) / size(RFs{1,1},2);
xVector = (0.5:size(RFs{1,1},2)) * xSpace + stimPosition(1);
ySpace = diff(stimPosition(3:4)) / size(RFs{1,1},1);
yVector = (0.5:size(RFs{1,1},1)) * ySpace + stimPosition(3);

figure('Position', [680 680 1110 420])

subplot(1,2,1)
hold on
for x = 1:size(RFs,2)
    for y = 1:size(RFs,1)
        contour(xVector, yVector, RFs{y,x}, [1 1] * contourLevel, 'Color', ...
            squeeze(colors(y,x,:)), 'LineWidth', 2);
    end
end
set(gca, 'YDir', 'reverse')
xlabel('Horizontal visual space (in degrees)')
ylabel('Vertical visual space (in degrees)')
title('RF outlines of brain regions')

subplot(1,2,2)
image(colors)
xlabel('Brain space (medial-lateral)')
ylabel('Brain space (anterior-posterior)')
title('Color code of brain regions')