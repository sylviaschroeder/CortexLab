function plotTuningOfBrainRegions(dataPatches, samplingRate, stimulusMatrix, ...
    stimulusSequence, plane)

% Define time points of signal to consider (in sec)
timeAfterOnset = 0;
timeAfterOffset = 0;

colors = zeros(size(dataPatches,1), size(dataPatches,2), 3);
colors(:,:,1) = repmat(flip(0:size(dataPatches,1)-1)' ./ ...
    (size(dataPatches,1)-1), 1, size(dataPatches,2)); % red
colors(:,:,2) = repmat((0:size(dataPatches,1)-1)' ./ ...
    (size(dataPatches,1)-1), 1, size(dataPatches,2)); % green
colors(:,:,3) = repmat((0:size(dataPatches,2)-1) ./ ...
    (size(dataPatches,2)-1), size(dataPatches,1), 1); % blue

limitsInFrames = round([timeAfterOnset timeAfterOffset] .* samplingRate);
% arrange all pixels of one frame in one row, and time varations in
% columns
d = reshape(dataPatches, size(dataPatches,1) * size(dataPatches,2), ...
    size(dataPatches,3))';

figHand = gratings.plotOrientationResponses(d, samplingRate, stimulusMatrix, ...
    stimulusSequence);
for f = 1:length(figHand)
    figure(figHand(f))
    if nargin < 5
        annotation('textbox', [0.4 0.93 0.2 0.05], 'String', ...
            sprintf('Subset %d', f), ...
            'HorizontalAlignment', 'center', 'EdgeColor', 'none', ...
            'FontWeight', 'bold');
    else
        annotation('textbox', [0.4 0.93 0.2 0.05], 'String', ...
        sprintf('Plane %d - subset %d', plane, f), ...
        'HorizontalAlignment', 'center', 'EdgeColor', 'none', ...
        'FontWeight', 'bold');
    end
end

stimTraces = ssLocal.getTracesPerStimulus(d, stimulusMatrix, ...
    limitsInFrames); % [patches x stimuli x repetitions x time]
tempMeans = nanmean(stimTraces, 4);

% Parse stimulus information
orientations = zeros(0,2); %1st col: orientation, 2nd col: ID of stimulus
for stim = 1:length(stimulusSequence.labels)
    description = stimulusSequence.labels{stim};
    parameters = sscanf(description, 'ori = %d, cg = %d');
    if parameters(2) ~= 0 % contrast is not 0
        orientations(end+1,:) = [parameters(1), stim];
    else
        blank = stim;
    end
end
orientations = sortrows(orientations);

% subtract response to blank from stimulus evoked responses, then divide by
% it
blankResponse = nanmedian(tempMeans(:,blank,:), 3);
tempMeans_blankSub = bsxfun(@rdivide, ...
    bsxfun(@minus, tempMeans, blankResponse), blankResponse);
medianResponses = nanmedian(tempMeans_blankSub, 3);
medianResponses(medianResponses < 0) = 0;
semResponses = nanstd(tempMeans_blankSub, 0, 3) ./ ...
    sqrt(sum(~isnan(tempMeans_blankSub), 3));
medianResponses = medianResponses(:, [orientations(:,2); blank]);
semResponses = semResponses(:, [orientations(:,2); blank]);

figure('Position', [700 140 1080 960], 'Color', 'w')
% the direction of the plotted vector represents the drift direction of the
% grating (0 degrees -> rightward, 90 degrees -> downward); therefore
% orientation vector has to be adapted (otherwise 90 degrees would point
% upward)
oris = 360 - orientations(:,1)';
max_x = bsxfun(@times, sin(oris ./ 180 .* pi), medianResponses(:,1:end-1) + ...
    semResponses(:,1:end-1));
max_x = max(abs(max_x(:)));
max_y = bsxfun(@times, cos(oris ./ 180 .* pi), medianResponses(:,1:end-1) + ...
    semResponses(:,1:end-1));
max_y = max(abs(max_y(:)));
maxi = max(max_x, max_y);

rimFigure = 0.05;
width = (1 - 2*rimFigure) / size(dataPatches,2);
height = (1 - 2*rimFigure) / size(dataPatches,1);
rimPlot = 0.01;

for col = 1:size(dataPatches,2)
    for row = 1:size(dataPatches,1)
        subplot('Position', [rimFigure+(col-1)*width+rimPlot ...
            rimFigure+(size(dataPatches,1)-row)*height+rimPlot ...
            width-2*rimPlot height-2*rimPlot])
        hold on
        patch = (col-1)*size(dataPatches,1)+row;
        h = polar(oris([1:end 1]) ./ 180 .* pi, ...
            medianResponses(patch, [1:end-1 1]), 'k');
        set(h, 'LineWidth', 2, 'Color', colors(row,col,:))
        h = polar(oris([1:end 1]) ./ 180 .* pi, ...
            medianResponses(patch, [1:end-1 1]) + ...
            semResponses(patch, [1:end-1 1]), 'k');
        set(h, 'LineWidth', 0.5, 'Color', colors(row,col,:))
        r = medianResponses(patch, [1:end-1 1]) - ...
            semResponses(patch, [1:end-1 1]);
        r(r<0) = 0;
        h = polar(oris([1:end 1]) ./ 180 .* pi, r, 'k');
        set(h, 'LineWidth', 0.5, 'Color', colors(row,col,:))
        
        plot([-maxi maxi], [0 0], 'k:')
        plot([0 0], [-maxi maxi], 'k:')
        
        axis([-maxi maxi -maxi maxi])
        axis off
        
%         set(gca, 'XTick', [], 'YTick', [])
        
        if col == 1 && row == 1
            text(0, 1.3*maxi, 'medial/top', ...
                'HorizontalAlignment', 'center', ...
                'FontWeight', 'bold')
            text(-1.3*maxi, 0, 'anterior/central', ...
                'HorizontalAlignment', 'center', 'Rotation', 90, ...
                'FontWeight', 'bold')
            text(-1.1*maxi, 0.1*maxi, sprintf('%.1f',maxi))
            text(-1.1*maxi, -0.1*maxi, '\DeltaF/F')
        elseif col == size(dataPatches,2) && row == 1
            text(0, 1.3*maxi, 'lateral/bottom', ...
                'HorizontalAlignment', 'center', ...
                'FontWeight', 'bold')
        elseif col == 1 && row == size(dataPatches,1)
            text(-1.3*maxi, 0, 'posterior/peripheral', ...
                'HorizontalAlignment', 'center', 'Rotation', 90, ...
                'FontWeight', 'bold')
        end
    end
end
if nargin > 4
    annotation('textbox', [0.4 0.95 0.2 0.05], 'String', ...
        sprintf('Plane %d', plane), ...
        'HorizontalAlignment', 'center', 'EdgeColor', 'none', ...
        'FontWeight', 'bold');
end