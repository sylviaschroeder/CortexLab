function plotPositionTuning(calciumTraces, samplingRate, stimulusMatrix, ...
    stimulusSequence)

% Number of cells per figure
cellsPerFigure = 5;
colors = lines(cellsPerFigure);
% Define time points of signal to consider (in sec)
limitsAfterOnset = [1 3];

limitsInFrames = round(limitsAfterOnset * samplingRate);
stimTraces = ssLocal.getTracesPerStimulus(calciumTraces, stimulusMatrix, ...
    limitsInFrames);
tempMeans = nanmean(stimTraces, 4);
meanResponses = nanmean(tempMeans, 3);
semResponses = nanstd(tempMeans, 0, 3) ./ sqrt(sum(~isnan(tempMeans), 3));

% Parse stimulus information
xPositions = zeros(0,2); %1st col: x-pos, 2nd col: ID of stimulus
yPositions = zeros(0,2);
for stim = 1:length(stimulusSequence.labels)
    description = stimulusSequence.labels{stim};
    if strcmp('blank', description)
        continue
    end
%     parameters = sscanf(description, 'ori = %d, xfocus = %d, yfocus = %d');
%     if parameters(1) == 0 % vertical orientation -> xpos
%         xPositions(end+1,:) = [parameters(2), stim];
%     elseif parameters(1) == 90 % horizontal orientation -> ypos
%         yPositions(end+1,:) = [parameters(3), stim];
%     end
    
    parameters = sscanf(description, 'len = %d, xfocus = %d, yfocus = %d');
    if parameters(1) == 750 % vertical orientation -> xpos
        xPositions(end+1,:) = [parameters(2), stim];
    elseif parameters(1) == 2700 % horizontal orientation -> ypos
        yPositions(end+1,:) = [parameters(3), stim];
    end
end
xPositions(:,1) = xPositions(:,1) ./ 10;
yPositions(:,1) = yPositions(:,1) ./ 10;
xPositions = sortrows(xPositions);
yPositions = sortrows(yPositions);

% Plot stimulus-triggered responses
numFigures = ceil(size(calciumTraces, 2) / cellsPerFigure);
rows = cellsPerFigure;
cols = 2;
for fig = 1:numFigures
    figure('Position', [624 50 560 1064])
    % Each row: tunings of one cell
    for neuron = 1:min(cellsPerFigure, size(calciumTraces, 2)-(fig-1)*cellsPerFigure)
        neuronID = (fig-1)*cellsPerFigure + neuron;
        subplot(rows, cols, neuron*cols-1)
        hold on
        fill(xPositions([1:end, end:-1:1],1), ...
            [meanResponses(neuronID,xPositions(:,2)) + ...
            semResponses(neuronID,xPositions(:,2)), ...
            meanResponses(neuronID,xPositions(end:-1:1,2)) - ...
            semResponses(neuronID,xPositions(end:-1:1,2))], ...
            'k', 'EdgeColor', 'none', 'FaceColor', colors(neuron,:), ...
            'FaceAlpha', 0.3)
        plot(xPositions(:,1), meanResponses(neuronID,xPositions(:,2)), ...
            'Color', colors(neuron,:), 'LineWidth', 1)
        
        subplot(rows, cols, neuron*cols)
        hold on
        fill(yPositions([1:end, end:-1:1],1), ...
            [meanResponses(neuronID,yPositions(:,2)) + ...
            semResponses(neuronID,yPositions(:,2)), ...
            meanResponses(neuronID,yPositions(end:-1:1,2)) - ...
            semResponses(neuronID,yPositions(end:-1:1,2))], ...
            'k', 'EdgeColor', 'none', 'FaceColor', colors(neuron,:), ...
            'FaceAlpha', 0.3)
        plot(yPositions(:,1), meanResponses(neuronID,yPositions(:,2)), ...
            'Color', colors(neuron,:), 'LineWidth', 1)
        
        maxi = max(meanResponses(neuronID,[xPositions(:,2);yPositions(:,2)]) + ...
            semResponses(neuronID,[xPositions(:,2);yPositions(:,2)]));
        mini = min(meanResponses(neuronID,[xPositions(:,2);yPositions(:,2)]) - ...
            semResponses(neuronID,[xPositions(:,2);yPositions(:,2)]));
        range = maxi - mini;
        mini = mini - 0.05*range;
        maxi = maxi + 0.05*range;
        for plots = 1:2
            subplot(rows, cols, ((neuron-1)*cols)+plots)
            ylim([mini maxi])
            if neuron == cellsPerFigure || ...
                    neuronID == size(calciumTraces, 2)
                if plots == 1
                    xlabel('X-position (degrees)')
                else
                    xlabel('Y-position (degrees)')
                end
            end
            if neuron == 1
                if plots == 1
                    title('X-position tuning')
                else
                    title('Y-position tuning')
                end
            end
            if plots == 1
                xlim(xPositions([1 end],1))
                ylabel('\DeltaF/F_0')
            else
                xlim(yPositions([1 end],1))
            end
        end
    end
end