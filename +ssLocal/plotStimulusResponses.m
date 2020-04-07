function plotStimulusResponses(calciumTraces, samplingRate, stimulusMatrix, ...
    stimulusSequence, framesBeforeOnset, framesAfterOffset)

% Number of cells per figure
cellsPerFigure = 10;
colors = jet(ceil(size(calciumTraces, 2) * 1.1));
colors = colors(end-size(calciumTraces, 2)+1 : end, :);

stimTraces = ssLocal.getTracesPerStimulus(calciumTraces, stimulusMatrix, ...
    [framesBeforeOnset, framesAfterOffset]);

% Parse stimulus information
% labels = strsplit(stimulusSequence.labels{1}, ', ');
% parameters = zeros(size(stimulusMatrix, 1), length(labels));
% for stim = 1:size(stimulusMatrix, 1)
%     parameters(stim,:) = sscanf(stimulusSequence.labels{stim}, '%*s = %d%*c')';
% end

% Plot stimulus-triggered responses
numFigures = ceil(size(calciumTraces, 2) / cellsPerFigure);
rows = cellsPerFigure;
cols = size(stimulusMatrix, 1);
kernelSize = size(stimTraces, 4);
kernelTimes = ((0:kernelSize-1)-framesBeforeOnset) / samplingRate;
kernelTimesPatch  = [kernelTimes, fliplr(kernelTimes)]';
screen = get(0, 'ScreenSize');
stimDuration = size(stimTraces,4) - framesBeforeOnset - framesAfterOffset;
for fig = 1:numFigures
    figure('Position', screen)    
    % Each row: kernels of one cell
    for neuron = 1:min(cellsPerFigure, size(calciumTraces, 2)-(fig-1)*cellsPerFigure)
        neuronID = (fig-1)*cellsPerFigure + neuron;
        maxi = -Inf;
        mini = Inf;
        for stim = 1:size(stimulusMatrix,1)
            resp = squeeze(stimTraces(neuronID, stim, :, :));
            kernel = nanmedian(resp, 1);
            
            subplot(rows, cols, (neuron-1)*cols+stim)
            hold on
            patch(repmat(kernelTimesPatch, 1, size(resp, 1)), ...
                [resp, fliplr(resp)]', 'k', 'EdgeColor', colors(neuronID,:), ...
                'EdgeAlpha', 0.3, 'FaceColor', 'none');

            plot(kernelTimes, kernel, 'Color', colors(neuronID,:), 'LineWidth', 2)
            xlim(kernelTimes([1 end]))
            m = max(resp(:));
            if maxi < m
                maxi = m;
            end
            m = min(resp(:));
            if mini > m
                mini = m;
            end
            if neuron == 1
                title(regexprep(stimulusSequence.labels{stim}, ', ', '\n'))
            end
            if stim == 1
                ylabel(sprintf('Cell %d', neuronID))
            end
            set(gca, 'Color', 'k')
        end
        
        range = maxi - mini;
        mini = mini - 0.01*range;
        maxi = maxi + 0.01*range;
        for stim = 1:size(stimulusMatrix,1)
            subplot(rows, cols, ((neuron-1)*cols)+stim)
            plot([0 0], [mini maxi], 'w:')
            plot([1 1]*stimDuration/samplingRate, [mini maxi], 'w:')
            ylim([mini maxi])
            if neuron < min(cellsPerFigure, size(calciumTraces, 2)-(fig-1)*cellsPerFigure)
                set(gca, 'XTick', [])
            else
                xlabel('Time (s)')
            end
            if stim > 1
                set(gca, 'YTick', [])
            end
        end
    end
end