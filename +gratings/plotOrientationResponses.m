function figHandles = plotOrientationResponses(calciumTraces, samplingRate, ...
    stimulusMatrix, stimulusSequence)

% Number of cells per figure
cellsPerFigure = 10;
% colors = jet(ceil(size(calciumTraces, 2) * 1.1));
% colors = colors(end-size(calciumTraces, 2)+1 : end, :);
% Define length of stimulus-triggered kernels
framesBeforeOnset = 5;
framesAfterOffset = 30;

stimTraces = ssLocal.getTracesPerStimulus(calciumTraces, stimulusMatrix, ...
    [framesBeforeOnset, framesAfterOffset]);

% Parse stimulus information
orientations = zeros(0,2); %1st col: orientation, 2nd col: ID of stimulus
for stim = 1:length(stimulusSequence.labels)
    description = stimulusSequence.labels{stim};
    parameters = sscanf(description, 'ori = %d, cg = %d');
    if parameters(2) == 0 % contrast is 0 -> blank
        blank = stim;
    else
        orientations(end+1,:) = [parameters(1), stim];
    end
end
[~, stimulusOrder] = sort(orientations(:,1));

% Plot stimulus-triggered responses
numFigures = ceil(size(calciumTraces, 2) / cellsPerFigure);
figHandles = zeros(1, numFigures);
rows = cellsPerFigure;
cols = size(stimulusMatrix, 1);
kernelSize = size(stimTraces, 4);
kernelTimes = ((0:kernelSize-1)-framesBeforeOnset) / samplingRate;
kernelTimesPatch  = [kernelTimes, fliplr(kernelTimes)]';
screen = get(0, 'ScreenSize');
stimDuration = size(stimTraces,4) - framesBeforeOnset - framesAfterOffset;
for fig = 1:numFigures
    figHandles(fig) = figure('Position', ...
        [screen(1)+30 screen(2)+50 screen(3)-60 screen(4)-150]);
    % Each row: kernels of one cell
    for neuron = 1:min(cellsPerFigure, size(calciumTraces, 2)-(fig-1)*cellsPerFigure)
        neuronID = (fig-1)*cellsPerFigure + neuron;
        maxi = -Inf;
        mini = Inf;
        for ori = 1:size(orientations,1)
            stimID = orientations(stimulusOrder(ori), 2);
            resp = squeeze(stimTraces(neuronID, stimID, :, :));
            kernel = nanmedian(resp, 1);
            sd = nanstd(resp, 0, 1);
            
            subplot(rows, cols, (neuron-1)*cols+ori)
            hold on
%             fill(kernelTimes([1:end,end:-1:1]), [kernel+sd, ...
%                 kernel(end:-1:1)-sd(end:-1:1)], 'k', 'EdgeColor', 'none', ...
%                 'FaceColor', colors(neuron,:), 'FaceAlpha', 0.2)

            patch(repmat(kernelTimesPatch, 1, size(resp, 1)), ...
                [resp, fliplr(resp)]', 'k', 'EdgeAlpha', 0.3, 'FaceColor', 'none');

            plot(kernelTimes, kernel, 'k', 'LineWidth', 2)
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
                title([num2str(orientations(stimID,1)) ' degrees'])
            end
            if ori == 1
                ylabel(sprintf('Cell %d', neuronID))
            end
%             set(gca, 'Color', 'k')
        end
        
        stimID = blank;
        resp = squeeze(stimTraces(neuronID, stimID, :, :));
        kernel = nanmedian(resp, 1);
        sd = nanstd(resp, 0, 1);
        subplot(rows, cols, neuron*cols)
        hold on
%         fill(kernelTimes([1:end,end:-1:1]), [kernel+sd, ...
%             kernel(end:-1:1)-sd(end:-1:1)], 'k', 'EdgeColor', 'none', ...
%             'FaceColor', colors(neuron,:), 'FaceAlpha', 0.2)
        patch(repmat(kernelTimesPatch, 1, size(resp, 1)), ...
                [resp, fliplr(resp)]', 'k', 'EdgeAlpha', 0.3, 'FaceColor', 'none');
        
        plot(kernelTimes, kernel, 'k', 'LineWidth', 2)
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
            title('Blank')
        end
%         set(gca, 'Color', 'k')
        
        range = maxi - mini;
        mini = mini - 0.05*range;
        maxi = maxi + 0.05*range;
        for stim = 1:size(stimulusMatrix,1)
            subplot(rows, cols, ((neuron-1)*cols)+stim)
            plot([0 0], [mini maxi], 'k:')
            plot([1 1]*stimDuration/samplingRate, [mini maxi], 'k:')
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
%     tightfig(gcf)
end