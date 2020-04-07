function mapRetinotopy(calciumTraces, frameTimes, stimulusMatrix, stimulusSequence)

% Number of cells per figure
cellsPerFigure = 10;
colors = lines(cellsPerFigure);
% Define length of stimulus-triggered kernels
framesBeforeOnset = 5;
framesAfterOffset = 10;

% Get frames indices for each repetition of all stimuli
onsets = max(0, [stimulusMatrix(:,1), diff(stimulusMatrix, 1, 2)]);
[stimIDs, onsets] = find(onsets);
repetitions = sum(stimIDs == 1);
stimDuration = round(sum(stimulusMatrix(1,:)) / repetitions);
kernelSize = framesBeforeOnset + stimDuration + framesAfterOffset;
stimFrames = NaN(repetitions, kernelSize, size(stimulusMatrix,1));
for stim = 1:size(stimulusMatrix, 1)
    stimFrames(:,:,stim) = repmat(onsets(stimIDs == stim) - ...
        framesBeforeOnset, 1, kernelSize) + ...
        repmat(0:kernelSize-1, repetitions, 1);
end
% If for some stimuli kernel is larger than calcium traces, add NaN frames
% to beginning or end of traces
minFrame = min(stimFrames(:));
if minFrame < 1
    addFrames = -minFrame + 1;
    stimFrames = stimFrames + addFrames;
    calciumTraces = [NaN(addFrames, size(calciumTraces, 2)); calciumTraces];
end
maxFrame = max(stimFrames(:));
if maxFrame > size(calciumTraces, 1)
    addFrames = maxFrame - size(calciumTraces, 1);
    calciumTraces = [calciumTraces; NaN(addFrames, size(calciumTraces, 2))];
end

% Parse stimulus information
xPositions = zeros(0,2); %1st col: x-pos, 2nd col: ID of stimulus
yPositions = zeros(0,2);
for stim = 1:length(stimulusSequence.labels)
    description = stimulusSequence.labels{stim};
    if strcmp('blank', description)
        blank = stim;
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
xPositions = sortrows(xPositions);
yPositions = sortrows(yPositions);

% Plot stimulus-triggered responses
numFigures = ceil(size(calciumTraces, 2) / cellsPerFigure);
rows = cellsPerFigure + 1;
cols = size(stimulusMatrix, 1);
kernelTimes = ((0:kernelSize-1)-framesBeforeOnset) * diff(frameTimes([1 2]));
screen = get(0, 'ScreenSize');
for fig = 1:numFigures
    figure('Position', screen)
    % First row: sketches of stimuli
    for xpos = 1:size(xPositions, 1)
        subplot(rows, cols, xpos);
        stim = ones(1, size(xPositions, 1));
        stim(xpos) = 0;
        imagesc(stim)
        colormap gray
        set(gca, 'XTick', [], 'YTick', [])
    end
    for ypos = 1:size(yPositions, 1)
        subplot(rows, cols, size(xPositions,1)+ypos)
        stim = ones(size(yPositions, 1), 1);
        stim(ypos) = 0;
        imagesc(stim)
        colormap gray
        set(gca, 'XTick', [], 'YTick', [])
    end
    subplot(rows, cols, cols)
    imagesc(0.5, [0 1])
    colormap gray
    set(gca, 'XTick', [], 'YTick', [])
    
    % Each following row: kernels of cells
    for neuron = 1:min(cellsPerFigure, size(calciumTraces, 2)-(fig-1)*cellsPerFigure)
        neuronID = (fig-1)*cellsPerFigure + neuron;
        maxi = -Inf;
        mini = Inf;
        for xpos = 1:size(xPositions,1)
            stimID = xPositions(xpos, 2);
            responses = calciumTraces(stimFrames(:,:,stimID), neuronID);
            responses = reshape(responses, repetitions, []);
            kernel = nanmean(responses, 1);
            sd = nanstd(responses, 0, 1);
            
            subplot(rows, cols, neuron*cols+xpos)
            hold on
            fill(kernelTimes([1:end,end:-1:1]), [kernel+sd, ...
                kernel(end:-1:1)-sd(end:-1:1)], 'k', 'EdgeColor', 'none', ...
                'FaceColor', colors(neuron,:), 'FaceAlpha', 0.2)
            plot(kernelTimes, kernel, 'Color', colors(neuron,:))
            xlim(kernelTimes([1 end]))
            m = max(kernel+sd);
            if maxi < m
                maxi = m;
            end
            m = min(kernel-sd);
            if mini > m
                mini = m;
            end
            if xpos == 1
                ylabel(sprintf('Cell %d', neuronID))
            end
        end
        for ypos = 1:size(yPositions,1)
            stimID = yPositions(ypos, 2);
            responses = calciumTraces(stimFrames(:,:,stimID), neuronID);
            responses = reshape(responses, repetitions, []);
            kernel = nanmean(responses, 1);
            sd = nanstd(responses, 0, 1);
            
            subplot(rows, cols, neuron*cols+size(xPositions,1)+ypos)
            hold on
            fill(kernelTimes([1:end,end:-1:1]), [kernel+sd, ...
                kernel(end:-1:1)-sd(end:-1:1)], 'k', 'EdgeColor', 'none', ...
                'FaceColor', colors(neuron,:), 'FaceAlpha', 0.2)
            plot(kernelTimes, kernel, 'Color', colors(neuron,:))
            xlim(kernelTimes([1 end]))
            m = max(kernel+sd);
            if maxi < m
                maxi = m;
            end
            m = min(kernel-sd);
            if mini > m
                mini = m;
            end
            if xpos == 1
                title(sprintf('Cell %d', neuron))
            end
        end
        stimID = blank;
        responses = calciumTraces(stimFrames(:,:,stimID), neuronID);
        responses = reshape(responses, repetitions, []);
        kernel = nanmean(responses, 1);
        sd = nanstd(responses, 0, 1);
        
        subplot(rows, cols, (neuron+1)*cols)
        hold on
        fill(kernelTimes([1:end,end:-1:1]), [kernel+sd, ...
            kernel(end:-1:1)-sd(end:-1:1)], 'k', 'EdgeColor', 'none', ...
            'FaceColor', colors(neuron,:), 'FaceAlpha', 0.2)
        plot(kernelTimes, kernel, 'Color', colors(neuron,:))
        xlim(kernelTimes([1 end]))
        m = max(kernel+sd);
        if maxi < m
            maxi = m;
        end
        m = min(kernel-sd);
        if mini > m
            mini = m;
        end
        if xpos == 1
            title(sprintf('Cell %d', neuron))
        end
        
        range = maxi - mini;
        if range < 0
            continue
        end
        mini = mini - 0.05*range;
        maxi = maxi + 0.05*range;
        for stim = 1:size(stimulusMatrix,1)
            subplot(rows, cols, (neuron*cols)+stim)
            plot([0 0], [mini maxi], 'k:')
            plot([1 1]*stimDuration*diff(frameTimes([1 2])), [mini maxi], 'k:')
            ylim([mini maxi])
            if neuron < min(cellsPerFigure, size(calciumTraces, 2)-fig*cellsPerFigure)
                set(gca, 'XTick', [])
            end
            if stim > 1
                set(gca, 'YTick', [])
            end
        end
    end
%     tightfig(gcf)
end