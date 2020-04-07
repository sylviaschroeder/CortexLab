function plotOrientationTuning(calciumTraces, samplingRate, stimulusMatrix, ...
    stimulusSequence, plotTag)

% Define time points of signal to consider (in sec)
timeAfterOnset = 0;
timeAfterOffset = 0;

if strcmp(plotTag, 'multi')
    % Number of cells per figure
    rows = 5;
    cols = 5;
    cellsPerFigure = rows * cols;
end

limitsInFrames = round([timeAfterOnset timeAfterOffset] .* samplingRate);
stimTraces = ssLocal.getTracesPerStimulus(calciumTraces, stimulusMatrix, ...
    limitsInFrames); % [neurons x stimuli x repetitions x time]
tempMeans = nanmean(stimTraces, 4);
medianResponses = nanmedian(tempMeans, 3);
semResponses = nanstd(tempMeans, 0, 3) ./ sqrt(sum(~isnan(tempMeans), 3));

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

[neuronOrientations, neuronDirections, ~, scaledResponsesOri, ...
    scaledResponsesDir, oris, dirs] = ...
    ssLocal.getPreferredOrientations(calciumTraces, samplingRate, ...
    stimulusMatrix, stimulusSequence, timeAfterOnset, timeAfterOffset);
if oris(1) == 0
    oris(end+1) = 180;
    scaledResponsesOri(:,end+1) = scaledResponsesOri(:,1);
end
dirs(end+1) = dirs(1);
scaledResponsesDir(:,end+1) = scaledResponsesDir(:,1);

for neuron = 1:size(calciumTraces,2)
    if strcmp(plotTag, 'multi')
        if mod(neuron, cellsPerFigure) == 1
            figure('Position', [440 260 1060 860])
            
        end
        sp = mod(neuron, cellsPerFigure);
        if sp == 0
            sp = cellsPerFigure;
        end
        subplot(rows, cols, sp)
    else
        figure('Position', [350 680 1560 420])
        subplot(1,3,1)
    end
    plot(orientations(:,1), permute(tempMeans(neuron, orientations(:,2), :), ...
        [3 2 1]), 'k.', 'MarkerSize', 5)
    hold on
    plot(orientations(:,1), medianResponses(neuron,orientations(:,2)), ...
        'k', 'LineWidth', 2)
    % plot response to blank stimulus
    plot(orientations([1 end],1), [1 1] * medianResponses(neuron, blank), 'k:')
    xlim(orientations([1 end],1))
    if strcmp(plotTag, 'single')
        xTicks = 0:45:orientations(end,1);
    else
        xTicks = 0:90:orientations(end,1);
    end
    set(gca, 'XTick', xTicks)
    if strcmp(plotTag, 'single') || neuron > (rows-1)*cols || ...
            neuron > size(calciumTraces, 2)-cols
        xlabel('Direction (in degrees)')
    end
    if strcmp(plotTag, 'single') || neuron == round(cols/2)
        title('Direction tuning')
    end
    if strcmp(plotTag, 'single') || mod(neuron,cols) == 1
        ylabel('Response (F)')
    end
    
    if strcmp(plotTag, 'single')
        subplot(1,3,2)
        h = polar(oris ./ 180 .* pi, scaledResponsesOri(neuron,:)', 'k');
        set(h, 'LineWidth', 3);
        hold on
        theta = neuronOrientations(neuron, 1) / 180 * pi;
        h = compass(cos(theta) * neuronOrientations(neuron,2), sin(theta) * ...
            neuronOrientations(neuron,2), 'k');
        set(h, 'LineWidth', 6)
        title(sprintf('Orientation tuning (strength: %.2f)', ...
            neuronOrientations(neuron,2)))
        
        subplot(1,3,3)
        h = polar([dirs; dirs(1)] ./ 180 .* pi, scaledResponsesDir(neuron,[1:end,1])', 'k');
        set(h, 'LineWidth', 3);
        hold on
        theta = neuronDirections(neuron, 1) / 180 * pi;
        h = compass(cos(theta) * neuronDirections(neuron,2), sin(theta) * ...
            neuronDirections(neuron,2), 'k');
        set(h, 'LineWidth', 6)
        title(sprintf('Direction tuning (strength: %.2f)', ...
            neuronDirections(neuron,2)))
    end
end