function responses = getTracesPerStimulus(calciumTraces, ...
    stimulusMatrix, offsets)
% responses: [cells x stimuli x repetitions x timePoints]

% Define length of stimulus-triggered kernels
framesBeforeOnset = offsets(1);
framesAfterOffset = offsets(2);

% Get frames indices for each repetition of all stimuli
onsets = max(0, [stimulusMatrix(:,1), diff(stimulusMatrix, 1, 2)]);
[stimIDs, onsets] = find(onsets);
repetitions = sum(stimIDs == 1); % assuming number of repetitions is the same for all stimuli
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

% Get response for each neuron to each stimulus and each repetition
responses = NaN(size(calciumTraces, 2), size(stimulusMatrix, 1), ...
    repetitions, kernelSize); % [neurons x stimuli x repetitions x timePoints]
for neuron = 1:size(calciumTraces, 2)
    for stim = 1:size(stimulusMatrix, 1)
        resp = calciumTraces(stimFrames(:,:,stim), neuron);
        resp = reshape(resp, repetitions, []);
        responses(neuron, stim, :, :) = resp;
    end
end