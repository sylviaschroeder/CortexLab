function [tuningCurves, orientations, stimTraces, trialConditions] = ...
    getConditionedOrientationTuning_raw(calciumTraces, traceTimes, ...
    conditionValues, stimulusMatrix, stimulusSequence)

% calciumTraces     [t x neurons]; calcium traces of neurons
% traceTimes        [t x 1]; time points of calcium trace samples
% conditionValues   [t x 1]; each entry specifies which condition is met at
%                   the time point; 0 ... no condition is met, 1,2,...
%                   specifies class
% stimulusMatrix    [numStim x t]; entries specify whether stimulus was 
%                   presented; time samples match those of calcium traces
% stimulusSequence  .labels: {numStim x 1}, contains string describing stimulus
%                   .seq: [numTrials x 1], contains stimulus ID presented
%                   in each trial

% tuningCurves      [neuron x stimulus x condition]; contains median
%                   response of each neuron to each stimulus in each
%                   condition; last condition: all data not divided into
%                   conditions; last stimulus is blank
% stimTraces        [neuron x stimulus x repetition x time]; contains response 
%                   averaged across each trial
% trialConditions   {stimulus x repetition}; states for each
%                   trial which condition it belongs to

% Fraction of trial that has to be within the condition so that whole trial
% is classified to be in that condition
criterium = 0; % if 0, trial falls into class with most occurrences
% Considered time before and after trial (in sec)
trialLimits = [0 0];

samplingRate = 1 / median(diff(traceTimes));
limitsInFrames = round(trialLimits * samplingRate);
stimTraces = ssLocal.getTracesPerStimulus(calciumTraces, stimulusMatrix, ...
    limitsInFrames); % [neurons x stimuli x repetitions x timePoints]
timeAverageResp = nanmean(stimTraces, 4);
% [neurons x stimulus x repetitions]

[orientations, blank] = gratings.getOrientations(stimulusSequence);
timeAverageResp = [timeAverageResp(:, orientations(:,2), :), ...
    timeAverageResp(:, blank, :)];

stimConditions = ssLocal.getTracesPerStimulus(conditionValues', ...
    stimulusMatrix, limitsInFrames); % [1 x stimuli x repetitions x timePoints]
stimConditions = [stimConditions(:, orientations(:,2), :, :), ...
    stimConditions(:, blank, :, :)];
[~, frequency, trialConditions] = mode(squeeze(stimConditions), 3);
% [stimuli x repetitions]

tuningCurves = NaN(size(calciumTraces,2), size(orientations,1)+1, max(conditionValues)+1);
tuningCurves(:,:,end) = nanmedian(timeAverageResp, 3);
for cond = 1:max(conditionValues)
    resp = timeAverageResp;
    ind = frequency./size(stimTraces,4) >= criterium & ...
        cellfun(@ismember, num2cell(ones(size(trialConditions))*cond), trialConditions);
    trialConditions(frequency./size(stimTraces,4) < criterium) = {0};
    resp(:,~ind) = NaN;
    tuningCurves(:,:,cond) = nanmedian(resp, 3);
end