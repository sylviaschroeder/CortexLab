function [timeCourses, tunings, offsets, baselines, filters, predictions, ...
    errors, traces_baselineSub, kernelTimes] = ...
    getConditionedOrientationTuning_separated(calciumTraces, traceTimes, ...
    conditionValues, conditionTimes, stimulusMatrix, stimulusSequence, lambda)

% calciumTraces     [t x neurons]; calcium traces of neurons
% traceTimes        [t x 1]; time points of calcium trace samples
% conditionValues   [t2 x numConditions]; each column specifies whether a
%                   specific condition is met or not (all entries are 0 or
%                   1); there should be at least 2 conditions, e.g. running
%                   and not-running (depending on the classification
%                   criteria the condition signals are not necessarily
%                   opposites of each other)
% conditionTimes    [t2 x 1]; time points of condition signal samples
% stimulusMatrix    [numStim x t]; entries specify whether stimulus was 
%                   presented; time samples match those of calcium traces
% stimulusSequence  .labels: {numStim x 1}, contains string describing stimulus
%                   .seq: [numTrials x 1], contains stimulus ID presented
%                   in each trial

% for all output variables holds: the order of conditions is: 1st:
% condition 1 as given by 1st column in conditionValues, 2nd: condition 2
% as given by 2nd column in conditionValues, 3rd: all data is used (not
% divided into conditions)
% timeCourses       [kernelTime x neuron x condition], contains 
%                   stimulus-time separated kernel for each neuron and
%                   condition
% tunings           [stimuli x neuron x condition], contains 
%                   stimulus dependent scalings of kernel for each neuron
%                   and condition
% offsets           [1 x neuron x condition]; contains offsets from zero of
%                   scaled time courses for each condition
% baselines         [1 x neuron]; contains mean response before 1st
%                   stimulus for each neuron
% filters           each field: [stimuli x kernelTimes x neurons], contains
%                   kernel representing response of each neuron to each
%                   stimulus (for conditioned filters, number of stimuli is
%                   doubled, one set for each condition)
%                   fields: .cond_regr (conditioned from ridge regression)
%                           .cond_sep (conditioned and stim-time separated)
%                           .nonCond_regr (unconditioned from ridge regr.)
%                           .nonCond_sep (uncond. and stim-time sep.)
% predictions       each field: [t x neurons]; responses predicted by
%                   kernels
%                   fields: .cond_regr (conditioned from ridge regression)
%                           .cond_sep (conditioned and stim-time separated)
%                           .nonCond_regr (unconditioned from ridge regr.)
%                           .nonCond_sep (uncond. and stim-time sep.)
% error             each field: [1 x neurons]; contains measured
%                   goodness of fit of the predicted to the recorded
%                   response
%                   fields: .cond_regr (conditioned from ridge regression)
%                           .cond_sep (conditioned and stim-time separated)
%                           .nonCond_regr (unconditioned from ridge regr.)
%                           .nonCond_sep (uncond. and stim-time sep.)
% traces_baselineSub [t x neurons]; same data as given in calciumTraces but
%                    baseline response (mean activity before presentation of 1st
%                    stimulus) is subtracted
% kernelTimes       [1 x kernelTime]; contains times of kernel relative to
%                   stimulus onset

%% Parameters
% Fraction of trial that has to be within the condition so that whole trial
% is classified to be in that condition
criterium = 0.7;
% Considered time before and after trial (in sec) for categorisation
trialLimits = [0 0];
% Time limits for stimulus triggered response filters (in s)
tMin = -2.5;
tMax = 7; %12;

if nargin < 7
    lambda = [];
end
%% Sort stimulus trials into two conditions
samplingRate = 1 / median(diff(traceTimes));
limitsInFrames = round(trialLimits * samplingRate);

orientations = gratings.getOrientations(stimulusSequence);

stimulusMatrix = stimulusMatrix(orientations(:,2),:); % only consider gratings (not blank)
stM = [zeros(1, size(stimulusMatrix,1)); max(0, diff(stimulusMatrix'))]; % only stimulus onsets

stM_nonCond = stM; % trials are NOT divided by conditions

% determine which trials fall into which condition
condValuesResampled = interp1(conditionTimes, conditionValues, ...
    traceTimes, 'linear', NaN);
condValuesResampled = condValuesResampled >= 0.5;
stimConditions = ssLocal.getTracesPerStimulus(condValuesResampled, ...
    stimulusMatrix, limitsInFrames); % [conditions x stimuli x repetitions x timePoints]
trialConditions = (sum(stimConditions, 4) ./ size(stimConditions, 4)) > criterium;

% double stimuli: one set for each condition
stM_cond = {zeros(size(stM)), zeros(size(stM))};
for stim = 1:size(orientations,1)
    trialTimes = find(stM(:,stim) == 1);
    for cond = 1:2
        condTrials = trialConditions(cond,stim,:) == 1;
        stM_cond{cond}(trialTimes(condTrials),stim) = 1;
    end
end
stM = cat(2, stM_cond{:});

%% Subtract baseline response from calcium traces
[rows, ~] = find(stM);
preStim = 1 : min(rows)-1; % initial blank
baselines = mean(calciumTraces(preStim,:),1);
traces_baselineSub = bsxfun(@minus, calciumTraces, baselines);

% set responses to unclassied trials to zero
nonClassTrials = squeeze(all(trialConditions == 0, 1)); % [stimulus x repetition]
for stim = 1:size(nonClassTrials)
    ind = find(nonClassTrials(stim,:));
    if isempty(ind)
        continue
    end
    for trial = ind
        start = find(stM_nonCond(:,stim) == 1, trial);
        start = start(end);
        finish = start + find(stimulusMatrix(stim,start:end) == 0, 1) - 1;
        traces_baselineSub(start:finish,:) = 0;
    end
end

%% find the filters
samplesMin = floor(tMin * samplingRate);
samplesMax = ceil(tMax * samplingRate);
% list of temporal shifts
kernelTimes = samplesMin : samplesMax;

% 1st: for stimuli divided according to condition
% stimulus matrix with all shifts
stM_shifted = zeros(size(stM,1), length(kernelTimes) * size(stM,2));
for iShift = 1:length(kernelTimes)
    stM_shifted(:, (iShift-1)*size(stM,2) + (1:size(stM,2))) = ...
        circshift(stM, kernelTimes(iShift));
end
% compute filters
if isempty(lambda)
    lambda = general.searchOptimalLambda(stM_shifted, traces_baselineSub, 1);
end
% lambda = 0.5; % equivalent to a half-contrast stimulus shown once
X = vertcat(stM_shifted, lambda * eye(length(kernelTimes) * size(stM,2)));
Y = vertcat(traces_baselineSub, zeros(length(kernelTimes) * size(stM,2), size(traces_baselineSub,2)));
filts = X \ Y;
% responses predicted from filters
predictions.cond_regr = stM_shifted * filts;
% goodness of fit
m = mean(traces_baselineSub);
s = std(traces_baselineSub);
tz = bsxfun(@rdivide, bsxfun(@minus, traces_baselineSub, m), s);
pz = bsxfun(@rdivide, bsxfun(@minus, predictions.cond_regr, m), s);
errors.cond_regr = mean((tz - pz) .^ 2, 1);
% rearrange filters: [stimuli x shifts x cells]
filters.cond_regr = reshape(filts, size(stM,2), length(kernelTimes), size(traces_baselineSub,2));
% set to zero the filters for stimuli that don't exist
nonExistStim = all(stM == 0, 1);
filters.cond_regr(nonExistStim,:,:) = 0;

% 2nd: for stimuli not divided by condition
% stimulus matrix with all shifts
stM_nonCond_shifted = zeros(size(stM_nonCond,1), length(kernelTimes) * size(stM_nonCond,2));
for iShift = 1:length(kernelTimes)
    stM_nonCond_shifted(:, (iShift-1)*size(stM_nonCond,2) + (1:size(stM_nonCond,2))) = ...
        circshift(stM_nonCond, kernelTimes(iShift));
end
% compute filters
if isempty(lambda)
    lambda = general.searchOptimalLambda(stM_nonCond_shifted, traces_baselineSub, 1);
end
% lambda = 0.5; % equivalent to a half-contrast stimulus shown once
X = vertcat(stM_nonCond_shifted, lambda * eye(length(kernelTimes) * size(stM_nonCond,2)));
Y = vertcat(traces_baselineSub, zeros(length(kernelTimes) * size(stM_nonCond,2), ...
    size(traces_baselineSub,2)));
filts = X \ Y;
% responses predicted from filters
predictions.nonCond_regr = stM_nonCond_shifted * filts;
% goodness of fit
pz = bsxfun(@rdivide, bsxfun(@minus, predictions.nonCond_regr, m), s);
errors.nonCond_regr = mean((tz - pz) .^ 2, 1);
% rearrange filters: [stimuli x shifts x cells]
filters.nonCond_regr = reshape(filts, size(stM_nonCond,2), length(kernelTimes), ...
    size(traces_baselineSub,2));
% set to zero the filters for stimuli that don't exist
nonExistStim = all(stM_nonCond == 0, 1);
filters.nonCond_regr(nonExistStim,:,:) = 0;

%% Separate time courses from tuning curves
doGraphics = 0; % change to 1 to see MakeSeparable results
timeCourses = zeros(size(filters.cond_regr,2), size(calciumTraces,2), 3);
tunings = zeros(size(orientations,1), size(calciumTraces,2), 3);
offsets = zeros(1, size(calciumTraces,2), 3);
for iCell = 1:size(traces_baselineSub,2)
    % approximate filters for conditioned responses with separable functions
    [timeCourse, tuning, bestScl, ~, ~, ~, const] = ...
        MakeSeparable(filters.cond_regr(:,:,iCell), doGraphics, 1);
    [~, maxInd] = max(abs(timeCourse));
    scaling = timeCourse(maxInd);
    timeCourses(:,iCell,1) = timeCourse / scaling;
    timeCourses(:,iCell,2) = timeCourse / scaling; % assume same time course for both conditions
    tun = tuning * bestScl * scaling;
    tun(nonExistStim) = NaN; % no data for those
    tunings(:,iCell,1) = tun(1:size(orientations,1));
    tunings(:,iCell,2) = tun(size(orientations,1)+1:end);
    offsets(1,iCell,1) = const;
    offsets(1,iCell,2) = const; % assume same offset for both conditions
    
    % approximate filters for non-conditioned responses with separable functions
    [timeCourse, tuning, bestScl, ~, ~, ~, const] = ...
        MakeSeparable(filters.nonCond_regr(:,:,iCell), doGraphics, 1);
    [~, maxInd] = max(abs(timeCourse));
    scaling = timeCourse(maxInd);
    timeCourses(:,iCell,3) = timeCourse / scaling;
    tun = tuning * bestScl * scaling;
    tun(nonExistStim) = NaN; % no data for those
    tunings(:,iCell,3) = tun;
    offsets(1,iCell,3) = const;
end

%% Predict response based on stim-time separated filters
% 1st: for stimuli divided according to condition
filters.cond_sep = zeros(size(filters.cond_regr));
for iCell = 1:size(traces_baselineSub,2)
    tun = squeeze(tunings(:,iCell,1:2));
    tun(isnan(tun)) = 0;
    filters.cond_sep(1:size(orientations,1),:,iCell) = tun(:,1) * ...
        timeCourses(:,iCell,1)' + offsets(1,iCell,1);
    filters.cond_sep(size(orientations,1)+1:end,:,iCell) = tun(:,2) * ...
        timeCourses(:,iCell,2)' + offsets(1,iCell,2);
end
filts = reshape(filters.cond_sep, size(stM,2) * length(kernelTimes), ...
    size(traces_baselineSub,2));
predictions.cond_sep = stM_shifted * filts;
% goodness of fit
pz = bsxfun(@rdivide, bsxfun(@minus, predictions.cond_sep, m), s);
errors.cond_sep = mean((tz - pz) .^ 2, 1);

% 2nd: for stimuli not divided according to condition
filters.nonCond_sep = zeros(size(filters.nonCond_regr));
for iCell = 1:size(traces_baselineSub,2)
    tun = tunings(:,iCell,3);
    tun(isnan(tun)) = 0;
    filters.nonCond_sep(:,:,iCell) = tun * ...
        timeCourses(:,iCell,3)' + offsets(1,iCell,3);
end
filts = reshape(filters.nonCond_sep, size(stM_nonCond,2) * length(kernelTimes), ...
    size(traces_baselineSub,2));
predictions.nonCond_sep = stM_nonCond_shifted * filts;
% goodness of fit
pz = bsxfun(@rdivide, bsxfun(@minus, predictions.nonCond_sep, m), s);
errors.nonCond_sep = mean((tz - pz) .^ 2, 1);