function [prediction, explainedVar] = predictResponse(trace, traceTimes, ...
    stimFrames, stimTimes, RFtimesInFrames, spatTempRF)

% generate toplitz matrix for stimulus
[stim, time, ~, stimBin] = ...
    whiteNoise.makeStimToeplitz(stimFrames, stimTimes, RFtimesInFrames);

% get neural response
traceBin = median(diff(traceTimes));
numBins = round(stimBin / traceBin);
trace = smoothdata(trace, 1, 'movmean', numBins, 'omitnan');
trace = interp1(traceTimes, trace, time);
% z-score neural response
zTrace = (trace - nanmean(trace,1)) ./ nanstd(trace,0,1);

% delete stim frames for which neuron has NaN
ind = isnan(zTrace);
stim(ind,:) = [];
zTrace(ind,:) = [];

% duplicate stimulus matrix to predict ON part (1st half) and OFF
% part (2nd half)
s = stim;
s(stim < 0) = 0;
stim2 = s;
s = stim;
s(stim > 0) = 0;
stim2 = [stim2, s];
stim2 = (stim2 - nanmean(stim2(:))) ./ nanstd(stim2(:)); % normalise each column of stimulus matrix
clear sdesl

prediction = NaN(size(ind));
prediction(~ind) = stim2 * spatTempRF(:);
explainedVar = 1 - sum((zTrace - prediction).^2, 1) ./ sum(zTrace.^2, 1);