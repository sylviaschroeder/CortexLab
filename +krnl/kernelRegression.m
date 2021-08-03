function [fitKernels, ETA, predictedSignals, predETA, windowTimes, ...
    alignedTraces] = ...
    kernelRegression(inSignal, t, eventTimes, eventWindows, ...
    vectors, vectorWindows, lambda, normalise)
%
% Fits the "toeplitz regression" from Kenneth. 
%
% -- inSignal is nS by nTimePoints, any number of signals to be fit
% -- t is 1 by nTimePoints
% -- eventTimes is a cell array of times of each event, where each vector
% is a different event
% -- eventWindows is a cell array of 2 by 1 windows, [startOffset endOffset]
% -- vectors is a cell array of continuous signals used as predictors (e.g.
% running speed). Each vector in the cell is a different predictor.
% -- vectorWindows is a cell array of 2 by 1 windows, [startOffset endOffset]
% -- lambda is a scalar, the regularization amount. 0 to do no
% regularization
% -- normalise is true or false (default); if true, STD of each predictor,
% i.e. each column of toeplitz matrix is set to 1
% 
% fitKernels is a cell array of fit kernels
% ETA is a cell array of event triggered averages (one per event type)
% predictedSignals is the prediction using the kernels (same size as
% inSignal)
% predETA is a cell array of event triggered averages using the
% predictedSignals
% windowTimes are the kernel times relative to the event
% alignedTraces are the traces aligned to even times

% To use ridge regression instead of Lasso, comment out lines 44+45 and
% remove comment of line 47.

if nargin<7
    lambda = 0; % if this is non-zero, does ridge regression
end
if nargin<8
    normalise = false;
end

[A, numSamples, windowTimes] = krnl.getToeplitz(t, eventTimes, ...
    eventWindows, vectors, vectorWindows, normalise);

X = NaN(sum(numSamples)+1, size(inSignal,2));
for iCell = 1:size(inSignal,2)
    [B, fitInfo] = lasso(A, inSignal(:,iCell), 'Lambda', lambda);
    X(:,iCell) = [fitInfo.Intercept; B];
%     X(:,iCell) = ridge(inSignal(:,iCell), A, lambda, 0);
end

predictedSignals = [ones(size(A,1),1), A] * X;

fitKernels = cell(1, length(eventTimes)+length(vectors));
for ev = 1:length(eventTimes)
    fitKernels{ev} = X(1+sum(numSamples(1:ev-1))+(1:numSamples(ev)),:);
end
for v = 1:length(vectors)
    k = length(eventTimes)+v;
    fitKernels{k} = X(1+sum(numSamples(1:k-1))+(1:numSamples(k)),:);
end

if nargout > 5
    [ETA, alignedTraces] = krnl.getETA(inSignal, A, numSamples);
else
    ETA = krnl.getETA(inSignal, A, numSamples);
end
predETA = krnl.getETA(predictedSignals, A, numSamples);