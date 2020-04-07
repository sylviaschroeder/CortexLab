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
% -- eventValues is a cell array of "values" for each event, like if you want
% different instances of the event to be fit with a scaled version of the
% kernel. E.g. contrast of stimulus or velocity of wheel movement.
% -- windows is a cell array of 2 by 1 windows, [startOffset endOffset]
% -- lambda is a scalar, the regularization amount. 0 to do no
% regularization
% -- cvFold is 2 by 1, [foldSize, nToCalculate]. So [5 5] does 5-fold CV
% and calculates the error on all five of the test sets. [5 1] still holds
% out 20% but only does this once. 
% 
% fit_kernels is a cell array of nS by nW fit kernels

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