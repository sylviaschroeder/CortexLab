function [fitKernels, ETA, predictedSignals, predETA, windowTimes] = ...
    kernelRegression(inSignal, t, eventTimes, eventValues, windows, lambda, cvFold)
% function [fitKernels, predictedSignals] = kernelRegression(inSignal, t, eventTimes,
%                                                 eventValues, windows, lambda, cvFold)
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
%
% Some future version of this function could allow uneven sampling of the
% input signal, but this one doesn't. 
% 
% TODO: 
% - for cases in which the events have values, should also fit an
% "intercept" for the same event with values 1
% - Some future version could also allow for fitting as a sum of basis
% functions rather than this simplified "1's" method
% - could hypothetically allow for a kernel that is super-sampled in time
% (would probably need to be using basis functions)
% - should add the ability to fit also a continuous signal, in the same way
% using a toeplitz version of the vector. Would first interpolate it to the
% frame times. 
% - need to take care of "unpredictable point" - out of range of anything -
% here rather than in inputs, since otherwise they artificially inflate cv
% scores (predict zero and get zero a lot). 
      
nSig = size(inSignal,1);

if nargin<6
    lambda = 0; % if this is non-zero, does ridge regression
    cvFold = [0 0]; % number of folds of cross-validation to do
end

% this is the function used to evaluate the cross validated error. Should
% return nSig x 1, the performance on each signal to be predicted.
% cvEvalFunc = @(pred, actual)1-mean(mean((pred-actual).^2))/mean(mean(actual.^2));
% cvEvalFunc = @(pred, actual)1- var(pred-actual); % here assume variance of actual is 1 -- 
                                                 % it is (or is close) if data were zscored. 
                                                 % Otherwise you'd want to divide by it
% cvEvalFunc = @(pred, actual)1- var(pred-actual)./var(actual);

[A, numSamples, windowTimes] = krnl.getToeplitz(t, eventTimes, windows);

% if cvFold(1)>0
%     
%     cvp = cvpartition(nT,'KFold', cvFold(1));
%     cvErr = zeros(nSig,cvFold(2));
%     for k = 1:cvFold(2)
%         fprintf(1, 'cvFold %d/%d\n', k, cvFold(2))
%         if lambda>0
%             % if using regularization, you want the regularization rows to
%             % always be part of training
%             trainInds = vertcat(cvp.training(k), true(size(inSignal,2)-nT,1));
%         else
%             trainInds = cvp.training(k);
%         end
%         
%         testInds = cvp.test(k);
%         
%         trainSetObservations = inSignal(:,trainInds);
%         trainSetPredictors = A(trainInds,:);
%         X = solveLinEq(trainSetPredictors,trainSetObservations'); % X becomes nWinSampsTotal by nS 
%     
%         predictedSignals = (A(testInds,:)*X);
%         testSetObservations = inSignal(:,testInds)';
%         cvErr(:,k) = cvEvalFunc(predictedSignals, testSetObservations);
%     end
%     
% else
% end

X = NaN(sum(numSamples)+1, size(inSignal,2));
for iCell = 1:size(inSignal,2)
%     [B, fitInfo] = lasso(A, inSignal(:,iCell), 'Lambda', lambda);
%     X(:,iCell) = [fitInfo.Intercept; B];
    B = ridge(
end

% X = solveLinEq(Areg,inSignalReg); % X becomes nWinSampsTotal by nS

predictedSignals = [ones(size(A,1),1), A] * X;

fitKernels = cell(1, length(eventTimes));
for ev = 1:length(eventTimes)
    fitKernels{ev} = X(sum(numSamples(1:ev-1))+(1:numSamples(ev)),:);
end

ETA = krnl.getETA(inSignal, A, numSamples);
predETA = krnl.getETA(predictedSignals, A, numSamples);