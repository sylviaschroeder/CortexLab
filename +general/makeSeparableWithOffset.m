function [alphas, betas, kernel, predictions, residuals] = ...
    makeSeparableWithOffset(data, time, baselineTime)

% data          [trials x time]; each row contains response in a single
%               trial
% time          [1 x time]; time points of trial samples relative to 
%               stimulus onset
% baselineTime  double, in sec; use samples before stimulus onset (between 
%               -baselineTime and 0) to estimate offset for each trial

% alphas        [trials x 1]; multiplicative component for each trial
% betas         [trials x 1]; additive component for each trial
% predictions   [trials x time]; model predictions

baselineInds = time >= -baselineTime & time <= 0;
betas = nanmedian(data(:, baselineInds), 2);

dataRest = bsxfun(@minus, data, betas);
[kernel, alphas, scalar, model, residuals] = MakeSeparable(dataRest, 0);
% make sure that majority of alphas is positive -> temporal kernel reflects
% whether neurons is excited or suppressed by stimulus
if sum(alphas(~isnan(alphas))) < 0
    kernel = -kernel;
    alphas = -alphas;
end
alphas = alphas .* scalar;
predictions = bsxfun(@plus, model, betas);