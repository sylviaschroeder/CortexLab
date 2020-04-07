function [modelParams, predictions, nonVisualPerTrial] = ...
    modelResponsePerTrial(neuralData, nonVisualData, framesBeforeOnset)
% This function estimates the additive and mutliplicative impact of a
% non-visual signal (running speed or pupil diameter) on the visually
% driven as well as the baseline responses fo one neuron.
% The following function is fit:
% R_tr(t) = kernel(t) * (alpha_stim*f(nv_tr) + g(nv_tr)) + delta*nv_tr + c
% where:
% R_tr(t):      neural response in trial tr at time t from stimulus onset
% kernel(t):    underlying stimulus response
% alpha_stim:   response strength due to stimulus
% nv_tr:        non-visual signal in trial tr
% f(x):         multiplicative influence of non-visual signal on stimulus
%               response; here f(x) = 1 + beta*x
% g(x):         additive influence of non-visual signal on stimulus
%               response; here g(x) = gamma*x
% delta:        influence of non-visual signal on baseline response
% c:            constant (baseline)
% rewritten:
% R_tr(t) = kernel(t) * (alpha + alpha*beta*nv_tr + gamma*nv_tr) + ...
%           delta*nv_tr + c

% neuralData        [time x trials x stimuli]
% nonVisualData     [time x trials x stimuli]
% framesBeforeOnset int;

% modelParams       .baselineEachTrial: [trials x stimuli]
%                   .kernel: [time x 1]
%                   .alphaEachTrial: [trials x stimuli]
%                   .alphasUncond: [1 x stimuli]
%                   .gammasUncond: [1 x stimuli]
%                   .alphas: [1 x stimuli]
%                   .gammaStim: [1 x stimuli]
%                   .beta: double
%                   .gamma: double
%                   .delta: double
%                   .const: double


% estimate baseline for each trial
baselines = median(neuralData(1:framesBeforeOnset,:,:),1);
modelParams.baselineEachTrial = squeeze(baselines);

% subtract baselines from trial responses
dataSubt = bsxfun(@minus, neuralData, baselines);

% fit kernel and response strength for each trial (on baseline subtracted
% data)
[alphas, kernel, scalar] = MakeSeparable(reshape(dataSubt, ...
    size(neuralData,1), []), 0);
% normalize kernel to maximum of 1
m = max(abs(kernel));
kernel = kernel / m;
alphas = alphas * scalar * m;
% make sure that kernel is positive going
if max(kernel) < abs(min(kernel))
    kernel = -kernel;
    alphas = -alphas;
end
% make sure that majority of alphaTrials is positive -> temporal kernel reflects
% whether neurons is excited or suppressed by stimulus
% if sum(alphas(~isnan(alphas))) < 0
%     kernel = -kernel;
%     alphas = -alphas;
% end
modelParams.kernel = kernel;
alphas = reshape(alphas, size(neuralData,2), size(neuralData,3)); % [trials x stimuli]
modelParams.alphaEachTrial = alphas;

% from all alphas (one for each trial), find the best alpha_stim, beta, and
% gamma
nonVis = squeeze(mean(nonVisualData,1)); % [trials x stimuli]
% (1) fit a separate line to each stimulus
alphaUncond = NaN(1,size(alphas,2));
gammaUncond = NaN(1,size(alphas,2));
for iStim = 1:size(alphas,2)
    x = [nonVis(:,iStim), ones(size(neuralData,2),1)] \ alphas(:,iStim);
    alphaUncond(iStim) = x(2);
    gammaUncond(iStim) = x(1);
end
% (2) fit full model
xx = mat2cell(nonVis, size(nonVis,1), ones(1, size(nonVis,2)));
yy = mat2cell(alphas, size(alphas,1), ones(1, size(alphas,2)));
[coeffs, x0y0] = general.bilinfit(xx, yy);
alphaStim = coeffs(:,2)';
beta = -1/x0y0(1);
gamma = -beta*x0y0(2);
alphaNonV = alphaStim*beta + gamma;



% % (1) Fit alphas = nv_tr * alpha_nonVisual + alpha_stim
% alphaStim = NaN(size(neuralData,3), 1);
% alphaNonV = NaN(size(neuralData,3), 1);
% for iStim = 1:size(neuralData,3)
%     x = [nonVis(:,iStim), ones(size(neuralData,2),1)] \ alphas(:,iStim);
%     alphaStim(iStim) = x(2);
%     alphaNonV(iStim) = x(1);
% end
% % (2) Fit alpha_nonVisual = alpha_stim * beta + gamma
% params = [alphaStim, ones(size(alphaStim))] \ alphaNonV;
% beta = params(1);
% gamma = params(2);
modelParams.alphasUncond = alphaUncond;
modelParams.gammasUncond = gammaUncond;
modelParams.alphas = alphaStim;
modelParams.gammaStim = alphaNonV;
modelParams.beta = beta;
modelParams.gamma = gamma;

% Find delta and c based on baselines
params = [nonVis(:), ones(numel(nonVis),1)] \ baselines(:);
delta = params(1);
const = params(2);
modelParams.delta = params(1);
modelParams.const = params(2);

% calculate model predictions
% (1) calculate response strength for each trial (factor multiplied with
% kernel)
predictions = bsxfun(@plus, alphaStim, bsxfun(@times, alphaStim, ...
    beta .* nonVis)) + gamma .* nonVis;
% (2) get temporal dynamics by multiplying with kernel
predictions = bsxfun(@times, kernel, permute(predictions, [3 1 2]));
% (3) add baseline for each trial
predictions = bsxfun(@plus, predictions, ...
    permute(delta .* nonVis, [3 1 2])) + const;

nonVisualPerTrial = nonVis;