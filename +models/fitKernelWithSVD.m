function modelParams = fitKernelWithSVD(neuralData, framesBeforeOnset)
% This function estimates the temporal dynamic (kernel) and response 
% amplitude in response to each stimulus presentation using SVD on the
% matrix containing all responses to all trials of all stimuli [time x
% trials]. Baselines (estimated from mean activity just before stimulus
% onset) are subtracted from all trials.

% neuralData        [time x trials x stimuli]
% framesBeforeOnset int;

% modelParams       .baselineEachTrial: [trials x stimuli]
%                   .kernel: [time x 1]
%                   .alphaEachTrial: [trials x stimuli]

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
modelParams.kernel = kernel;
alphas = reshape(alphas, size(neuralData,2), size(neuralData,3)); % [trials x stimuli]
modelParams.alphaEachTrial = alphas;