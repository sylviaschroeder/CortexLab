function [receptiveFields, explainedVariance, predictions, time] = ...
    getReceptiveField(traces, traceTimes, ...
    stimFrames, stimTimes, RFtimesInFrames, ...
    lambdas, crossFolds)

%GETRECEPTIVEFIELD Returns spatiotemporal receptive field.
%   [receptiveFields, explainedVariance, predictions, time] = ...
%    GETRECEPTIVEFIELD(traces, traceTimes, ...
%    stimFrames, stimTimes, RFtimesInFrames, ...
%    lambdas, crossFolds) calculates the linear RF of the neuron.
%
%   receptiveFields     [rows x cols x RFframes x RFtype x neuron]
%                       containing linear regression solution for x in Ax=B
%                       where A is stimulus [rows x cols x time] and B is
%                       calcium response, for each neuron and stimulus 
%                       model; ridge regression is performed on all data
%                       using the optimal lambda value found with
%                       cross-validation
%   explainedVariance   [neuron x lambdaStim x crossFold], each entry:
%                       explained variance for fitted RF for
%                       each neuron, lambda, and cross val. fold
%   predictions         [t x neuron], each column contains
%                       prediction based on RF for test 
%                       responses of specific neuron (using optimal
%                       lambda)
%   time                [t x 1]; time points for predictions
%
%   traces              [trTime x neuron]; calcium traces of neurons
%   traceTimes          [trTime x 1]; sample times of calcium traces
%   stimFrames          [time x rows x cols]; noise stimulus
%   stimTimes           [time x 1]; times of stimulus frames
%   RFtimesInFrames     [1 x RFframes]; frames of receptive field relative
%                       to stimulus frames
%   lambdas             [1 x lambda]; values of lambda
%   crossFolds          ind; number of cross val. folds

% generate toeplitz matrix for stimuli: [time x pixels]
% each row holds all pixels at current and previous time points:
% [[all pixels at t=0], [all pixels at t=-1], ...]
% each column is time series of that particular pixel

% find time gaps in stimulus presentation (usually when same visual noise
% stimulus was repeated several times)
stimBin = median(diff(stimTimes));
indGap = find(diff(stimTimes) > 2 * stimBin);
time = stimTimes;
% fill gaps with zeros in stimulus matrix
for g = 1:length(indGap)
    add = round(diff(stimTimes(indGap(g) + [0 1])) ./ stimBin);
    stimFrames = [stimFrames(1:indGap(g),:,:); ...
        zeros(add, size(stimFrames,2), size(stimFrames,3)); ...
        stimFrames(indGap(g)+1:end,:,:)];
    time = [time(1:indGap(g)); ...
        time(indGap(g)) + (1:add)' .* stimBin; ...
        time(indGap(g)+1:end)];
end
% reshape stimulus frames to [time x px]; this represents a single
% "stimulus block", i.e. the pixels to estimate a single time point of the
% receptive field
stim = reshape(stimFrames, size(stimFrames,1), []);
% now concatinate time shifted stimulus blocks; for each time point there
% is a stimulus block for lag=0, another for lag=-1, another for lag=-2,...
st = [];
for t = 1:length(RFtimesInFrames)
    st = [st, ...
        [zeros(max(0,RFtimesInFrames(1)-1+t), size(stim,2)); ...
        stim(max(1,2-RFtimesInFrames(1)-t) : end-RFtimesInFrames(1)-t+1, :)]];
end
stim = st;
clear st

% get neural response
traceBin = median(diff(traceTimes));
numBins = round(stimBin / traceBin);
traces = smoothdata(traces, 1, 'movmean', numBins, 'omitnan');
traces = interp1(traceTimes, traces, time);
% z-score neural response
zTraces = (traces - nanmean(traces,1)) ./ nanstd(traces,0,1);

% delete stim frames for which all neurons have NaN
ind = all(isnan(zTraces),2);
stim(ind,:) = [];
zTraces(ind,:) = [];
time(ind) = [];
% if NaN values < 5% in a neuron, exchange NaNs for 0
ind = any(isnan(zTraces),1) & sum(isnan(zTraces),1)/size(zTraces,1) <= 0.05;
if sum(ind) > 0
    zTraces(:,ind) = fillmissing(zTraces(:,ind),'constant',0);
end
% skip neurons that have only NaN values
valid = ~all(isnan(zTraces),1)';

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

% scale lamdas according to number of samples and number of predictors
lamStim = sqrt(lambdas .* size(stim,1) .* size(stim,2));

% construct spatial smoothing lambda matrix
lamMatrix_stim = krnl.makeLambdaMatrix([size(stimFrames,2), size(stimFrames,3), ...
    length(RFtimesInFrames)], [1 1 0]);
lamMatrix_stim = blkdiag(lamMatrix_stim, lamMatrix_stim);

nPerFold = ceil(size(stim,1) / crossFolds);

explainedVariance = NaN(size(traces,2), length(lamStim), crossFolds);
predictions = NaN(nPerFold, crossFolds, size(traces,2), length(lamStim));

% get variances explained
fprintf('  Folds (of %d) to get expl. var. of RF: ', crossFolds)
for fold = 1:crossFolds
    fprintf('%d ',fold)
    ind = (1:nPerFold) + (fold-1)*nPerFold;
    ind(ind > size(zTraces,1)) = [];
    j = true(size(zTraces,1),1);
    j(ind) = false;
    
    y_train = gpuArray(padarray(zTraces(j,valid), size(lamMatrix_stim,1), 'post'));
    y_test = zTraces(~j,valid);
    x_train = stim2(j,:);
    x_test = stim2(~j,:);

    for lamS = 1:length(lamStim)
        lms = lamMatrix_stim .* lamStim(lamS);
        
        A = gpuArray([x_train; lms]);

        B = gather(A \ y_train);
        pred = x_test * B; % get prediction
        predictions(1:sum(~j), fold, valid, lamS) = pred;
        explainedVariance(valid, lamS, fold) = 1 - ...
            sum((y_test - pred) .^ 2,1) ./ ...
            sum((y_test - mean(zTraces(j, valid),1)) .^2, 1);
    end
end
fprintf('\n')

% determine RFs using all data and optimal lambdas
receptiveFields = NaN(size(stim2,2), size(traces,2));

[~, bestStimLams] = max(mean(explainedVariance, 3), [], 2);
fprintf('  Optimal lambdas (of %d) to get RFs: ', length(lamStim))
for lamS = 1:length(lamStim)
    fprintf('%d ', lamS)
    ind = bestStimLams == lamS & valid;
    if sum(ind) == 0
        continue
    end
    A = [stim2; lamMatrix_stim .* lamStim(lamS)];
    tr = padarray(zTraces(:,ind), size(lamMatrix_stim,1), 'post');
    
    B = gather(gpuArray(A) \ gpuArray(tr));
    receptiveFields(:,ind) = B; % get RF kernel
    predictions(:,:,ind,1) = predictions(:,:,ind,lamS);
end
fprintf('\n')

receptiveFields = reshape(receptiveFields, size(stimFrames,2), ...
    size(stimFrames,3), length(RFtimesInFrames), 2, size(traces,2));

predictions = reshape(predictions(:,:,:,1), [], size(traces,2));
predictions = predictions(1:length(time),:);