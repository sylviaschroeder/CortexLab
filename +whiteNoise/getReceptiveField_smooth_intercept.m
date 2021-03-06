function [receptiveFields, intercepts, explainedVariance, traces_test, ...
    predictions] = getReceptiveField_smooth_intercept(traces, traceTimes, stimFrames, ...
    stimFrameTimes, repetitionTimes, RFframes, lambdas, crossValFolds, ...
    models)
%GETRECEPTIVEFIELD Returns spatiotemporal receptive field.
%   [receptiveField, frameTimes] = GETRECEPTIVEFIELD(trace, ...
%    traceTimes, stimFrames, stimFrameTimes, repetitionTimes, RFtype, ...
%    plotResponseTraces, method) calculates the linear RF of the neuron.
%
%   receptiveFields     [rows x cols x RFframes x neuron x crossFold x lambda x model]
%                       containing linear regression solution for x in Ax=B
%                       where A is stimulus [rows x cols x time] and B is
%                       calcium response, for each neuron, cross validation
%                       fold, regularisation value of lambda, and model
%   explainedVariance   [neuron x crossFold x lambda x model], each entry:
%                       explained variance for fitted RF for each neuron,
%                       cross val. fold, lambda, and model
%   traces_test         [t x neuron x crossFold], each column contains test
%                       data of neuron and cross val. fold
%   predictions         [t x neuron x crossFold x lambda x model], each
%                       column contains prediction for test set of specific
%                       neuron, cross val. fold, lambda, and model
%
%   traces              [trTime x neuron]; calcium traces of neurons
%   traceTimes          [trTime x 1]; sample times of calcium traces
%   stimFrames          [rows x cols x time]; noise stimulus
%   stimFrameTimes      [1 x time]; times of stimulus frames
%   repetitionTimes     struct; .onset of all stimulus repetitions
%   RFframes            [1 x RFframes]; frames of receptive field relative
%                       to stimulus frames
%   lambdas             [1 x lambda]; values of lambda
%   crossValFolds       ind; number of cross val. folds
%   models              {1 x m}; list of strings: 'linear', 'absolute',
%                       'white', 'black'

% generate toeplitz matrix: [[all pixels at t=0];[all pixels at t=-1];...]
% each column is time series of that particular pixel
stim = reshape(stimFrames, [], size(stimFrames, 3))'; % [time x px], single stimulus block containing time series of all pixels
stim = [stim; NaN(max(-RFframes(1),RFframes(end)), size(stim,2))]; % add NaNs at end; concatinated stimulus blocks at t<0 will add non-NaN entries in these rows
st = []; % concatinate time shifted stimulus blocks
for t = 1:length(RFframes)
    st = [st, ...
        [NaN(max(0,RFframes(1)-1+t), size(stim,2)); ...
        stim(max(1,2-RFframes(1)-t) : end-t+1, :)]];
end
stim = st;
stim(isnan(stim)) = 0; % set NaNs to zero as if no stimulus was presented before and after each repetition
stim = repmat(stim, length(repetitionTimes.onset), 1); % repeat stimulus matrix for each repetition of presentation

% get neural response for each repetition
trialTraces = whiteNoise.getTrialTraces(traces, traceTimes, repetitionTimes, ...
    stimFrameTimes, 1); % [time x repetition x neuron]
trialTraces(size(st,1)+1:end,:,:) = []; % cut each time series to length of stimulus presentation

% z-score neural response
trialTraces = reshape(trialTraces, [], size(trialTraces,3));
zTrace = trialTraces;
% zTrace = (trialTraces - nanmean(trialTraces,1)) ./ nanstd(trialTraces,0,1);

% delete stim frames for which all neurons have NaN
ind = all(isnan(zTrace),2);
stim(ind,:) = [];
zTrace(ind,:) = [];
% if NaN values < 5% in a neuron, exchange NaNs for 0
ind = any(isnan(zTrace),1) & sum(isnan(zTrace),1)/size(zTrace,1) <= 0.05;
if sum(ind) > 0
    zTrace(:,ind) = fillmissing(zTrace(:,ind),'constant',nanmean(zTrace(:,ind)));
end

stim_models = cell(1, length(models));
for m = 1:length(models)
    switch models{m}
        case 'linear'
            stim_models{m} = stim;
        case 'absolute'
            stim_models{m} = abs(stim);
        case 'white'
            stim_models{m} = stim;
            stim_models{m}(stim < 0) = 0;
        case 'black'
            stim_models{m} = stim;
            stim_models{m}(stim > 0) = 0;
        otherwise
            disp(['Error: Model ' models{m} ' not known'])
            return
    end
    % normalise stimulus matrix
    stim_models{m} = (stim_models{m} - nanmean(stim_models{m}(:))) ./ ...
        nanstd(stim_models{m}(:));
    % add ones for intercept estimation
    stim_models{m} = [stim_models{m} ones(size(stim,1),1)];
end

% scale lamdas according to number of predictors
lambdas = lambdas .* size(stim,2);

% construct spatial smoothing lambda matrix
lamMatrix = whiteNoise.makeLambdaMatrix(size(stimFrames,1), ...
    size(stimFrames,2), length(RFframes));
% add column of zeros for intercept (not regularised)
lamMatrix = [lamMatrix zeros(size(lamMatrix,1),1)];

% calculate RFs using cross-validation and ridge regression
receptiveFields = NaN(size(stim,2)+1, size(traces,2), crossValFolds, ...
    length(lambdas), length(models));
explainedVariance = NaN(size(traces,2), crossValFolds, length(lambdas), ...
    length(models));
n = ceil(size(stim,1) / crossValFolds);
traces_test = NaN(n, size(traces,2), crossValFolds);
predictions = NaN(n, size(traces,2), crossValFolds, length(lambdas), length(models));
fprintf('  Folds: ')
for fold = 1:crossValFolds
    fprintf('%d ',fold)
    ind = (1:n) + (fold-1)*n;
    ind(ind>size(zTrace,1)) = [];
    j = true(size(zTrace,1),1);
    j(ind) = false;
    traces_test(1:sum(~j),:,fold) = zTrace(~j,:); % testing neural trace
    tr_train = [zTrace(j,:); zeros(size(lamMatrix,1), size(traces,2))]; % training neural trace with padded zeros for ridge regression
    tr_train = gpuArray(tr_train);
    for m = 1:length(models)
        st_train = stim_models{m}(j,:); % training stimulus
        st_train = gpuArray(st_train);
        st_test = stim_models{m}(~j,:);
        for lam = 1:length(lambdas)
            st_lam = [st_train; gpuArray(lamMatrix .* lambdas(lam))]; % add lambdas for ridge regression
            
            receptiveFields(:,:,fold,lam,m) = gather(st_lam \ tr_train); % get RF kernel
            pred = st_test * receptiveFields(:,:,fold,lam,m); % get prediction based on RF and testing neural response
            predictions(1:sum(~j),:,fold,lam,m) = pred;
            explainedVariance(:,fold,lam,m) = 1 - sum((zTrace(~j,:)-pred).^2,1) ./ ...
                sum((zTrace(~j,:)-mean(zTrace(j,:),1)).^2,1);
        end
    end
    clear tr_train
end
fprintf('\n')

intercepts = squeeze(receptiveFields(end,:,:,:,:));
receptiveFields = reshape(receptiveFields(1:end-1,:,:,:,:), size(stimFrames,1), ...
    size(stimFrames,2), length(RFframes), size(traces,2), crossValFolds, ...
    length(lambdas), length(models));