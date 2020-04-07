function [receptiveFields, explainedVariance, traces_test, ...
    predictions] = getReceptiveField_smooth(traces, traceTimes, stimFrames, ...
    stimFrameTimes, repetitionTimes, RFtimesInFrames, lambdas, crossFolds, ...
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

addBins = max(-RFtimesInFrames(1),RFtimesInFrames(end));

timeBin = median(diff(stimFrameTimes));
time = repetitionTimes.onset(:)'+(0:size(stimFrames, 3)+addBins-1)'.*timeBin;

% generate toeplitz matrix for sitmuli: [[all pixels at t=0];[all pixels at t=-1];...]
% each column is time series of that particular pixel
stim = reshape(stimFrames, [], size(stimFrames, 3))'; % [time x px], single stimulus block containing time series of all pixels
stim = [stim; NaN(addBins, size(stim,2))]; % add NaNs at end; concatinated stimulus blocks at t<0 will add non-NaN entries in these rows
st = []; % concatinate time shifted stimulus blocks
for t = 1:length(RFtimesInFrames)
    st = [st, ...
        [NaN(max(0,RFtimesInFrames(1)-1+t), size(stim,2)); ...
        stim(max(1,2-RFtimesInFrames(1)-t) : end-t+1, :)]];
end
stim = st;
% set NaNs to zero as if no stimulus was presented before and after each repetition
stim(isnan(stim)) = 0;
% repeat stimulus matrix for each repetition of presentation
stim = repmat(stim, 1, 1, length(repetitionTimes.onset));
% for time bins overlapping between trials, add stimuli to next trial and
% delete those extra time bins
for tr = 1:length(repetitionTimes.onset)-1
    ind = time(:,tr) > time(1,tr+1);
    stim(1:sum(ind),:,tr+1) = stim(1:sum(ind),:,tr+1) + stim(ind,:,tr);
    time(ind,tr) = NaN;
end

% concatenate trials (times and stim), but disregard double counted time
% bins
time = time(:);
stim = permute(stim, [1 3 2]);
stim = reshape(stim, [], size(stim,3));
ind = isnan(time);
time(ind) = [];
stim(ind,:) = [];

% get neural response
traceBin = median(diff(traceTimes));
numBins = round(timeBin / traceBin);
traces = smoothdata(traces, 1, 'movmean', numBins, 'omitnan');
traces = interp1(traceTimes', traces, time);
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
end

% scale lamdas according to number of samples and number of predictors
lambdas = sqrt(lambdas .* size(stim,1) .* size(stim,2));

% construct spatial smoothing lambda matrix
lamMatrix = krnl.makeLambdaMatrix([size(stimFrames,1), size(stimFrames,2), ...
    length(RFtimesInFrames)], [1 1 0]);

% calculate RFs using cross-validation and ridge regression
receptiveFields = NaN(size(stim,2), size(traces,2), crossFolds, ...
    length(lambdas), length(models));
explainedVariance = NaN(size(traces,2), crossFolds, length(lambdas), ...
    length(models));
n = ceil(size(stim,1) / crossFolds);
traces_test = NaN(n, size(traces,2), crossFolds);
predictions = NaN(n, size(traces,2), crossFolds, length(lambdas), length(models));
fprintf('  Folds: ')
for fold = 1:crossFolds
    fprintf('%d ',fold)
    ind = (1:n) + (fold-1)*n;
    ind(ind>size(zTraces,1)) = [];
    j = true(size(zTraces,1),1);
    j(ind) = false;
    traces_test(1:sum(~j),:,fold) = zTraces(~j,:); % testing neural trace
    % training neural trace with padded zeros for ridge regression
    tr_train = [zTraces(j,:); zeros(size(lamMatrix,1), size(traces,2))];
    tr_train = gpuArray(tr_train);
    for m = 1:length(models)
        st_train = stim_models{m}(j,:); % training stimulus
        st_train = gpuArray(st_train);
        st_test = stim_models{m}(~j,:);
        for lam = 1:length(lambdas)
            % add lambdas for ridge regression
            A = [st_train; gpuArray(lamMatrix .* lambdas(lam))];
            
            receptiveFields(:,:,fold,lam,m) = gather(A \ tr_train); % get RF kernel
            pred = st_test * receptiveFields(:,:,fold,lam,m); % get prediction based on RF and testing neural response
            predictions(1:sum(~j),:,fold,lam,m) = pred;
            explainedVariance(:,fold,lam,m) = 1 - sum((zTraces(~j,:)-pred).^2,1) ./ ...
                sum((zTraces(~j,:)-mean(zTraces(j,:),1)).^2,1);
        end
    end
    clear tr_train
end
fprintf('\n')

receptiveFields = reshape(receptiveFields, size(stimFrames,1), ...
    size(stimFrames,2), length(RFtimesInFrames), size(traces,2), crossFolds, ...
    length(lambdas), length(models));