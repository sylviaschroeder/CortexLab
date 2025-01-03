function [receptiveField, explainedVariance, traces, ...
    predictions] = getReceptiveField(trace, traceTimes, stimFrames, ...
    stimFrameTimes, repetitionTimes, RFframes, ...
    lambdas, crossValFolds, plotResponseTraces)
%GETRECEPTIVEFIELD Returns spatiotemporal receptive field.
%   [receptiveField, frameTimes] = GETRECEPTIVEFIELD(trace, ...
%    traceTimes, stimFrames, stimFrameTimes, repetitionTimes, RFtype, ...
%    plotResponseTraces, method) calculates the linear RF of the neuron.
%
%   receptiveField      {1 x RFtypes}; in each entry: [rows x cols x time]
%                       containing Pearson's correlation between the pixel
%                       of the stimulus with the response at the specified 
%                       time lag 
%   frameTimes          in sec
%
%   trace               [n x 1]; calcium trace of neuron
%   traceTimes          [n x 1]; sample times of calcium trace
%   stimFrames          [rows x cols x m]; noise stimulus
%   stimFrameTimes      [1 x m]; times of stimulus frames
%   repetitionTimes     struct; contains stimulus on- and offsets of all
%                       stimulus repetitions
%   plotResponseTraces  0 or 1

% generate toeplitz matrix: [time x pixels]
% each row holds all pixels at current and previous time points:
% [[all pixels at t=0], [all pixels at t=-1], ...]
% each column is time series of that particular pixel
stim = reshape(stimFrames, [], size(stimFrames, 3))'; % single stimulus block containing time series of all pixels
stim = [stim; NaN(length(RFframes)-1,size(stim,2))]; % add NaNs at end; concatinated stimulus blocks at t<0 will add non-NaN entries in these rows
st = []; % concatinate time shifted stimulus blocks
for t = 1:length(RFframes)
    st = [st, [NaN(RFframes(1)-2+t, size(stim,2)); ...
        stim(max(2-RFframes(1)-t,1):end-t+1,:)]];
end
stim = st;
stim(isnan(stim)) = 0; % set NaNs to zero as if no stimulus was presented before and after each repetition
stim = repmat(stim, length(repetitionTimes.onset), 1); % repeat stimulus matrix for each repetition of presentation

% get neural response for each repetition
trialTraces = whiteNoise.getTrialTraces(trace, traceTimes, repetitionTimes, ...
    stimFrameTimes, 1);
trialTraces(size(st,1)+1:end,:) = []; % cut each time series to length of stimulus presentation

% z-score neural response
zTrace = (trialTraces - nanmean(trialTraces(:))) ./ nanstd(trialTraces(:));
zTrace = zTrace(:);

if plotResponseTraces == 1
    screenSize = get(0, 'ScreenSize');
    figure('Position', [10 695 screenSize(3)-20 420])
    hold on
    for rep = 1:length(repetitionTimes.onset)
        plot(stimFrameTimes, trialTraces(1:length(stimFrameTimes),rep), 'Color', ...
            ones(1, 3) * (rep - 1) / length(repetitionTimes.onset))
    end
    plot(stimFrameTimes, median(trialTraces(1:length(stimFrameTimes),:),2), ...
        'r', 'LineWidth', 2)
    axis tight
    xlabel('Time (in s)')
    ylabel('Raw Ca-response')
end

% delete entries where response to NaN
ind = isnan(zTrace);
stim(ind,:) = [];
zTrace(ind) = [];

% duplicate stimulus matrix to predict linear part (1st half) and complex
% part (2nd half, where stimulus was set to absolute values)
stim = [stim, abs(stim)];
stim = (stim - nanmean(stim(:))) ./ nanstd(stim(:)); % normalise each column of stimulus matrix

% scale lamdas according to number of predictors
lambdas = lambdas .* size(stim,2);

% calculate RFs using cross-validation and ridge regression
receptiveField = NaN(size(stim,2), crossValFolds, length(lambdas));
explainedVariance = NaN(crossValFolds, length(lambdas));
n = ceil(size(stim,1) / crossValFolds);
traces = NaN(n, crossValFolds);
predictions = NaN(n, crossValFolds, length(lambdas));
fprintf('  Folds: ')
for fold = 1:crossValFolds
    fprintf('%d ',fold)
    ind = (1:n) + (fold-1)*n;
    ind(ind>length(zTrace)) = [];
    j = true(length(zTrace),1);
    j(ind) = false;
    st_fold = stim(j,:); % training stimulus
    st_fold = gpuArray(st_fold);
    traces(1:sum(~j),fold) = zTrace(~j); % testing neural trace
    tr_fold = [zTrace(j); zeros(size(st_fold,2),1)]; % training neural trace with padded zeros for ridge regression
    tr_fold = gpuArray(tr_fold);
    for lam = 1:length(lambdas)
        st_lam = [st_fold; gpuArray(eye(size(st_fold,2)) .* lambdas(lam))]; % add lambdas for ridge regression
%         st_gpu = gpuArray(st_lam);
        
        receptiveField(:,fold,lam) = gather(st_lam \ tr_fold); % get RF kernel
        pred = stim(~j,:) * receptiveField(:,fold,lam); % get predicion based on RF and testing neural response
        predictions(1:sum(~j),fold,lam) = pred;
        explainedVariance(fold,lam) = 1 - sum((zTrace(~j)-pred).^2) / ...
            sum((zTrace(~j)-mean(zTrace(j))).^2);
    end
    clear tr_fold
end
fprintf('\n')