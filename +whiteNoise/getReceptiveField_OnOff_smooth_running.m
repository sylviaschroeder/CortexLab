function [receptiveFields, runKernels, runWin, ...
    explainedVariance, explainedVariance_runOnly, ...
    explainedVariance_stimOnly, predictions, predictions_runOnly, time] = ...
    getReceptiveField_OnOff_smooth_running(traces, traceTimes, ...
    stimFrames, stimFrameTimes, stimTimes, RFtimesInFrames, ...
    runSpeed, runTime, runKrnlLimits, lambdas, crossFolds, ...
    fitSeparate)
%GETRECEPTIVEFIELD Returns spatiotemporal receptive field.
%   [receptiveField, frameTimes] = GETRECEPTIVEFIELD(trace, ...
%    traceTimes, stimFrames, stimFrameTimes, repetitionTimes, RFtype, ...
%    plotResponseTraces, method) calculates the linear RF of the neuron.
%
%   receptiveFields     [rows x cols x RFframes x RFtype x neuron]
%                       containing linear regression solution for x in Ax=B
%                       where A is stimulus [rows x cols x time] and B is
%                       calcium response, for each neuron and stimulus 
%                       model; ridge regression is performed on all data
%                       using the optimal lambda value found with
%                       cross-validation
%   runKernels          [t x neuron x (model)]; containing linear
%                       regression kernel fitting calcium based on running
%                       speed; if run and stimulus kernel are fit
%                       simutaneously, there is a different run kernel for
%                       each stimulus model, otherwise there is only one
%                       run kernel
%   runWin              [1 x t]; time of run kernel relative to neural
%                       response
%   explainedVariance   [neuron x crossFold x lambda x model], each entry:
%                       explained variance for fitted RF and run kernel for
%                       each neuron, cross val. fold, lambda, and model
%   explainedVariance_runOnly   [neuron x crossFold x lambda (x model)], 
%                       each entry: explained variance by run kernel only, 
%                       for each neuron, cross val. fold, lambda, and model
%   predictions         [t x neuron x lambda x model], each column contains
%                       prediction based on RF and run kernel for test 
%                       responses of specific neuron, cross val. fold, 
%                       lambda, and model
%   predictions_runOnly [t x neuron x lambda (x model)], each column contains
%                       prediction based on run kernel only for test 
%                       responses of specific neuron, cross val. fold, 
%                       lambda, (and model if run and stimulus kernel
%                       fitted simultaneously)
%   time                [t x 1]; time points for predictions
%
%   traces              [trTime x neuron]; calcium traces of neurons
%   traceTimes          [trTime x 1]; sample times of calcium traces
%   stimFrames          [rows x cols x time]; noise stimulus
%   stimFrameTimes      [1 x time]; times of stimulus frames (from start of
%                       a single repetition)
%   repetitionTimes     struct; .onset of all stimulus repetitions
%   RFtimesInFrames     [1 x RFframes]; frames of receptive field relative
%                       to stimulus frames
%   runSpeed            [rTime x 1]; running speed
%   runTime             [rTime x 1]; time points of running speed
%   runKrnlLimits       [1 x 2]; time limits of running kernel in s
%   lambdas             [1 x lambda]; values of lambda
%   crossValFolds       ind; number of cross val. folds
%   models              {1 x m}; list of strings: 'linear', 'absolute',
%                       'white', 'black'
%   fitSeparate         logical; if true (default), running kernel is fit
%                       first, then RF using separate lambdas; if false,
%                       running kernel and RF are fit together with same
%                       lambda
%   shiftTest           logical; if true (not default), circ-shift traces
%                       and determine explained variance for these traces

if nargin < 13
    fitSeparate = false;
end
if ~iscell(lambdas)
    lambdas = {lambdas, lambdas};
end

if length(stimTimes) > 1
    addBins = max(-RFtimesInFrames(1),RFtimesInFrames(end));
else
    addBins = 0;
end

timeBin = median(diff(stimFrameTimes));
time = stimTimes.onset(:)'+(0:size(stimFrames, 3)+addBins-1)'.*timeBin;

% generate toeplitz matrix for sitmuli: [[all pixels at t=0];[all pixels at t=-1];...]
% each column is time series of that particular pixel
stim = reshape(stimFrames, [], size(stimFrames, 3))'; % [time x px], single stimulus block containing time series of all pixels
stim = [stim; NaN(addBins, size(stim,2))]; % add NaNs at end; concatinated stimulus blocks at t<0 will add non-NaN entries in these rows
st = []; % concatinate time shifted stimulus blocks
for t = 1:length(RFtimesInFrames)
    st = [st, ...
        [NaN(max(0,RFtimesInFrames(1)-1+t), size(stim,2)); ...
        stim(max(1,2-RFtimesInFrames(1)-t) : end-RFtimesInFrames(1)-t+1, :)]];
end
stim = st;
% set NaNs to zero as if no stimulus was presented before and after each repetition
stim(isnan(stim)) = 0;
% repeat stimulus matrix for each repetition of presentation
stim = repmat(stim, 1, 1, length(stimTimes.onset));
% for time bins overlapping between trials, add stimuli to next trial and
% delete those extra time bins
for tr = 1:length(stimTimes.onset)-1
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

% generate Toeplitz matrix for running speed
runSpeed = medfilt1(runSpeed, 5);
runBin = median(diff(runTime));
numBins = round(timeBin / runBin);
runSpeed = smooth(runSpeed, numBins);
runSpeed = interp1(runTime, runSpeed, time, 'pchip');
[runToepl, ~, runWin] =  krnl.getToeplitz(time, [], [], {runSpeed}, ...
    {runKrnlLimits}, true);
runWin = runWin{1};
% z-score
runToepl = (runToepl - mean(runToepl(:))) ./ std(runToepl);

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
runToepl(ind,:) = [];
time(ind) = [];
% if NaN values < 5% in a neuron, exchange NaNs for 0
ind = any(isnan(zTraces),1) & sum(isnan(zTraces),1)/size(zTraces,1) <= 0.05;
if sum(ind) > 0
    zTraces(:,ind) = fillmissing(zTraces(:,ind),'constant',0);
end
% skip neurons that have only NaN values
valid = ~all(isnan(zTraces),1)';
validInd = find(valid);

% duplicate stimulus matrix to predict ON part (1st half) and OFF
% part (2nd half)
s = stim;
s(stim < 0) = 0;
stim2 = s;
s = stim;
s(stim > 0) = 0;
stim2 = [stim2, s];
stim2 = (stim2 - nanmean(stim2(:))) ./ nanstd(stim2(:)); % normalise each column of stimulus matrix

% scale lamdas according to number of samples and number of predictors
lamStim = sqrt(lambdas{1} .* size(stim,1) .* size(stim,2));
lamRun = sqrt(lambdas{2} .* size(stim,1) .* size(runToepl,2));
% lambdas = sqrt(lambdas .* size(stim,1) .* (size(stim,2)+size(runToepl,2)));

% construct spatial smoothing lambda matrix
lamMatrix_stim = krnl.makeLambdaMatrix([size(stimFrames,1), size(stimFrames,2), ...
    length(RFtimesInFrames)], [1 1 0]);
lamMatrix_stim = blkdiag(lamMatrix_stim, lamMatrix_stim);
lamMatrix_run = krnl.makeLambdaMatrix(size(runToepl,2), 1);

nPerFold = ceil(size(stim,1) / crossFolds);

runKernels = NaN(length(runWin), size(traces,2));
if fitSeparate
    explainedVariance = NaN(size(traces,2), crossFolds, length(lamStim));
    predictions = NaN(nPerFold, crossFolds, size(traces,2), length(lamStim));
    explainedVariance_runOnly = NaN(size(traces,2), crossFolds, length(lamRun));
    predictions_runOnly = NaN(nPerFold, crossFolds, size(traces,2), length(lamRun));
    explainedVariance_stimOnly = NaN(size(traces,2), crossFolds, length(lamStim));
else
    explainedVariance = NaN(size(traces,2), crossFolds, length(lamStim), length(lamRun));
    predictions = NaN(nPerFold, crossFolds, size(traces,2), length(lamStim), length(lamRun));
    explainedVariance_runOnly = NaN(size(traces,2), crossFolds, length(lamStim), length(lamRun));
    predictions_runOnly = NaN(nPerFold, crossFolds, size(traces,2), length(lamStim), length(lamRun));
    explainedVariance_stimOnly = NaN(size(traces,2), crossFolds, length(lamStim), length(lamRun));
end

if fitSeparate
    % get variances explained by running alone
    fprintf('  Folds (of %d) to get expl. var. of running: ', crossFolds)
    for fold = 1:crossFolds
        fprintf('%d ',fold)
        ind = (1:nPerFold) + (fold-1)*nPerFold;
        ind(ind>size(zTraces,1)) = [];
        j = true(size(zTraces,1),1);
        j(ind) = false;
        % training neural trace with padded zeros for ridge regression
        y_train = [zTraces(j,valid); zeros(size(runToepl,2), sum(valid))];
        %         tr_train = gpuArray(tr_train);
        run_test = runToepl(~j,:);
        run_train = runToepl(j,:);
        %         run_train = gpuArray(run_train);
        for lamR = 1:length(lamRun)
            % add lambdas for ridge regression
            A = [run_train; lamMatrix_run .* lamRun(lamR)];
            %             A = [run_train; gpuArray(lamMatrix_run .* lambdas(lam))];
            
            B = A \ y_train;
            %             B = gather(A \ tr_train);
            pred = run_test * B; % get prediction
            predictions_runOnly(1:sum(~j),fold,valid,lamR) = pred;
            explainedVariance_runOnly(valid,fold,lamR) = 1 - ...
                sum((zTraces(~j,valid)-pred).^2,1) ./ ...
                sum((zTraces(~j,valid)-mean(zTraces(j,valid),1)).^2,1);
        end
    end
    fprintf('\n')
    
    % get running kernel based on optimal lambda using all data
    [~, bestLamsRun] = max(mean(explainedVariance_runOnly,2), [], 3);
    fprintf('  Optimal lambdas (of %d) to get running kernels: ', length(lamRun))
    for lamR = 1:length(lamRun)
        fprintf('%d ', lamR)
        ind = bestLamsRun == lamR & valid;
        A = [runToepl; lamMatrix_run .* lamRun(lamR)];
%         A = [runToepl; gpuArray(lamMatrix_run .* lambdas(lam))];
        tr = [zTraces(:,ind); zeros(size(runToepl,2), sum(ind))];
%         tr = gpuArray(tr);
        B = A \ tr;
%         B = gather(A \ tr);
        runKernels(:,ind) = B;
        predictions_runOnly(:,:,ind,1) = predictions_runOnly(:,:,ind,lamR);
    end
    fprintf('\n')
    predRun = runToepl * runKernels;
    residuals = zTraces - predRun;
    residuals = (residuals - nanmean(residuals)) ./ nanstd(residuals);
else
    residuals = zTraces;
end

% get variances explained by stimulus alone and by stimulus and running
fprintf('  Folds (of %d) to get expl. var. of RF: ', crossFolds)
for fold = 1:crossFolds
    fprintf('%d ',fold)
    ind = (1:nPerFold) + (fold-1)*nPerFold;
    ind(ind>size(zTraces,1)) = [];
    j = true(size(zTraces,1),1);
    j(ind) = false;
    
    y_train = gpuArray(padarray(residuals(j,valid), size(lamMatrix_stim,1), 'post'));
    y_test = residuals(~j,valid);
    x_train = stim2(j,:);
    x_test = stim2(~j,:);
    if ~fitSeparate
        y2_train = padarray(y_train, size(lamMatrix_run,1), 'post');
        x2_train = [x_train, runToepl(j,:)];
        x2_test = [x_test, runToepl(~j,:)];
    end
    for lamS = 1:length(lamStim)
        lms = lamMatrix_stim .* lamStim(lamS);
        
        if fitSeparate
            % add lambdas for ridge regression
            % A = [st_train; lamMatrix_stim .* lamStim(lam)];
            A = gpuArray([x_train; lms]);
            % B = A \ tr_train;
            B = gather(A \ y_train);
            pred = x_test * B; % get prediction
            explainedVariance_stimOnly(valid,fold,lamS) = 1 - ...
                sum((y_test - pred).^2,1) ./ ...
                sum((y_test - mean(residuals(j,valid),1)).^2,1);
            for iCell = 1:size(pred,2)
                pred(:,iCell) = pred(:,iCell) + predictions_runOnly(1:sum(~j),fold, ...
                    validInd(iCell),bestLamsRun(validInd(iCell)));
            end
            predictions(1:sum(~j),fold,valid,lamS) = pred;
            explainedVariance(valid,fold,lamS) = 1 - ...
                sum((zTraces(~j,valid) - pred).^2,1) ./ ...
                sum((zTraces(~j,valid) - mean(zTraces(j,valid),1)).^2,1);
        else
            for lamR = 1:length(lamRun)
                lmr = lamMatrix_run .* lamRun(lamR);
                A = gpuArray([x2_train; ...
                    [[lms, zeros(size(lms,1), size(lmr,2))]; ...
                    [zeros(size(lmr,1), size(lms,2)), lmr]]]);
                
                B = gather(A \ y2_train);
                pred = x2_test * B; % get prediction
                predictions(1:sum(~j),fold,valid,lamS,lamR) = pred;
                explainedVariance(valid,fold,lamS,lamR) = 1 - ...
                    sum((y_test-pred).^2,1) ./ ...
                    sum((y_test-mean(zTraces(j,valid),1)).^2,1);
                
                pred = x_test * B(1:size(stim2,2),:);
                explainedVariance_stimOnly(valid,fold,lamS,lamR) = 1 - ...
                    sum((y_test-pred).^2,1) ./ ...
                    sum((y_test-mean(zTraces(j,valid),1)).^2,1);
                
                pred = runToepl(~j,:) * B(size(stim2,2)+1:end,:);
                predictions_runOnly(1:sum(~j),fold,valid,lamS,lamR) = pred;
                explainedVariance_runOnly(valid,fold,lamS,lamR) = 1 - ...
                    sum((y_test-pred).^2,1) ./ ...
                    sum((y_test-mean(zTraces(j,valid),1)).^2,1);
            end
        end
    end
end
fprintf('\n')

% determine RFs (and running kernels if not fitSeparate) using all data and optimal lambdas
receptiveFields = NaN(size(stim2,2), size(traces,2));

[evRun, bestStimLams] = max(mean(explainedVariance, 2), [], 3);
evRun = permute(evRun, [1 4 2 3]); % [neuron x 1/lambdasRun]
bestStimLams = permute(bestStimLams, [1 4 2 3]); % [neuron x 1/lambdasRun]
[~, bestRunLams] = max(evRun, [], 2); % [neuron x 1]
ind = sub2ind(size(bestStimLams), (1:size(bestStimLams,1))', bestRunLams);
bestStimLams = bestStimLams(ind); % [neuron x 1]
fprintf('  Optimal lambdas (of %d) to get RFs: ', length(lamStim))
for lamS = 1:length(lamStim)
    fprintf('%d ', lamS)
    ind1 = bestStimLams == lamS & valid;
    A = [stim2; lamMatrix_stim .* lamStim(lamS)];
    tr = padarray(residuals, size(lamMatrix_stim,1), 'post');
    if fitSeparate
        %             B = A \ tr;
        B = gather(gpuArray(A) \ gpuArray(tr(:,ind1)));
        receptiveFields(:,ind1) = B; % get RF kernel
        predictions(:,:,ind1,1) = predictions(:,:,ind1,lamS);
    else
        for lamR = 1:length(lamRun)
            ind2 = ind1 & bestRunLams == lamR;
            A2 = [[A, padarray(runToepl, size(lamMatrix_stim,1), 'post')]; ...
                padarray(lamMatrix_run .* lamRun(lamR), [0 size(A,2)], 'pre')];
            tr2 = padarray(tr(:,ind2), size(lamMatrix_run,1), 'post');
            B = gather(gpuArray(A2) \ gpuArray(tr2));
            receptiveFields(:,ind2) = B(1:size(stim2,2),:); % get RF kernel
            runKernels(:,ind2) = B(size(stim2,2)+1:end,:);
            predictions(:,:,ind2,1,1) = predictions(:,:,ind2,lamS,lamR);
            predictions_runOnly(:,:,ind2,1,1) = predictions_runOnly(:,:,ind2,lamS,lamR);
        end
    end
end
fprintf('\n')

receptiveFields = reshape(receptiveFields, size(stimFrames,1), ...
    size(stimFrames,2), length(RFtimesInFrames), 2, size(traces,2));

predictions = reshape(predictions(:,:,:,1,1), [], size(traces,2));
predictions = predictions(1:length(time),:);
predictions_runOnly = reshape(predictions_runOnly(:,:,:,1,1), [], size(traces,2));
predictions_runOnly = predictions_runOnly(1:length(time),:);