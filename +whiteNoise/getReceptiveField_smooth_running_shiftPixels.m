function [receptiveFields, runKernels, runWin, ...
    explainedVariance, explainedVariance_runOnly, ...
    explainedVariance_stimOnly, explainedVariance_stimOnly_shifted, ...
    predictions, predictions_runOnly, time] = ...
    getReceptiveField_smooth_running_shiftPixels(traces, traceTimes, ...
    stimFrames, stimFrameTimes, stimTimes, RFtimesInFrames, ...
    runSpeed, runTime, runKrnlLimits, lambdas, crossFolds, models, ...
    fitSeparate, shiftTest)
%GETRECEPTIVEFIELD Returns spatiotemporal receptive field.
%   [receptiveField, frameTimes] = GETRECEPTIVEFIELD(trace, ...
%    traceTimes, stimFrames, stimFrameTimes, repetitionTimes, RFtype, ...
%    plotResponseTraces, method) calculates the linear RF of the neuron.
%
%   receptiveFields     [rows x cols x RFframes x neuron x model]
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
    fitSeparate = true;
end
if nargin < 14
    shiftTest = false;
end
if ~iscell(lambdas)
    lambdas = {lambdas, lambdas};
end
numShifts = 200;
batchSize = 200;
explainedVariance_stimOnly_shifted = [];

addBins = max(-RFtimesInFrames(1),RFtimesInFrames(end));

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
lamStim = sqrt(lambdas{1} .* size(stim,1) .* size(stim,2));
lamRun = sqrt(lambdas{2} .* size(stim,1) .* size(runToepl,2));
% lambdas = sqrt(lambdas .* size(stim,1) .* (size(stim,2)+size(runToepl,2)));

% construct spatial smoothing lambda matrix
lamMatrix_stim = krnl.makeLambdaMatrix([size(stimFrames,1), size(stimFrames,2), ...
    length(RFtimesInFrames)], [1 1 1]);
lamMatrix_run = krnl.makeLambdaMatrix(size(runToepl,2), 1);

nPerFold = ceil(size(stim,1) / crossFolds);

if fitSeparate
    explainedVariance_runOnly = NaN(size(traces,2), crossFolds, length(lamRun));
    predictions_runOnly = NaN(nPerFold, crossFolds, size(traces,2), length(lamRun));
    runKernels = NaN(length(runWin), size(traces,2));
    
    fprintf('  Folds (of %d) to get expl. var. of running: ', crossFolds)
    for fold = 1:crossFolds
        fprintf('%d ',fold)
        ind = (1:nPerFold) + (fold-1)*nPerFold;
        ind(ind>size(zTraces,1)) = [];
        j = true(size(zTraces,1),1);
        j(ind) = false;
        % training neural trace with padded zeros for ridge regression
        tr_train = [zTraces(j,valid); zeros(size(runToepl,2), sum(valid))];
%         tr_train = gpuArray(tr_train);
        run_test = runToepl(~j,:);
        run_train = runToepl(j,:);
%         run_train = gpuArray(run_train);
        for lam = 1:length(lamRun)
            % add lambdas for ridge regression
            A = [run_train; lamMatrix_run .* lamRun(lam)];
%             A = [run_train; gpuArray(lamMatrix_run .* lambdas(lam))];
            
            B = A \ tr_train;
%             B = gather(A \ tr_train);
            pred = run_test * B; % get prediction
            predictions_runOnly(1:sum(~j),fold,valid,lam) = pred(:,:,1);
            explainedVariance_runOnly(valid,fold,lam) = 1 - ...
                sum((zTraces(~j,valid)-pred).^2,1) ./ ...
                sum((zTraces(~j,valid)-mean(zTraces(j,valid),1)).^2,1);
        end
    end
    fprintf('\n')
    [~, bestLamsRun] = max(mean(explainedVariance_runOnly,2), [], 3);
    fprintf('  Optimal lambdas (of %d) to get running kernels: ', length(lamRun))
    for lam = 1:length(lamRun)
        fprintf('%d ', lam)
        ind = bestLamsRun == lam & valid;
        A = [runToepl; lamMatrix_run .* lamRun(lam)];
%         A = [runToepl; gpuArray(lamMatrix_run .* lambdas(lam))];
        tr = [zTraces(:,ind); zeros(size(runToepl,2), sum(ind))];
%         tr = gpuArray(tr);
        B = A \ tr;
%         B = gather(A \ tr);
        runKernels(:,ind) = B;
    end
    fprintf('\n')
    predRun = runToepl * runKernels;
    residuals = zTraces - predRun;
    residuals = (residuals - nanmean(residuals)) ./ nanstd(residuals);
else
    explainedVariance_runOnly = NaN(size(traces,2), crossFolds, length(lamRun), ...
        length(models));
    predictions_runOnly = NaN(nPerFold, crossFolds, size(traces,2), length(lamStim), length(models));
    runKernels = NaN(length(runWin), size(traces,2), length(models));
end

% find best lambda values using cross-validation and ridge regression, get
% explainedVariances and predictions for all lambda values
explainedVariance = NaN(size(traces,2), crossFolds, length(lamStim), ...
    length(models));
explainedVariance_stimOnly = NaN(size(traces,2), crossFolds, length(lamStim), ...
    length(models));
predictions = NaN(nPerFold, crossFolds, size(traces,2), length(lamStim), length(models));

% fit RF (and running kernel if not fit before)
fprintf('  Folds (of %d) to get expl. var. of RF: ', crossFolds)
for fold = 1:crossFolds
    fprintf('%d ',fold)
    ind = (1:nPerFold) + (fold-1)*nPerFold;
    ind(ind>size(zTraces,1)) = [];
    j = true(size(zTraces,1),1);
    j(ind) = false;
    % training neural trace with padded zeros for ridge regression
    if fitSeparate
        tr_train = [residuals(j,valid); zeros(size(stim,2), sum(valid))];
    else
        tr_train = [zTraces(j,valid); zeros(size(stim,2) + size(runToepl,2), sum(valid))];
    end
    tr_train = gpuArray(tr_train);
    run_test = runToepl(~j,:);
    run_train = runToepl(j,:);
    run_train = gpuArray(run_train);
    for m = 1:length(models)
        st_train = stim_models{m}(j,:); % training stimulus
        st_train = gpuArray(st_train);
        st_test = stim_models{m}(~j,:);
        for lam = 1:length(lamStim)
            if fitSeparate
                % add lambdas for ridge regression
%                 A = [st_train; lamMatrix_stim .* lamStim(lam)];
                A = [st_train; gpuArray(lamMatrix_stim .* lamStim(lam))];
                
%                 B = A \ tr_train;
                B = gather(A \ tr_train);
                pred = st_test * B; % get prediction
                explainedVariance_stimOnly(valid,fold,lam,m) = 1 - ...
                    sum((residuals(~j,valid)-pred).^2,1) ./ ...
                    sum((residuals(~j,valid)-mean(residuals(j,valid),1)).^2,1);
                for iCell = 1:size(pred,2)
                    pred(:,iCell) = pred(:,iCell) + predictions_runOnly(1:sum(~j),fold, ...
                        validInd(iCell),bestLamsRun(validInd(iCell)));
                end
                predictions(1:sum(~j),fold,valid,lam,m) = pred;
                explainedVariance(valid,fold,lam,m) = 1 - ...
                    sum((zTraces(~j,valid)-pred).^2,1) ./ ...
                    sum((zTraces(~j,valid)-mean(zTraces(j,valid),1)).^2,1);
            else
                % add lambdas for ridge regression (using lamStim for stim
                % and running!!)
                A = [[st_train, run_train]; ...
                    gpuArray([[lamMatrix_stim .* lamStim(lam), ...
                    zeros(size(lamMatrix_stim,1), size(runToepl,2))]; ...
                    [zeros(size(lamMatrix_run,1), size(stim,2)), ...
                    lamMatrix_run .* lamStim(lam)]])];
                
                B = gather(A \ tr_train);
                pred = [st_test, run_test] * B; % get prediction
                predictions(1:sum(~j),fold,:,lam,m) = pred;
                explainedVariance(:,fold,lam,m) = 1 - sum((zTraces(~j,:)-pred).^2,1) ./ ...
                    sum((zTraces(~j,:)-mean(zTraces(j,:),1)).^2,1);
                pred = run_test * B(size(stim,2)+1:end,:); % get prediction for running only
                predictions_runOnly(1:sum(~j),fold,:,lam,m) = pred;
                explainedVariance_runOnly(:,fold,lam,m) = 1 - sum((zTraces(~j,:)-pred).^2,1) ./ ...
                    sum((zTraces(~j,:)-mean(zTraces(j,:),1)).^2,1);
            end
        end
    end
end
fprintf('\n')

% determine RFs (and running kernels if not before) using all data and optimal lambdas
receptiveFields = NaN(size(stim,2), size(traces,2), length(models));

[~, bestLams] = max(mean(explainedVariance, 2), [], 3);
bestLams = squeeze(bestLams); % [neuron x model]
fprintf('  Optimal lambdas (of %d) to get RFs: ', length(lamStim))
for lam = 1:length(lamStim)
    fprintf('%d ', lam)
    for m = 1:length(models)
        ind = bestLams(:,m) == lam & valid;
        if fitSeparate
            A = [stim_models{m}; lamMatrix_stim .* lamStim(lam)];
            A = gpuArray(A);
            tr = [residuals(:,ind); zeros(size(stim,2), sum(ind))];
            tr = gpuArray(tr);
%             B = A \ tr;
            B = gather(A \ tr);
            receptiveFields(:,ind,m) = B; % get RF kernel
        else %(using lamStim for stim and running!!)
            A = gpuArray([[stim_models{m}, runToepl]; ...
                [[lamMatrix_stim .* lamStim(lam), ...
                zeros(size(lamMatrix_stim,1), size(runToepl,2))]; ...
                [zeros(size(lamMatrix_run,1), size(stim,2)), ...
                lamMatrix_run .* lamStim(lam)]]]);
            tr = gpuArray([zTraces(:,ind); zeros(size(stim,2) + size(runToepl,2), sum(ind))]);
            B = gather(A \ tr);
            receptiveFields(:,ind,m) = B(1:size(stim,2),:); % get RF kernel
            runKernels(:,ind,m) = B(size(stim,2)+1:end,:);
        end
    end
end
fprintf('\n')

receptiveFields = reshape(receptiveFields, size(stimFrames,1), ...
    size(stimFrames,2), length(RFtimesInFrames), size(traces,2), length(models));

if shiftTest && fitSeparate % shiftTest not implemented when fitting together
    shifts = randi(size(zTraces,1), size(stim,2), numShifts);
    explainedVariance_stimOnly_shifted = NaN(size(traces,2), crossFolds, ...
        length(models), numShifts);
    gDev = gpuDevice;
    reset(gDev);
    fprintf('  Shifts (of %d) to get expl. var. of RF on shifted pixels: \n', numShifts)
    for sh = 1:numShifts
        if mod(sh,10) == 0
            fprintf('%d ', sh)
        end
        shStim = NaN(size(stim,1), size(stim,2));
        for px = 1:size(stim,2)
            shStim(:,px) = circshift(stim(:,px), shifts(px,sh));
        end
        shSt = cell(1, length(models));
        for m = 1:length(models)
            switch models{m}
                case 'linear'
                    shSt{m} = shStim;
                case 'absolute'
                    shSt{m} = abs(shStim);
                case 'white'
                    shSt{m} = shStim;
                    shSt{m}(shSt{m} < 0) = 0;
                case 'black'
                    shSt{m} = shStim;
                    shSt{m}(shSt{m} > 0) = 0;
            end
            % normalise stimulus matrix
            shSt{m} = (shSt{m} - nanmean(shSt{m}(:))) ./ nanstd(shSt{m}(:));
        end
        
        % calculate explained variances on shifted residuals using optimal lambda
        for fold = 1:crossFolds
            ind = (1:nPerFold) + (fold-1)*nPerFold;
            ind(ind>size(zTraces,1)) = [];
            j = true(size(zTraces,1),1);
            j(ind) = false;
            % training neural trace with padded zeros for ridge regression
            tr_train = [residuals(j,:); zeros(size(stim,2), size(residuals,2))];
            tr_train = gpuArray(tr_train);
            for m = 4%1:length(models)
                st_train = shSt{m}(j,:); % training stimulus
                st_test = shSt{m}(~j,:);
                for lam = 3%1:length(lamStim)
                    indCells = find(bestLams(:,m) == lam & valid);
                    A = gpuArray([st_train; lamMatrix_stim .* lamStim(lam)]);
                    B = gather(A \ tr_train(:,indCells));
                    pred = st_test * B;
                    explainedVariance_stimOnly_shifted(indCells,fold,m,sh) = ...
                        1 - sum((residuals(~j,indCells)-pred).^2,1) ./ ...
                        sum((residuals(~j,indCells)-mean(residuals(j,indCells),1)).^2,1);
                end
                fprintf(' ')
            end
        end
        fprintf('\n')
    end
end

predictions = reshape(predictions, [], size(traces,2), length(lamStim), ...
    length(models));
predictions = predictions(1:length(time),:,:,:);
if fitSeparate
    predictions_runOnly = reshape(predictions_runOnly, [], size(traces,2), ...
        length(lamRun));
    predictions_runOnly = predictions_runOnly(1:length(time),:,:);
else
    predictions_runOnly = reshape(predictions_runOnly, [], size(traces,2), ...
        length(lamStim), length(models));
    predictions_runOnly = predictions_runOnly(1:length(time),:,:,:);
end