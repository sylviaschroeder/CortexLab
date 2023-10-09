function [explainedVariance, explainedVariance_shifted] = ...
    receptiveFieldShiftTest_inclRunning(traces, traceTimes, ...
    stimFrames, stimFrameTimes, stimTimes, RFtimesInFrames, ...
    runSpeed, runTime, runKernels, runWin, receptiveFields, lambdas, numShifts)
%GETRECEPTIVEFIELD Returns spatiotemporal receptive field.
%   [receptiveField, frameTimes] = GETRECEPTIVEFIELD(trace, ...
%    traceTimes, stimFrames, stimFrameTimes, repetitionTimes, RFtype, ...
%    plotResponseTraces, method) calculates the linear RF of the neuron.
%
%   explainedVariance   [neuron x 1], each entry:
%                       explained variance for fitted RF for each neuron
%   explainedVariance_shifted   [neuron x numShifts], 
%                       each entry: explained variance for each neuron and
%                       shift of neural response to stimulus
%
%   traces              [trTime x neuron]; calcium traces of neurons
%   traceTimes          [trTime x 1]; sample times of calcium traces
%   stimFrames          [rows x cols x time]; noise stimulus
%   stimFrameTimes      [1 x time]; times of stimulus frames (from start of
%                       a single repetition)
%   stimTimes           struct; .onset of all stimulus repetitions
%   RFtimesInFrames     [1 x RFframes]; frames of receptive field relative
%                       to stimulus frames
%   runSpeed            [rTime x 1]; running speed
%   runTime             [rTime x 1]; time points of running speed
%   runKernels          [t x neuron]; containing linear
%                       regression kernel fitting calcium based on running
%                       speed
%   runWin              [1 x t]; time of run kernel relative to neural
%                       response
%   receptiveFields     [rows x columns x t x 2 x neuron]; fitted receptive
%                       fields
%   lambdas             [neuron x 1]; optimal lambda for RF fitting found
%                       with cross-validation
%   numShifts           int; number of shifts to generate null
%                       distributions

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
runToepl =  krnl.getToeplitz(time, [], [], {runSpeed}, {runWin([1 end])}, true);
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
% if NaN values < 5% in a neuron, exchange NaNs for 0
ind = any(isnan(zTraces),1) & sum(isnan(zTraces),1)/size(zTraces,1) <= 0.05;
if sum(ind) > 0
    zTraces(:,ind) = fillmissing(zTraces(:,ind),'constant',0);
end
% skip neurons that have only NaN values
valid = ~all(isnan(zTraces),1)';

% subtract running predictions from traces
pred = runToepl * runKernels;
zTraces = zTraces - pred;

% duplicate stimulus matrix to predict ON part (1st half) and OFF part (2nd half)
s = stim;
s(stim < 0) = 0;
stim2 = s;
s = stim;
s(stim > 0) = 0;
stim2 = [stim2, s];
stim2 = (stim2 - nanmean(stim2(:))) ./ nanstd(stim2(:)); % normalise each column of stimulus matrix

% scale lamdas according to number of samples and number of predictors
lamStim = sqrt(lambdas .* size(stim,1) .* size(stim,2));

% construct spatial smoothing lambda matrix
lamMatrix_stim = krnl.makeLambdaMatrix([size(stimFrames,1), size(stimFrames,2), ...
    length(RFtimesInFrames)], [1 1 0]);
lamMatrix_stim = blkdiag(lamMatrix_stim, lamMatrix_stim);

% reshape receptive fields
receptiveFields = reshape(receptiveFields, [], size(receptiveFields, 5));

% get explained variances on original and shifted data
explainedVariance = NaN(size(traces,2), 1);
explainedVariance_shifted = NaN(size(traces,2), numShifts);

lamValues = unique(lamStim);
shifts = randi(size(zTraces,1), numShifts, 1);
shiftedTraces = NaN(size(zTraces,1), numShifts, size(zTraces,2));
for sh = 1:numShifts
    shiftedTraces(:,sh,:) = circshift(zTraces, shifts(sh), 1);
end
for lam = 1:length(lamValues)
    indNeurons = find((lamStim == lamValues(lam))' & valid);
    if isempty(indNeurons)
        continue
    end
    lms = lamMatrix_stim .* lamValues(lam);
    A = gpuArray([stim2; lms]);
    
%     B = gather(A \ gpuArray(zTraces(:,indNeurons)));
%     pred = stim2 * B; % get prediction
    pred = stim2 * receptiveFields(:,indNeurons);
    explainedVariance(indNeurons) = 1 - ...
        sum((zTraces(:,indNeurons) - pred).^2,1) ./ ...
        sum((zTraces(:,indNeurons) - mean(zTraces(:,indNeurons),1)).^2,1);
    
    for iCell = 1:length(indNeurons)
        tr = shiftedTraces(:,:,indNeurons(iCell));
        B = gather(A \ gpuArray(padarray(tr, size(lms,1), 'post'))); 
        pred = stim2 * B;
        explainedVariance_shifted(indNeurons(iCell), :) = 1 - ...
            sum((tr - pred).^2,1) ./ sum((tr - mean(tr,1)).^2,1);
    end
end