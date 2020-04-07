function [fitResults, adjR2] = ...
    fitKernelIteratively(calciumTrace, stimMatrix, ...
    blankStim, kernelLength, basisFunLength, framesBeforeOnset, doPlot, ...
    doShuffle, stimLength)

% calciumTrace          [time x 1], trace of one neuron
% stimMatrix            [stimuli x time], 1 when stimulus is presented, 0
%                       otherwise
% blankStim             [1 x numBlanks]; IDs of stimuli that are gray
%                       screens
% kernelLength          int; number of samples defining duration of kernel
% basisFunLength        int; number of samples defining duration of basis
%                       functions used to model the continuous baseline
%                       activity (set to 0 or longer than kernelLength)
% framesBeforeOnset     int; number of samples before stimulus onset to
%                       estimate baseline (only used when basisFunLength>0)
% doPlot                if 1, data and kernels are plotted
% doShuffle             if 1, calcium trace is shifted randomly to estimate
%                       null hypothesis of kernels
% stimLength            int, stimulus duration in samples (time *
%                       samplingRate); only important if kernelLength much
%                       longer than stimLength

testRep = 500;
sigThreshold = 0.05;
batch = 20;
numBatches = ceil(testRep/batch);

if nargin < 7
    doPlot = 0;
end
if nargin < 8
    doShuffle = 0;
end
if nargin < 9
    stimLength = kernelLength;
end

% get stimulus onsets (exclude blank)
nonBlanks = setdiff(1:size(stimMatrix,1), blankStim);
stimOnsetsAll = [0 diff(sum(stimMatrix,1))];
stimOnsetsAll = find(stimOnsetsAll == 1);
stimIDsAll = sum(bsxfun(@times, stimMatrix, (1:size(stimMatrix,1))'),1);
stimIDsAll = stimIDsAll(stimOnsetsAll);
stimOnsets = stimOnsetsAll;
stimOnsets(stimIDsAll == blankStim) = [];
stimIDs = stimIDsAll;
stimIDs(stimIDs==blankStim) = [];

[alphas, kernel, baseline, prediction, numBasisFuns] = ...
    models.SumOfPulses(calciumTrace, stimOnsets, kernelLength, ...
    basisFunLength, doPlot, stimLength);
mse = nansum((calciumTrace-prediction).^2);

% shift calcium trace and redo kernel to test whether cell has stimulus
% response
if doShuffle == 1
    shuffled = mod(bsxfun(@plus, randi(length(calciumTrace), 1, testRep), ...
        (0:length(calciumTrace)-1)'), length(calciumTrace));
    shuffled(shuffled==0) = length(calciumTrace);
    shuffled = calciumTrace(shuffled);
    mseShuffled = NaN(testRep,1);
    for iBatch = 1:numBatches
        n = min(batch, testRep-(iBatch-1)*batch);
        mseSh = zeros(1, n);
        sh = shuffled(:,(1:n)+(iBatch-1)*batch);
        parfor iRep = 1:n
            [~,~,~,pred] = models.SumOfPulses(sh(:,iRep), stimOnsets, ...
                kernelLength, basisFunLength, 0);
            mseSh(iRep) = nansum((sh(:,iRep)-pred).^2);
        end
        mseShuffled((1:n)+(iBatch-1)*batch) = mseSh;
        if iBatch<numBatches && sum(mseShuffled<mse)/testRep > sigThreshold
            pValue = 2;
            break
        end
    end
    if ~any(isnan(mseShuffled))
        pValue = find(sort(mseShuffled) > mse, 1) / testRep;
    end
    
    if pValue >= 0.05
        close gcf
        models.SumOfPulses(calciumTrace, ...
            stimOnsets, 0, basisFunLength, doPlot);
        fitResults.kernel = [];
        fitResults.alphaEachTrial = [];
        fitResults.baselineEachTrial = [];
        fitResults.baseline = [];
        fitResults.prediction = [];
        adjR2 = [];
        return
    end
end

if ~isempty(alphas)
    alphaMat = NaN(length(stimOnsetsAll)/size(stimMatrix,1), size(stimMatrix,1));
    % [trials x stimuli]
    for iStim = nonBlanks
        alphaMat(:,iStim) = alphas(stimIDs == iStim);
    end
    alphaMat(:,blankStim) = 0;
else
    alphaMat = [];
end

if isempty(baseline) && ~isempty(prediction)
    baseline = calciumTrace - prediction;
end
if ~isempty(baseline)
    indBaseline = bsxfun(@plus, stimOnsetsAll(:), -framesBeforeOnset:-1);
    baselines = mean(baseline(indBaseline),2);
    baseMat = NaN(length(stimOnsetsAll)/size(stimMatrix,1), size(stimMatrix,1));
    % [trials x stimuli]
    for iStim = 1:size(stimMatrix,1)
        baseMat(:,iStim) = baselines(stimIDsAll == iStim);
    end
else
    baseMat = [];
end

fitResults.kernel = kernel(:);
fitResults.alphaEachTrial = alphaMat;
fitResults.baselineEachTrial = baseMat;
fitResults.baseline = baseline;
fitResults.prediction = prediction;

adjR2 = 1 - (mse / ...
    (length(calciumTrace)-sum(isnan(calciumTrace))-length(kernel)- ...
    numel(alphaMat)-numBasisFuns-1)) / ...
    (nansum((calciumTrace - nanmean(calciumTrace)).^2) / ...
    (length(calciumTrace)-sum(isnan(calciumTrace))-1));