function [orientationSelectivity, directionSelectivity] = ...
    getTuningSelectivity_raw(calciumTraces, samplingRate, stimulusMatrix, ...
    stimulusSequence, timeAfterOnset, timeAfterOffset)

limitsInFrames = round([timeAfterOnset timeAfterOffset] * samplingRate);
stimTraces = ssLocal.getTracesPerStimulus(calciumTraces, stimulusMatrix, ...
    limitsInFrames);
timeMeans = nanmean(stimTraces, 4);
trialMedians = nanmedian(timeMeans, 3); % [neurons x stimuli]

orientationSelectivity = NaN(size(calciumTraces,2), 1);
directionSelectivity = NaN(size(calciumTraces,2), 1);

% Parse stimulus information
[orientations, blank] = gratings.getOrientations(stimulusSequence);
trialMedians = trialMedians(:,[orientations(:,2); blank]);
oris = mod(orientations(:,1), 360);

for neuron = 1:size(calciumTraces,2)
    [maxi, maxInd] = max(trialMedians(neuron,1:end-1));
    blankResp = trialMedians(neuron,end);
    if blankResp > maxi
        orientationSelectivity(neuron) = NaN;
        directionSelectivity(neuron) = NaN;
        continue
    end
    orthOriInds = [find(oris == mod(oris(maxInd)-90, 360)), ...
        find(oris == mod(oris(maxInd)+90, 360))];
    if isempty(orthOriInds)
        orthOriResp = NaN;
    else
        orthOriResp = max(0, mean(trialMedians(neuron,orthOriInds)) - ...
            blankResp);
    end
    nullDirInd = find(oris == mod(oris(maxInd)+180, 360));
    if isempty(nullDirInd)
        nullDirResp = NaN;
    else
        nullDirResp = max(0, trialMedians(neuron, nullDirInd) - blankResp);
    end
    orientationSelectivity(neuron) = (maxi - orthOriResp) / (maxi + orthOriResp);
    directionSelectivity(neuron) = (maxi - nullDirResp) / (maxi + nullDirResp);
end