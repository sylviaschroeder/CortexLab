function [receptiveField, RFtime] = getRF_spikes(spikeTimes, stimFrames, ...
    stimFrameTimes, stimOnsets, RFtime, RFtype)

frameDur = median(diff(stimFrameTimes));
allFrameTimes = reshape(stimOnsets(:)' + stimFrameTimes(:), [], 1);
sp = spikeTimes(spikeTimes >= allFrameTimes(1) & ...
    spikeTimes <= allFrameTimes(end) + frameDur);
useStimTime = false;

receptiveField = cell(1, length(RFtype));

if length(RFtime) > 2 % all time bins relative to spike onset are given -> timing is given by RFtime
    timeToSpike = sp(:) + RFtime(:)'; % [spikes x RFtime]
    [~,~,frameIDToSpike] = histcounts(timeToSpike(:), ...
        [allFrameTimes; allFrameTimes(end)+frameDur]);
    frameIDToSpike = mod(frameIDToSpike-1, size(stimFrames,3))+1;
    frameIDToSpike = reshape(frameIDToSpike, size(timeToSpike));
    frameIDToSpike(timeToSpike < allFrameTimes(1)) = NaN;
    bins = (1:size(stimFrames,3)+1) - 0.5;
else % timing is defined by timing of stimulus frames, RFtime just specifies the limits 
    spikecount = histcounts(sp, [allFrameTimes; allFrameTimes(end)+frameDur]);
    binCount = floor(RFtime(1)/frameDur):ceil(RFtime(2)/frameDur);
    RFtime = binCount .* frameDur;
    useStimTime = true;
end

for type = 1:length(RFtype)
    stim = stimFrames;
    switch RFtype{type}
        case {'ON','White'}
            stim(stim < 0) = 0;
        case {'OFF','Black'}
            stim(stim > 1) = 0;
            stim = -stim;
        case 'Absolute'
            stim = abs(stim);
    end
    stim = reshape(stim, size(stim,1) * size(stim,2), []); % [pixels x frameID]
    stim = stim - mean(stim(:));    
    STA = NaN(size(stim,1), length(RFtime));
    for t = 1:length(RFtime)
        if useStimTime
            n = zeros(size(spikecount));
            n(max(1,1+binCount(t)) : min(length(n), length(n)+binCount(t))) = ...
                spikecount(max(1,1-binCount(t)) : min(length(n), length(n)-binCount(t)));
            n = sum(reshape(n, size(stimFrames,3), length(stimOnsets)), 2)';
        else
            n = histcounts(frameIDToSpike(:,t), bins);
        end
        STA(:,t) = sum(stim .* n, 2) ./ sum(n);
    end
%     stimMatrix = stim(frameIDToSpike,:); % [frameID relative to spike x pixels]
%     stimMatrix(timeToSpike < allFrameTimes(1),:) = NaN;
%     stimMatrix = reshape(stimMatrix, size(timeToSpike,1), ...
%         size(timeToSpike,2), size(stimMatrix,2)); % [spikes x RFtime x pixels]
%     STA = permute(nanmean(stimMatrix, 1), [3 2 1]); % [pixels x RFtime]
    STA = reshape(STA, size(stimFrames,1), size(stimFrames,2), []); % [pixRows, pixCols, RFtime]
    receptiveField{type} = STA;
end