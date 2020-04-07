function [trialTraces, trialTimes] = getTrialTraces(calciumTrace, traceTime, ...
    stimTimes, stimFrameTimes, convertToStimulusTime)

traceSR = 1 / median(diff(traceTime));
frameDur = median(diff(stimFrameTimes));
% add 5 sec after last frame of each trial
stimFrameTimes = [stimFrameTimes, stimFrameTimes(end) + ...
    (1:round(5/frameDur)).*frameDur];
if convertToStimulusTime
    trialTraces = NaN(length(stimFrameTimes), length(stimTimes.onset));
    trialTimes = stimFrameTimes;
    % smooth calcium trace before interpolation
    calciumTrace = smooth(calciumTrace, ceil(frameDur * traceSR));
else
    stimDur = diff(stimFrameTimes([1 end]));
    trialTraces = NaN(round(stimDur * traceSR), length(stimTimes.onset));
    trialTimes = (0:size(trialTraces,1)-1) ./ traceSR;
end

for rep = 1:length(stimTimes.onset)
    if convertToStimulusTime
        if length(calciumTrace) < 15 * length(stimFrameTimes)
            % this is probably calcium signal -> interpolate trace to times
            % of stimulus frames
            t = stimTimes.onset(rep) + stimFrameTimes;
            last = find(t <= traceTime(end), 1, 'last');
            trace = NaN(1, length(stimFrameTimes));
            trace(1:last) = interp1(traceTime, calciumTrace, t(1:last), 'pchip');
        else % this is probably spike counts -> sum counts within each stimulus frame
            frameDur = median(diff(stimFrameTimes));
            binSize = median(diff(traceTime));
            t = stimTimes.onset(rep) + stimFrameTimes;
            t = [t t(end)+frameDur];
            [n,~,b] = histcounts(traceTime,t);
            nt_first = find(b,1);
            nt_last = find(b,1,'last');
            nt = nt_first : nt_last;
            b = b(nt);
            binIndices = sparse(b, nt, ones(1,length(b)));
            trace = binIndices * calciumTrace(1:nt_last)'./ (n' .* binSize);
        end
        trialTraces(:, rep) = trace;
    else
        [~, startInd] = min(abs(traceTime - stimTimes.onset(rep)));
        trace = calciumTrace(startInd : ...
            min(size(trialTraces,1)-1 + startInd, length(calciumTrace)));
        trialTraces(1:length(trace), rep) = trace;
    end
end