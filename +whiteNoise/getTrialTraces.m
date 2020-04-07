function [trialTraces, trialTimes] = getTrialTraces(traces, traceTime, ...
    stimTimes, stimFrameTimes, convertToStimulusTime)

traceSR = 1 / median(diff(traceTime));
frameDur = median(diff(stimFrameTimes));
% add 5 sec after last frame of each trial
stimFrameTimes = [stimFrameTimes, stimFrameTimes(end) + ...
    (1:round(5/frameDur)).*frameDur];
if convertToStimulusTime
    trialTraces = NaN(length(stimFrameTimes), length(stimTimes.onset), size(traces,2));
    trialTimes = stimFrameTimes;
    % smooth calcium trace before interpolation
    traces = smoothdata(traces, 1, 'movmean', ceil(frameDur * traceSR), 'omitnan');
else
    stimDur = diff(stimFrameTimes([1 end]));
    trialTraces = NaN(round(stimDur * traceSR), length(stimTimes.onset), size(traces,2));
    trialTimes = (0:size(trialTraces,1)-1) ./ traceSR;
end

if convertToStimulusTime
    t = reshape(stimTimes.onset,1,[]) + stimFrameTimes(:);
%     if length(traces) < 15 * length(stimFrameTimes)
        % this is probably calcium signal -> interpolate trace to times
        % of stimulus frames
    last = find(t <= traceTime(end), 1, 'last');
    trialTraces = NaN(numel(t), size(traces,2));
    for n = find(~all(isnan(traces),1))
        trialTraces(1:last,n) = interp1(traceTime, traces(:,n), t(1:last), 'pchip');
    end
    trialTraces = reshape(trialTraces, size(t,1), size(t,2), []);
%     else % this is probably spike counts -> sum counts within each stimulus frame
%         frameDur = median(diff(stimFrameTimes));
%         binSize = median(diff(traceTime));
%         t = [t; t(end,:) + frameDur];
%         [n,~,b] = histcounts(traceTime,t);
%         nt_first = find(b,1);
%         nt_last = find(b,1,'last');
%         nt = nt_first : nt_last;
%         b = b(nt);
%         binIndices = sparse(b, nt, ones(1,length(b)));
%         tr = binIndices * traces(1:nt_last)'./ (n' .* binSize);
%     end
else
%     [~, startInd] = min(abs(traceTime - stimTimes.onset(rep)));
%     tr = traces(startInd : ...
%         min(size(trialTraces,1)-1 + startInd, length(traces)));
%     trialTraces(1:length(tr), rep, :) = tr;
end