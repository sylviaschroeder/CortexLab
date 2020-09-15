function [aligned, t] = getAlignedTraces(signal, time, eventTimes, window)
% aligned   % [t x trials x numSignals]

aligned = [];
t = [];

dt = median(diff(time));
w = round(window / dt);
mini = min(w(:,1));
maxi = max(w(:,2));
tmp = mini : maxi;
if size(window,1) == 1
    tInd = repmat(tmp, length(eventTimes), 1);
else
    if size(window,1) ~= length(eventTimes)
        disp('No. events and windows have to be the same!')
        return
    end
    tInd = repmat(tmp, size(window,1), 1);
    for j = 1:size(window,1)
        tInd(j, tmp<w(j,1) | tmp>w(j,2)) = NaN;
    end
end
t = tmp .* dt;

aligned = NaN(length(t), length(eventTimes), size(signal, 2));

for ev = 1:length(eventTimes)
    evInd = find(time > eventTimes(ev), 1);
    if isempty(evInd)
        continue
    end
    if evInd>1 && eventTimes(ev)-time(evInd-1) < time(evInd)-eventTimes(ev)
        evInd = evInd - 1;
    end
    inds = evInd + tInd(ev,:);
    valid = inds > 0 & inds <= length(time);
    aligned(valid,ev,:) = signal(inds(valid),:);
end