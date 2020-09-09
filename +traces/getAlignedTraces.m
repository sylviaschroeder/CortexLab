function [aligned, t] = getAlignedTraces(signal, time, eventTimes, window)

dt = median(diff(time));
tInd = round(window(1)/dt) : round(window(2)/dt);
t = tInd .* dt;

aligned = NaN(length(t), length(eventTimes), size(signal, 2));

for ev = 1:length(eventTimes)
    evInd = find(time > eventTimes(ev), 1);
    if isempty(evInd)
        continue
    end
    if evInd>1 && eventTimes(ev)-time(evInd-1) < time(evInd)-eventTimes(ev)
        evInd = evInd - 1;
    end
    inds = evInd + tInd;
    valid = inds > 0 & inds <= length(time);
    aligned(valid,ev,:) = signal(inds(valid),:);
end