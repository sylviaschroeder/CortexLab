function [times_aligned, trials] = alignData(dataTimes, alignTimes, limits)

times_aligned = [];
trials = [];

for j = 1:length(alignTimes)
    if isnan(alignTimes(j))
        continue
    end
    dt = dataTimes - alignTimes(j);
    dt(dt < limits(1) | dt > limits(2)) = [];
    times_aligned = [times_aligned; dt];
    trials = [trials; ones(length(dt),1).*j];
end