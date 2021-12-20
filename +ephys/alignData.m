function [times_aligned, trials] = alignData(dataTimes, alignTimes, limits)

% Aligns times of some events (e.g. spike times) to other events (e.g. stimulus presentation times).
% Events that will be aligned are called ev1.
% Events that are the "anchors" are called ev2.
%
% times_aligned     [ev1_al x 1]; the time of each ev1 relative to any ev2
%                   +-limits
% trials            [ev1_al x 1]; the identity of ev2 that each ev1 was
%                   aligned to
%
% dataTimes         [ev1 x 1]; the absolute time of each ev1
% alignTimes        [ev2 x 1]; the absolute time of each ev2
% limits            [pre post]; relative time spans before and after each
%                   ev2, in which aligned ev1 can fall

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