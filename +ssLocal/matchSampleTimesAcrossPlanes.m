function [traces, times] = matchSampleTimesAcrossPlanes(infos)
% Resamples neural responses so that sample times are equal across planes.
% infos     structure containing meta structures of all planes

% traces    {1 x planes}, each entry contains neural responses [time x
%           neurons]
% times     [time x 1], time for new resposnes

traces = cell(1, length(infos));
planeTimes = cell(1, length(infos));
planeTraces = cell(1, length(infos));
for iPlane = 1:length(infos)
    planeTimes{iPlane} = ppbox.getFrameTimes(infos(iPlane));
    planeTraces{iPlane} = infos(iPlane).F;
%     planeTraces{iPlane} = infos(iPlane).chData(2).F;
end
refPlane = round(length(infos) / 2);
times = planeTimes{refPlane};
for iPlane = 1:length(infos)
    traces{iPlane} = NaN(length(times), size(planeTraces{iPlane},2));
    ind = ~all(isnan(planeTraces{iPlane}), 1);
    traces{iPlane}(:,ind) = interp1(planeTimes{iPlane}, ...
        planeTraces{iPlane}(:,ind), times, 'pchip');
end