function periods = getThresholdedPeriods(data, time, thresholds, ...
    minInterval)

% periods       [1 x length(data)], classifies data according to given
%               thresholds, each entry is 1, 2, ..., or
%               length(thresholds)+1, or 0 if interval of same class is
%               shorter than minInterval
%
% data          [1 x t], data values to be thresholded
% time          [1 x t], time points of data (in sec)
% thresholds    [1 x threshs]
% minInterval   in sec; minimum length of threshold excess to be considered

periods = zeros(1, length(data));
data = data(:);
time = time(:);

threshs = [min(data)-1; thresholds(:); max(data)+1];
for k = 1:length(threshs)-1
    periods(data > threshs(k) & data <= threshs(k+1)) = k;
end

periodStartInds = find([1, diff(periods(1:end-1))~=0, 1]);
periodStartTimes = time(periodStartInds);
intervalLengths = diff(periodStartTimes);
exclude = find(intervalLengths < minInterval);
for k = exclude'
    periods(periodStartInds(k) : periodStartInds(k+1)) = 0;
end

% periodLimits = diff(periods);
% periodStarts = find(periodLimits == 1) + 1;
% periodEnds = find(periodLimits == -1);
% if periodStarts(1) > periodEnds(1)
%     periodStarts = [1; periodStarts];
% end
% if periodEnds(end) < periodStarts(end)
%     periodEnds(end+1) = length(periods);
% end
% 
% intervalLengths = time(periodEnds) - time(periodStarts);
% exclude = find(intervalLengths < minInterval);
% for i = exclude'
%     periods(periodStarts(i):periodEnds(i)) = 0;
% end
% 
% periods = double(periods);