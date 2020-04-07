function traces = getDeltaFOverFTraces(rawTraces, samplingRate)

windowLimits = [-15 15]; % in sec; window size (start and end relative to 
                         % current time point) for calculation of baseline
baselineFun = @(x) prctile(x, 8, 1); % function to calculate F_0 in (F - F_0)/F_0;
                                  % take care about NaN-values!

windowLimitsInFrames = round(windowLimits * samplingRate);
movingWindows = repmat((windowLimitsInFrames(1):windowLimitsInFrames(2))', ...
    1, size(rawTraces, 1)) + repmat(0:size(rawTraces,1)-1, ...
    diff(windowLimitsInFrames)+1, 1);
nanInd = movingWindows < 1 | movingWindows > size(rawTraces, 1);
movingWindows(nanInd) = 1;
traces = NaN(size(rawTraces));
for neuron = 1:size(rawTraces, 2)
    baseline = rawTraces(movingWindows, neuron);
    baseline = reshape(baseline, size(movingWindows));
    baseline(nanInd) = NaN;
    baseline = baselineFun(baseline)';
    traces(:,neuron) = (rawTraces(:,neuron) - baseline) ./ baseline;
end