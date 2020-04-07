function f_zeros = extractBaselines(traces, samplingRate, percentile, ...
    movingWindow)

movingWindow = round(movingWindow * samplingRate);
halfWindow = ceil(movingWindow / 2);

f_zeros = zeros(size(traces));

indices = false(size(traces, 1), size(traces, 1));
indices(1:movingWindow, 1:halfWindow) = true;
indices(end-movingWindow+1:end, end-halfWindow+1:end) = true;
for k = 1 : size(traces, 1)-2*halfWindow
    indices(k+(1:movingWindow), k+halfWindow) = true;
end

for roi = 1:size(traces, 2)
    trace = repmat(traces(:,roi), 1, size(traces, 1));
    f_zeros(:,roi) = prctile(reshape(trace(indices), movingWindow, ...
        size(traces,1)), percentile, 1)';
end