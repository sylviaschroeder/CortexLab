function [traces, bins] = tracesFromPSTH(times, trials, groups, limits, ...
    binSize, sigma)

if prod(limits) < 0
    edges = [-flip(0:binSize:-limits(1)) binSize:binSize:limits(2)];
else
    edges = limits(1):binSize:limts(2);
end
bins = edges(1:end-1) + 0.5*binSize;

% bin spikes of each trial
binned = zeros(length(groups), length(bins));
for tr = 1:length(groups)
    ind = trials == tr;
    if sum(ind) < 1
        continue
    end
    binned(tr,:) = histcounts(times(ind), edges);
end

% calculate spikes per sec
binned = binned ./ binSize;

% average trials of same group
uniGr = unique(groups(~isnan(groups)))';
traces = NaN(max(uniGr), length(bins));
for g = uniGr
    tr = groups == g;
    if sum(tr) < 1
        continue
    end
    traces(g,:) = nanmean(binned(tr,:), 1);
end

% smooth traces with causal half-gaussian
sigBins = sigma/binSize;
x = floor(-5*sigBins):ceil(5*sigBins);
win = normpdf(x, 0, sigBins);
win(x < 0) = [];
win = win ./ sum(win);
for g = uniGr
    tr = conv([ones(1,length(win)) .* nanmean(traces(g,1:sigBins)), ...
        traces(g,:)], win);
    traces(g,:) = tr(length(win)+(1:length(bins)));
end