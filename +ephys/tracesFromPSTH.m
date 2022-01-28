function [traces, bins] = tracesFromPSTH(times, trials, groups, limits, ...
    binSize, sigma, isCausal)

% Calculate smoother PSTHs from spike times.
%
% traces    [group x t]; One trace per group in spikes per second
% bins      [1 x t]; time of traces
%
% times     [st x 1]; spike times aligned to event of interest
% trials    [st x 1]; trial ID of spike
% groups    [trial x 1]; group ID of each trial
% limits    [pre post]; time spans before and after each aligned event,
%           which defines the length of the traces
% binSize   binSize
% sigma     sigma of smoothing Gaussian
% isCausal  logical; if true, smoothing kernel is causal

if nargin < 7
    isCausal = true;
end

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
if isCausal
    win(x < 0) = [];
end
win = win ./ sum(win);
w = ones(1,length(win));
for g = uniGr
    if isCausal
        tr = conv([w .* nanmean(traces(g,1:sigBins)), ...
            traces(g,:)], win);
    else
        tr = conv([w .* nanmean(traces(g,1:sigBins)), ...
            traces(g,:), w .* nanmean(traces(g,end-sigBins:end))], win, 'same');
    end
    tr = tr(length(win)+(1:length(bins)));
    traces(g,:) = tr;
end