function [acceptedSpikes, discardedSpikes] = ...
    separateDistributions(st, sa, OT_times, numSTDs, doPlot)

numBins = 40;
span = 20;

if nargin < 5
    doPlot = false;
end

edges = {st(1:100:end), linspace(min(sa),prctile(sa,99.9),numBins)};

binCentres_time = edges{1}(1:end-1) + diff(edges{1})./2;
binCentres_amps = edges{2}(1:end-1) + diff(edges{2})./2;

% bin spikes in time so that each bin contains 100 spikes + bin spikes in
% each time bin by amplitude (using numBins bins); then smooth across
% neighbouring amplitudes for eah time bin
n = hist3([st, sa], 'Edges', edges); % [time x amplitude]; spike counts for each bin
n_sm = zeros(size(n)); % spike counts after smoothing
for t = 1:size(n,1)
    n_sm(t,:) = smooth(n(t,:), span, 'rloess');
end

if doPlot
    figure('Position', [3 694 1916 420])
    imagesc(edges{1}([1 end]), edges{2}([1 end]), n_sm')
    set(gca, 'YDir', 'normal')
    colormap hot
end

% for each time point, find amplitude bin with the highest spike count
n_diff = diff(n_sm,1,2); % 1st derivative for amplitude for each time bin
n_max = diff(sign(n_diff),1,2) == -2; % find peak (where derivative switches from positive to negative)
inds = all(n_max > -2, 2); % exclude time bins with NaNs (I think)
n_max(inds & n(:,1)>n(:,end), 1) = true; % detect when peak is in lowest...
n_max(inds & n(:,1)<n(:,end), end) = true; % ... or highest amplitude bin
[peaks,t] = find(n_max');
t_unique = unique(t); % find time bins where more than one peak was detected
m = hist(t, t_unique);
doubles = find(m > 1);
discard = [];
for d = 1:length(doubles)
    inds = find(t == t_unique(doubles(d)));
    [~,ind] = max(peaks(inds)); % if multiple peaks, take the largest
    disc = setdiff(1:length(inds),ind);
    discard = [discard; inds(disc)];
end
peaks(discard) = [];
t(discard) = [];

% get a smooth trace of most frequent spike amplitude
peaks = binCentres_amps(peaks+1); % get spike amplitude in each peak (spike count) bin
t = binCentres_time(t);
peaks = medfilt1(peaks, 3, 'omitnan', 'truncate');
peaks = interp1(t, peaks, binCentres_time, 'pchip');
peaks = smooth(peaks, 80, 'loess');

if doPlot
    figure('Position', [3 185 1916 420])
    plot(st, sa, 'k.', 'MarkerSize', .1)
    hold on
    plot(binCentres_time, peaks, 'LineWidth', 2)
    axis tight
    xlabel('Time(s)')
    ylabel('Spike amplitude')
end

% exclude bad time periods
if nargin <= 2 || isempty(OT_times)
    inTime = true(size(st));
else
    inTime = false(size(st));
    for j = 1:size(OT_times,1)
        inTime(st >= OT_times(j,1) & st <= OT_times(j,2)) = true;
    end
end


peakPerSpike = interp1(binCentres_time, peaks, st,'pchip'); % get estimated center amplitude for each spike
spikesAbovePeaks = sa(inTime) - peakPerSpike(inTime); % difference of spike amplitude to center amplitude
spikesAbovePeaks(spikesAbovePeaks < 0) = []; % use spike amplitudes above center amplitude to estimate spread of distribution
[counts,bins] = hist(spikesAbovePeaks, 100);
model = fit(bins', counts', @(a,b,x)(a*exp(-(x/b).^2)), ...
    'StartPoint', [counts(1), bins(end)/3]); % use half-Gaussian to fit distribution
sigma = sqrt(model.b / 2);
threshPerSpike = peakPerSpike - numSTDs * sigma; % set lower threshold for acceptable amplitudes to center amplitude minus numSTDs STDs

ind = st > 2100;
threshPerSpike(ind) = 60;

acceptedSpikes = inTime & sa >= threshPerSpike;
discardedSpikes = inTime & sa < threshPerSpike;