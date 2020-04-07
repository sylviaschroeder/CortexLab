%% Load data for all units in recording
db_ephys_opticTract
k = 3;
probe = 1;

serverFolder = 'J:\Ephys';
dataFolder = 'J:\Ephys';

fprintf('\nProcessing: %s %s\n', db(k).subject, db(k).date);

sp = loadAllKsDir(db(k).subject, db(k).date);
b = [1 0];
if db(k).OTprobe > 1
    alignDir = fullfile(serverFolder, db(k).subject, db(k).date, 'alignments');
    b = readNPY(fullfile(alignDir, ...
        sprintf('correct_ephys_%s_to_ephys_%s.npy', sp(probe).name, db(k).probeNames{1})));
end

%% Separate good from bad spikes of specific unit
unit = 14; %21; %34;
samplingRate = 30000;
cols = lines(2);

st = sp(probe).st(sp(probe).clu == unit);
sa = sp(probe).spikeAmps(sp(probe).clu == unit);

% Plot all spike amplitudes
figure
plot(st, sa, 'k.', 'MarkerSize', .1)
xlabel('Time (s)')
ylabel('Spike amplitude')
axis tight
set(gca, 'box', 'off')

OT_ID = find(db(k).OTunits{probe} == unit);
[accepted, discarded] = general.separateDistributions(st, sa, ...
    db(k).OTtimes{probe}{OT_ID}, db(k).OTampSTDs{probe}(OT_ID), true);
% [accepted, discarded] = general.separateDistributions(st, sa, ...
%     [15 600], 10);
% [accepted, discarded] = general.separateDistributions(st, sa, ...
%     [0 570;770 2600;4820 5500], 10);
ind = st > st(end)-.002;
accepted(ind) = false;
discarded(ind) = false;

% Plot amplitude of each spike in time
figure('Position',[4 815 1350 300])
hold on
plot(st(accepted), sa(accepted), '.', 'Color', cols(1,:), 'MarkerSize', .1)
plot(st(discarded), sa(discarded), '.', 'Color', cols(2,:), 'MarkerSize', .1)
axis tight
set(gca, 'box', 'off')
xlabel('Time (s)')
ylabel('Spike amplitude')

% Plot wave shape of good vs bad spikes
params.dataDir = fullfile(dataFolder, db(k).subject, db(k).date, ...
    ['ephys_' sp(probe).name]);
params.fileName = sprintf('%s_%s_%s_g0_t0.imec.ap_CAR.bin', ...
    db(k).subject, db(k).date, sp(probe).name);
params.dataType = 'int16';
params.nCh = 385;
params.wfWin = round([-.0005 .001] .* samplingRate);
params.nWf = 40000;
% for spikeTimes: convert spike times (in time of probe 1) to time of
% selected probe, then multiply by sampling rate
params.spikeTimes = round(([st(accepted); st(discarded)] - b(2)) ./ ...
    b(1) .* samplingRate);
params.spikeClusters = [ones(sum(accepted),1); ones(sum(discarded),1).*2];

wf = getWaveForms(params);
figure('Position',[1359 42 560 1074])
hold on
for w = 1:2
    thisTemp = squeeze(wf.waveFormsMean(w,:,:))';
    if w == 1
        [~,peakChan] = max(max(abs(thisTemp),[],1),[],2);
    end
    plotChans = max(1,peakChan-5) : min(length(sp(probe).xcoords),peakChan+5);
    arrayfun(@(x)plot(sp(probe).xcoords(x)+0.6*(1:size(thisTemp,1))', ...
        sp(probe).ycoords(x)+thisTemp(:,x), 'Color', cols(w,:)), plotChans)
end
ycoords = unique(sp(probe).ycoords(plotChans));
mini = min(ycoords);
maxi = max(ycoords);
chanDist = median(diff(ycoords));
ylim([mini - 0.5*chanDist, maxi + 0.5 * chanDist])
ylabel('Depth (um)')
set(gca, 'XTick', [])

numWF = 200;
wfs = squeeze(wf.waveForms(:,:,peakChan,:)); % [neurons x spikes x time]
maxi = -Inf;
mini = Inf;
f = [0 0];
for w = 1:2
    f(w) = figure('Position',[4 815-w*385 1350 300]);
    waves = squeeze(wfs(w,:,:));
    waves(all(isnan(waves),2),:) = [];
    freq = max(1, floor(size(waves,1) / numWF));
    waves = waves((1 : min(numWF,size(waves,1))) .* freq, :);
    maxi = max(maxi,max(waves(:)));
    mini = min(mini,min(waves(:)));
    
    plot((params.wfWin(1):params.wfWin(end))./samplingRate,waves', 'Color', cols(w,:), 'LineWidth', .1)
    set(gca, 'box', 'off')
end
for w = 1:2
    figure(f(w))
    xlim(params.wfWin([1 end])./samplingRate)
    ylim([mini maxi])
end

%% Overwrite data files (spike_clusters, cluster_groups)
clu = readNPY(fullfile(dataFolder, db(k).subject, db(k).date, ...
    ['ephys_' sp(probe).name], 'spike_clusters.npy'));
spInd = find(clu==unit);
clu(spInd(discarded)) = max(clu)+1;
writeNPY(clu, fullfile(dataFolder, db(k).subject, db(k).date, ...
    ['ephys_' sp(probe).name], 'spike_clusters.npy'))