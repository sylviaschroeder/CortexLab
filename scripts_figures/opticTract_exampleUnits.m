db_ephys_opticTract

%% Example units
dataset = 3; %(SS096, 2018-03-08)
probe = 1;
unit = 70;

dataset = 5; %(SS098, 2018-03-16)
probe = 1;
unit = 65;

samplingRate = 30000;

%% Parameters
runThreshold = 1;
runningSigma = 0.25;
sigma = 1;
binSizeRun = 1/7.5;
minRunTime = 3;

binEdges = db(dataset).darkTime(1) : binSizeRun : db(dataset).darkTime(2);
timeBins = binEdges(1:end-1) + binSizeRun/2;

sig = round(sigma / binSizeRun);
win = normpdf(-5*sig : 5*sig, 0, sig);
highPassWindow = 180; % in sec; to high-pass filter firing rates
prctileFilter = 8;

binSizeFlicker = 0.001;

%% Folders
ephysFolder = 'J:\Ephys\';
plotFolder = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\OpticTract\plots_runningCorr\spikeWaveforms_running-notRunning';
subjectsFolder = '\\ZUBJECTS.cortexlab.net\Subjects';

%% Load spike times
sp = loadAllKsDir(db(dataset).subject, db(dataset).date);
st = sp(probe).st(sp(probe).clu == unit);

%% Load sync information
[expNums, ~, ~, ~, ~, tl, hasTimeline] = ...
    dat.whichExpNums(db(dataset).subject, db(dataset).date);
Tlexp = find(hasTimeline, 1, 'last');
tl = tl{end};
alignDir = fullfile(ephysFolder, db(dataset).subject, db(dataset).date, 'alignments');
if ~strcmp(sp(probe).name, db(dataset).probeNames{1})
    bEphysToMaster = readNPY(fullfile(alignDir, ...
        sprintf('correct_ephys_%s_to_ephys_%s.npy', sp(probe).name, ...
        db(dataset).probeNames{1})));
else % this one is master, so use a dummy conversion
    bEphysToMaster = [1; 0];
end
bTLtoMaster = readNPY(fullfile(alignDir, ...
    sprintf('correct_timeline_%d_to_ephys_%s.npy', Tlexp, db(dataset).probeNames{1})));

%% Load running data
tlTime = applyCorrection(tl.rawDAQTimestamps, bTLtoMaster);

% Load and prepare data for running correlation
rotary = double(tl.rawDAQData(:, strcmp({tl.hw.inputs.name}, ...
    'rotaryEncoder')));
runSpeed = nonVis.getRunningSpeed_wheel(rotary, tlTime, runningSigma);
runningTime = runSpeed.t;
runSpeed = runSpeed.total;
cmPerUnit = 2*pi * 8.75 / (4 * 1024);
runSpeed = runSpeed * cmPerUnit;
run = interp1(runningTime, runSpeed, timeBins, 'pchip'); % run speed during darkness only
run = [ones(1,length(win)) .* mean(run(1:min(length(run),length(win)))), ...
    run, ones(1,length(win)) .* mean(run(end-min(length(run),length(win))+1:end))];
run = conv(run, win, 'same');
run = run(length(win)+1 : end-length(win));
isRunning = run > runThreshold;

sta = find(diff(isRunning)==1);
sto = find(diff(isRunning)==-1);
if sta(1)>sto(1)
    sta = [1, sta];
end
if sta(end)>sto(end)
    sto(end+1) = length(beh);
end
ind = find(sta(2:end) - sto(1:end-1) < minRunTime/binSizeRun);
sta(ind+1) = [];
sto(ind) = [];
ind = (sto - sta) >= minRunTime/binSizeRun;
sta = sta(ind);
sto = sto(ind);

%% Plot firing rate and running speed and mark times
iCell = db(dataset).OTunits{probe} == unit;
validTime = db(dataset).OTtimes{probe}{iCell};
st_dark = st;
st_dark(st_dark<binEdges(1) | st_dark>binEdges(end)) = [];
[sr,~,b] = histcounts(st_dark, binEdges);
sr = [ones(1,length(win)) .* mean(sr(1:min(length(sr),length(win)))), ...
    sr, ones(1,length(win)) .* ...
    mean(sr(end-min(length(sr),length(win))+1:end))] ./ binSizeRun;
sr = conv(sr, win, 'same');
sr = sr(length(win)+1 : end-length(win));
[sr, smoothed] = preproc.removeSlowDrift(sr', 1/binSizeRun, ...
    highPassWindow, prctileFilter);
sr = sr' + mean(smoothed);
if ~isempty(validTime)
    sr2 = NaN(size(sr));
    run2 = NaN(size(run));
    for k = 1:size(validTime,1)
        ind = timeBins>=validTime(k,1) & timeBins<=validTime(k,2);
        sr2(ind) = sr(ind);
        run2(ind) = run(ind);
    end
    sr = sr2;
    run = run2;
end

figure('Position', [4 678 1914 420]);
subplot(2,1,1)
hold on
subplot(2,1,2)
hold on
for k = 1:length(sta)
    subplot(2,1,1)
    fill(timeBins([sta(k) sto(k) sto(k) sta(k)]), [[1 1].*min(run), [1 1].*max(run)], 'k', ...
        'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.1)
    subplot(2,1,2)
    fill(timeBins([sta(k) sto(k) sto(k) sta(k)]), [[1 1].*min(sr), [1 1].*max(sr)], 'k', ...
        'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.1)
end
ax = [0 0];
subplot(2,1,1)
plot(timeBins, run, 'k')
ylabel('Running speed (cm/s)')
title(sprintf('%s %s unit %03d', db(dataset).subject, ...
    db(dataset).date, unit), 'Interpreter', 'none')
axis tight
ax(1) = gca;
subplot(2,1,2)
plot(timeBins, sr, 'k')
xlabel('Time (s)')
ylabel('Firing rate (sp/s)')
axis tight
ax(2) = gca;
linkaxes(ax, 'x')
xlim(timeBins([1 end]))

% savefig(figs, fullfile(plotFolder, sprintf('%s_%s_unit%03d.fig', ...
%     db(dataset).subject, db(dataset).date, unit)))

%% Load stimulus information: flickering screens
data = load(fullfile(subjectsFolder, db(dataset).subject, db(dataset).date, ...
    num2str(db(dataset).expFlicker), sprintf('%s_%d_%s_parameters.mat', ...
    db(dataset).date, db(dataset).expFlicker, db(dataset).subject)));
parsFlicker = data.parameters.Protocol;
data = load(fullfile(subjectsFolder, db(dataset).subject, db(dataset).date, ...
    num2str(db(dataset).expFlicker), sprintf('%s_%d_%s_hardwareInfo.mat', ...
    db(dataset).date, db(dataset).expFlicker, db(dataset).subject)));
monitorRefreshRate = data.myScreenInfo.FrameRate;
flickerDurs = parsFlicker.pars(strcmp(parsFlicker.parnames, ...
    'nfr'),:) / monitorRefreshRate;

stimOnTL = readNPY(fullfile(alignDir, ...
    sprintf('mpep_%d_onsets_in_timeline_%d.npy', db(dataset).expFlicker, Tlexp))); % times of stim on- and offsets
% stimOffTL = readNPY(fullfile(alignDir, ...
%     sprintf('mpep_%d_offsets_in_timeline_%d.npy', db(dataset).expFlicker, Tlexp)));
% flickerWhite = applyCorrection(stimOnTL, bTLtoMaster);
% flickerBlack = applyCorrection(stimOffTL, bTLtoMaster);
% flickerAll = reshape([flickerWhite, flickerBlack]', [], 1);
flickerTimes = cell(size(parsFlicker.seqnums));
flickerAll = applyCorrection(stimOnTL, bTLtoMaster);

longestFlicker = max(flickerDurs);
durs = diff(flickerAll);
flickerStimStartInds = [1; find(durs > ...
    longestFlicker * 1.5) + 1; length(flickerAll)+1];
for stim = 1:size(flickerTimes,1)
    for trial = 1:size(flickerTimes,2)
        ind = parsFlicker.seqnums(stim,trial);
        fl = flickerAll(flickerStimStartInds(ind) : ...
            flickerStimStartInds(ind+1)-1);
        fl = fl(1:floor(length(fl)/2)*2);
        flickerTimes{stim,trial} = reshape(fl, 2, [])';
    end
end

%% Plot response to flickering screens
psth_f = cell(1, size(flickerTimes,2));
stim = 1;
lum = 1; % white
limits = [-.1 .2];
onToOn = [];
onToOff = [];
for rep = 1:size(flickerTimes,2)
    flicks = flickerTimes{stim,rep};
    [psth_f{1,rep}, bins_f] = ...
        psthAndBA(st, flicks(:,lum), 2.*limits, binSizeFlicker);
    onToOn = [onToOn; diff(flicks(:,1))];
    onToOff = [onToOff; diff(flicks,1,2)];
end

bins_f = bins_f - 0.5*binSizeFlicker;
onToOn = median(onToOn);
onToOff = median(onToOff);
flicks = (-5:10)'.*onToOn;
flicks = [flicks, flicks + onToOff];
stim = [ones(size(flicks,1),1), zeros(size(flicks,1),1)];
ax = [0 0];
figure
subplot(4,1,1)
stairs(reshape(flicks',[],1), reshape(stim',[],1), 'k')
ax(1) = gca;
ylabel('Luminance')
title(sprintf('Flickering monitor (%s %s unit %d)', db(dataset).subject, ...
    db(dataset).date, unit))
subplot(4,1,2:4)

stairs(bins_f, mean(cat(1,psth_f{1,:}),1), 'k')
% hold on
% y = 0;
% for rep = 1:size(psth_f,2)
%     stairs(bins_f, psth_f{1,rep} - y - 1.1*max(psth_f{1,rep}), 'k')
%     y = y + 1.1*max(psth_f{1,rep});
% end

ax(2) = gca;
linkaxes(ax, 'x')
xlim([0 4].*onToOff-0.02)
set(ax(1), 'box', 'off')
set(ax(2), 'box', 'off')
xlabel('Time (s)')
ylabel('Firing rate (sp/s)')

%% Load stimulus information: visual noise
data = load(fullfile(subjectsFolder, db(dataset).subject, db(dataset).date, ...
    num2str(db(dataset).expNoise), sprintf('%s_%d_%s_parameters.mat', ...
    db(dataset).date, db(dataset).expNoise, db(dataset).subject)));
parsNoise = data.parameters.Protocol;
stimFile = str2func(strtok(parsNoise.xfile, '.'));
% load myScreenInfo
load(fullfile(subjectsFolder, db(dataset).subject, db(dataset).date, ...
    num2str(db(dataset).expNoise), sprintf('%s_%d_%s_hardwareInfo.mat', ...
    db(dataset).date, db(dataset).expNoise, db(dataset).subject)));
myScreenInfo.windowPtr = NaN;

stimOnTL = readNPY(fullfile(alignDir, ...
    sprintf('mpep_%d_onsets_in_timeline_%d.npy', db(dataset).expNoise, Tlexp)));
noiseOn = applyCorrection(stimOnTL, bTLtoMaster);
% call x-file to create stimuli
SS = stimFile(myScreenInfo, parsNoise.pars);
stimFrames = cat(3, SS.ImageTextures{:});

framesPerImage = parsNoise.pars(6,1);
frameTimes = (0 : size(stimFrames, 3)-1) * framesPerImage / myScreenInfo.FrameRate;
noiseFrameTimes = (frameTimes + noiseOn)';
noisePosition = parsNoise.pars(2:5)./10;

%% Plot receptive field
respPerFrame = histcounts(st, [noiseFrameTimes(:); ...
    noiseFrameTimes(end)+ framesPerImage / myScreenInfo.FrameRate]) / ...
    (framesPerImage / myScreenInfo.FrameRate);
respPerFrame = respPerFrame - mean(respPerFrame);
stims = reshape(stimFrames, [], size(stimFrames,3))';
stims = repmat(stims, size(noiseFrameTimes,2), 1);
kernels = zeros(size(stims,2), 6);
for delay = 0:5
    kernels(:,delay+1) = abs(stims(1:end-delay,:)) \ respPerFrame(1+delay:end)';
%     kernels(:,delay+1) = stims(1:end-delay,:) \ respPerFrame(1+delay:end)';
end
[RFmaxi,ind] = max(max(abs(kernels),[],1));
kernels = kernels(:,ind);
pix = reshape(kernels, size(stimFrames,1), size(stimFrames,2));

grad = linspace(0,1,40)';
reds = [ones(40,1),grad,grad];
blues = [grad,grad,ones(40,1)];
cm = [blues; flip(reds(1:end-1,:),1)];
stdResp = std(respPerFrame);
figure
imagesc(noisePosition([1 2]), noisePosition([3 4]), ...
    pix, [-2 2].*stdResp)
% imagesc(noisePosition([1 2]), noisePosition([3 4]), ...
%     pix, [-1 1].*RFmaxi)
c = colorbar;
c.Label.String = 'spikes/s';
title('Visual noise (absolute RF)')
xlabel('Azimuth (vis. deg.)')
ylabel('Altitude (vis. deg.)')
colormap(gca,cm)
axis image
set(gca, 'box', 'off')

%% Parameters for raw ephys data
params.dataDir = fullfile(ephysFolder, db(dataset).subject, ...
    db(dataset).date, ['ephys_' sp(probe).name]);
params.fileName = sprintf('%s_%s_%s_g0_t0.imec.ap_CAR.bin', ...
    db(dataset).subject, db(dataset).date, sp(probe).name);
params.dataType = 'int16';
params.nCh = 385;
params.wfWin = round([-0.001 .002] .* samplingRate);
params.nWf = 10000;
% for spikeTimes: first convert stimulus times (in master time) to time of
% probe, then multiply by sampling rate
params.spikeTimes = round((st - bEphysToMaster(2)) ./ ...
    bEphysToMaster(1) .* samplingRate);
clusters = isRunning(b)'+1;
params.spikeClusters = clusters;

wf = getWaveForms(params);

%% Plot waveforms
numChans = 10;
cols = 'kr';

% find channel with largest spike
[~, peakChan] = max(max(abs(mean(wf.waveFormsMean,1)),[],3),[],2);
plotChans = max(1,peakChan-numChans) : min(length(sp(probe).xcoords),peakChan+numChans);
ycoords = unique(sp(probe).ycoords(plotChans));
mini = min(ycoords);
maxi = max(ycoords);
chanDist = median(diff(ycoords));

figure('Position', [1145 42 754 1074]);
hold on
h = [0 0];
for c = 1:2
    H = arrayfun(@(x)plot(sp(probe).xcoords(x)+0.3*(1:size(wf.waveFormsMean,3))', ...
        sp(probe).ycoords(x)+.5.*squeeze(wf.waveFormsMean(c,x,:)), cols(c)), plotChans);
    h(c) = H(1);
end
ylim([mini - 0.5*chanDist, maxi + 1.5 * chanDist])
ylabel('Depth (um)')
set(gca, 'XTick', [])
legend(h, {'stationary','running'}, 'Location', 'NorthEast')
title(sprintf('%s %s %s, unit %d', db(dataset).subject, ...
    db(dataset).date, sp(probe).name, unit))