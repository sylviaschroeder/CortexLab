%% TODO: Check which direction to subtract when determining 2nd derivative!!!


%% Parameters
% following values taken from metadata, explained in 
% https://billkarsh.github.io/SpikeGLX/Support/Metadata_3A.html
samplingRate = 2500; % in ....imec.lf.meta file
numChans = 385;
Imax = 512;
Vmax = 0.6;
lfp_gain = 250;

chanDistance = 0.020; % millimeters

smoothLFP = {5 21};

stimWin = [-.05 .25];
baseWin = [-.05 0];

colors = colmaps.getBlueWhiteRedMap(255);

%% Folders
folderData = 'D:\Data';
folderSave = 'D:\Results\SC-LFP';

folderScript = 'C:\dev\workspaces\CortexLab';
folderTools = 'C:\dev\toolboxes';

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')));
addpath(genpath(fullfile(folderTools, 'spikes')));
addpath(genpath(fullfile(folderScript)));

%% Define dataset
subject = 'FG005';
date = '2023-08-21';

%% Load data
f = fullfile(folderData, subject, date);

% load visual noise data
t_stim = readNPY(fullfile(f, "sparse.startTime.npy"));
stim = readNPY(fullfile(f, "sparse.map.npy")); % [time x columns x rows]
stim = (stim - 0.5) .* 2; % new values: -1 , 0, 1
stimFlat = reshape(stim, size(stim,1), []);

% load LFP data
fileLFP = dir(fullfile(f, '*.lf.bin'));
lfp = memmapfile(fullfile(f, fileLFP.name), ...
    'Format',  {'int16', [numChans fileLFP.bytes/(numChans*2)], 'x'}); % [channels x time]
t_lfp = (1:fileLFP.bytes/(numChans*2))' ./ samplingRate;
chanLab = readNPY(fullfile(f, "channel_labels.npy"));

% probe channel coordinates
x = repmat([43; 11; 59; 27], 96, 1);
y = reshape(repmat(20:20:3840, 2, 1), [], 1);
channelCoords = [x, y];

% load stimulus period of LFP
start = find(t_lfp > t_stim(1)-5, 1);
stop = find(t_lfp > t_stim(end)+1, 1);
lfp = lfp.Data.x(1:end-1,start:stop); % last channel is external input -> not LFP

% conversion to Volt explained in 
% https://billkarsh.github.io/SpikeGLX/Support/Metadata_3A.html
lfp = double(lfp) .* Vmax  .* 1000 ./ Imax ./ lfp_gain; % mV

t_lfp = t_lfp(start:stop);
% subtract median from each channel
lfp = lfp - median(lfp,2);

% divide channels into 2 columns, replace noisy channel data with 
% interpolation, then smooth across channels and time
noiseChans = chanLab > 0;
chans1 = 1:2:numChans-1;
lfp = lfp(1:2:end,:);
if any(noiseChans(1:2:end))
    [X1, X2] = ndgrid(find(~noiseChans(chans1)), t_lfp);
    F = griddedInterpolant(X1, X2, ...
        lfp(~noiseChans(chans1), :), 'linear');
    [X1, X2] = ndgrid(find(noiseChans(chans1)), t_lfp);
    lfp(noiseChans(chans1),:) = F(X1,X2);
end
lfp = smoothdata2(lfp, "gaussian", smoothLFP);
% lfp2 = smoothdata2(lfp(2:2:end,:), "gaussian", smoothLFP);

%% Get stimulus triggered LFP + most driving pixel
winSamples = stimWin(1)*samplingRate : stimWin(2)*samplingRate;
baseSamples = baseWin(1)*samplingRate : baseWin(2)*samplingRate;
numTrials = max(max(sum(stimFlat > 0,1)), max(sum(stimFlat < 0,1)));
lfpEvokedW = NaN(size(lfp,1), size(stimFlat,2), ...
    length(winSamples)); % [channels x pixels x trials x time]
lfpEvokedB = NaN(size(lfp,1), size(stimFlat,2), ...
    length(winSamples)); % [channels x pixels x trials x time]

[~,~,stimTBins] = histcounts(t_stim, t_lfp);

for pix = 1:size(stimFlat,2)
    % response to white
    ind = stimFlat(:,pix) > 0;
    base = lfp(:,stimTBins(ind) + baseSamples);
    base = reshape(base, size(lfp,1), sum(ind), length(baseSamples)); % [channels x trials x time]
    base = mean(base, 3); % [channels x trials]
    ev = lfp(:, stimTBins(ind) + winSamples);
    ev = reshape(ev, size(lfp,1), sum(ind), length(winSamples)); % [channels x trials x time]
    lfpEvokedW(:,pix,:) = mean(ev - base, 2);

    % response to black
    ind = stimFlat(:,pix) < 0;
    base = lfp(:,stimTBins(ind) + baseSamples);
    base = reshape(base, size(lfp,1), sum(ind), length(baseSamples)); % [channels x trials x time]
    base = mean(base, 3); % [channels x trials]
    ev = lfp(:, stimTBins(ind) + winSamples);
    ev = reshape(ev, size(lfp,1), sum(ind), length(winSamples)); % [channels x trials x time]
    lfpEvokedB(:,pix,:) = mean(ev - base, 2);
end

% find best driving pixel and time of max response
maxW = max(abs(lfpEvokedW), [], "all");
maxB = max(abs(lfpEvokedB), [], "all");
if maxW > maxB
    lfpEvoked = lfpEvokedW;
else
    lfpEvoked = lfpEvokedB;
end
[~,maxInd] = max(abs(lfpEvoked), [], "all");
maxChan = mod(maxInd, size(lfp,1));
maxPix = mod(ceil(maxInd / size(lfp,1)), size(stimFlat,2));
maxT = ceil(maxInd / size(lfp,1) / size(stimFlat,2));

% Find superficial border of SC
profile = lfpEvoked(:,maxPix,maxT);
[minAmp, chanMin] = min(profile);
ind = find(profile(chanMin+1:end) > minAmp/20, 1);
chanTop = ceil(interp1(profile(chanMin:chanMin+ind), chanMin:chanMin+ind, ...
    minAmp/2)) * 2;

% Calculate CSD
if maxW > maxB
    ind = stimFlat(:,maxPix) > 0;
else
    ind = stimFlat(:,maxPix) < 0;
end
base = lfp(:,stimTBins(ind) + baseSamples);
base = reshape(base, size(lfp,1), sum(ind), length(baseSamples)); % [channels x trials x time]
base = mean(base, 3); % [channels x trials]
ev = lfp(:, stimTBins(ind) + winSamples);
ev = reshape(ev, size(lfp,1), sum(ind), length(winSamples)); % [channels x trials x time]
csd = squeeze(mean(diff(diff(ev,1,1),1,1), 2)) ./ (chanDistance^2);
csd = padarray(csd, 1, 0);

[~, chanSource] = max(csd(:,maxT));
chanSink = chanSource - find(csd(chanSource-1:-1:1, maxT) < 0, 1);
chan_SGS_SO = floor(interp1(csd(chanSink:chanSource, maxT), ...
    chanSink:chanSource, 0)) * 2;

% Make plots
figure('Position', [410 50 820 946])
tiledlayout(1,2)
nexttile
m = max(abs(profile),[],"all");
imagesc(stimWin, [1 numChans-1], squeeze(lfpEvoked(:,maxPix,:)), [-m m])
colormap(colors)
c = colorbar;
c.Label.String ='mV';
hold on
plot(stimWin, [1 1].*chanTop, 'k', 'LineWidth', 2)
set(gca, 'YDir', 'normal')
xlabel('Time (ms)')
ylabel('#Channel')
title('Evoked LFP')
legend(sprintf('Top of SC (chan: %d)', chanTop))

nexttile
m = max(abs(csd),[],"all");
imagesc(stimWin, [1 numChans-1], csd, [-m m])
colormap(colors)
c = colorbar;
c.Label.String = 'mV/mm^2';
hold on
plot(stimWin, [1 1].*chan_SGS_SO, 'k:', 'LineWidth', 2)
set(gca, 'YDir', 'normal')
xlabel('Time (ms)')
ylabel('#Channel')
title('Current-source density')
legend(sprintf('SGS-SO border (chan: %d)', chan_SGS_SO))