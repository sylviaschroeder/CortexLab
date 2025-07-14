
%% Parameters
% following values taken from metadata, explained in 
% https://billkarsh.github.io/SpikeGLX/Support/Metadata_3A.html
samplingRate = 2500; % in ....imec.lf.meta file
numChans = 385;
Imax = 512;
Vmax = 0.6;
lfp_gain = 250;
chanDistance = 0.020; % millimeters

smoothLFP = 21;

baseWin = [-.03 0];

RFtimesInFrames = 0;
lambda = 0.002;
% lambda = 0.01;
crossFolds = 1;

%% Folders
folderData = 'D:\Data';
folderSave = 'D:\Results\SC-LFP\RFs';
subject = 'RT010';
lfpFolders = {'catgt_RT010_2025-04-16_g2\RT010_2025-04-16_g2_imec0', ...
    'catgt_RT010_2025-04-16_g2\RT010_2025-04-16_g2_imec1'};

folderScript = 'C:\dev\workspaces\CortexLab';
folderTools = 'C:\dev\toolboxes';

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')));
addpath(genpath(fullfile(folderTools, 'spikes')));
addpath(genpath(fullfile(folderScript)));

%% Define dataset
% subject = 'V1_insertion2';
% date = '2023-08-21';

%% Load data
colors = colmaps.getBlueWhiteRedMap(255);

f = fullfile(folderData, subject);

% load visual noise data
t_stim = readNPY(fullfile(f, "sparse.startTime.npy"));
stim = readNPY(fullfile(f, "sparse.map.npy")); % [time x rows x columns]
stim = (stim - 0.5) .* 2; % new values: -1 , 0, 1
stimFlat = reshape(stim, size(stim,1), []);

if ~isfolder(fullfile(folderSave, subject))
    mkdir(fullfile(folderSave, subject))
end

for f1 = 1:length(lfpFolders)
    txt = split(lfpFolders{f1}, filesep);
    fSave = fullfile(folderSave, subject, txt{end});

    % load LFP data
    fileLFP = dir(fullfile(f, lfpFolders{f1}, '*.lf.bin'));
    lfp = memmapfile(fullfile(f, lfpFolders{f1}, fileLFP.name), ...
        'Format',  {'int16', [numChans fileLFP.bytes/(numChans*2)], 'x'}); % [channels x time]
    t_lfp = (1:fileLFP.bytes/(numChans*2))' ./ samplingRate;
    % chanLab = readNPY(fullfile(f, "channel_labels.npy"));

    % load stimulus period of LFP
    start = find(t_lfp > t_stim(1)-5, 1);
    stop = find(t_lfp > t_stim(end)+1, 1);
    if isempty(stop)
        ind = find(t_stim < t_lfp(end), 1, "last");
        stim(ind:end,:,:) = [];
        stimFlat(ind:end,:) = [];
        t_stim(ind:end) = [];
        stop = length(t_lfp);
    end
    lfp = double(lfp.Data.x(1:end-1,start:stop)); % last channel is external input -> not LFP
    t_lfp = t_lfp(start:stop);
    % zscore each channel
    lfp = zscore(lfp, 0, 2);
    % smooth lfp in time
    lfp = smoothdata(lfp, 2, "movmean", smoothLFP);

    %% Get stimulus triggered LFP
    stimDur = median(diff(t_stim));
    staSamples = floor(RFtimesInFrames(1) * stimDur * samplingRate) : ...
        ceil((RFtimesInFrames(end)+1) * stimDur * samplingRate);
    baseSamples = floor(baseWin(1)*samplingRate) : ceil(baseWin(2)*samplingRate);
    [~,~,stimTBins] = histcounts(t_stim, t_lfp);

    % lfp_trial: [channels x trial x time]
    % average LFP across time after stimulus onset for each stimulus frame
    % (trial) and each channel, then subtract baseline for each trial
    lfp_trial = reshape( lfp(:, stimTBins + staSamples), ...
        size(lfp,1), length(t_stim), length(staSamples) );
    lfp_trial = lfp_trial - ...
        mean( reshape( lfp(:, stimTBins + baseSamples), ...
        size(lfp,1), length(t_stim), length(baseSamples) ), 3);
    
    mx = max(abs(mean(lfp_trial, 2)), [], "all");
    figure
    imagesc(staSamples([1 end]) ./ samplingRate, [1 size(lfp_trial,1)], ...
        squeeze(mean(lfp_trial, 2)), [-mx mx])
    colormap(colors)
    colorbar
    set(gca,'YDir','normal')
    saveas(gcf, [fSave '_respAfterStimOnset.jpg'])
    close(gcf)

    %% Get receptive fields
    % STA RFs
    rfs = whiteNoise.getReceptiveField(lfp', t_lfp, stim, t_stim, ...
        RFtimesInFrames, lambda, crossFolds); % [rows x columns x t x ON/OFF x chan]
    rfs = squeeze(mean(rfs, 3)); % [rows x columns x ON/OFF x chan]
    m = max(abs(rfs),[],"all");
    % OFF
    figure('WindowState', 'maximized')
    colormap(colors)
    tiledlayout("flow", "TileSpacing", "tight")
    for ch = 1:3:size(lfp,1)
        nexttile
        imagesc(squeeze(-rfs(:,:,2,ch)), [-m m])
        axis image off
        title(sprintf('%d', ch))
    end
    colorbar
    sgtitle(sprintf("%s: OFF fields", subject), "Interpreter", 'none')
    saveas(gcf, [fSave '_OFF_STA.jpg'])
    close(gcf)
    % ON
    figure('WindowState', 'maximized')
    colormap(colmaps.getBlueWhiteRedMap(255))
    tiledlayout("flow", "TileSpacing", "tight")
    for ch = 1:3:size(lfp,1)
        nexttile
        imagesc(squeeze(rfs(:,:,1,ch)), [-m m])
        axis image off
        title(sprintf('%d', ch))
    end
    colorbar
    sgtitle(sprintf("%s: ON fields", subject), "Interpreter", 'none')
    saveas(gcf, [fSave '_ON_STA.jpg'])
    close(gcf)
end
