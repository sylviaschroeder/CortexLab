%% Defintion of dataset
subject = 'SS082';
date = '2017-11-25';
ephysTag = 'K1';
ephysMaster = 'K1';
exp = 2;
TLexp = 1;
samplingRate = 30000;

%% Folders
ephysFolder = 'J:\Ephys\';
subjectsFolder = '\\zubjects.cortexlab.net\Subjects\';
expInfoFolder = '\\ZSERVER.cortexlab.net\Data\expInfo';
stimulusFolder = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\experimentCode';

alignFolder = fullfile(subjectsFolder, subject, date, 'alignments');

%% Load sync information
if ~strcmp(ephysTag, ephysMaster)
    bEphysToMaster = readNPY(fullfile(alignFolder, ...
        sprintf('correct_ephys_%s_to_ephys_%s.npy', ephysTag, ephysMaster)));
else % this one is master, so use a dummy conversion
    bEphysToMaster = [1; 0];
end

bTLtoMaster = readNPY(fullfile(alignFolder, ...
    sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, ephysMaster)));
stimTimesTL = readNPY(fullfile(alignFolder, ...
    sprintf('mpep_%d_onsets_in_timeline_%d.npy', exp, TLexp)));
% stimTimes = readNPY(fullfile(alignFolder, ...
%     sprintf('correct_block_%d_to_timeline_%d.npy', exp, TLexp)));

%% Get stimulus data
% Load stimulus parameters
data = load(fullfile(expInfoFolder, subject, date, num2str(exp), ...
    sprintf('%s_%d_%s_parameters.mat', date, exp, subject)));
parsNoise = data.parameters.Protocol;
stimFile = str2func(strtok(parsNoise.xfile, '.'));
% load myScreenInfo
load(fullfile(expInfoFolder, subject, date,num2str(exp), ...
    sprintf('%s_%d_%s_hardwareInfo.mat', date, exp, subject)));
myScreenInfo.windowPtr = NaN;
SS = stimFile(myScreenInfo, parsNoise.pars);
stimFrames = cat(3, SS.ImageTextures{:});

noiseOn = applyCorrection(stimTimesTL, bTLtoMaster);
framesPerImage = parsNoise.pars(6,1);
frameTimesRelative = (0 : size(stimFrames, 3)-1) * framesPerImage / myScreenInfo.FrameRate;
frameTimes = (frameTimesRelative + noiseOn)';

%% Parameters for raw ephys data
parts = 1;
partitioning = [round(numel(frameTimes)/parts) .* (1:parts-1), numel(frameTimes)];
frameIndices = cell(1,parts);
k = 1;
for p = 1:parts
    frameIndices{p} = k:partitioning(p);
    k = partitioning(p)+1;
end
params.dataDir = fullfile(ephysFolder, subject, date);
params.fileName = sprintf('%s_%s_K1_g0_t0.imec.ap_CAR.bin', subject, date);
params.dataType = 'int16';
params.nCh = 385;
params.wfWin = round([0 .06] .* samplingRate);
params.nWf = length(frameIndices{1});
% for spikeTimes: first convert stimulus times (in master time) to time of
% probe, then multiply by sampling rate
params.spikeTimes = round((frameTimes(frameIndices{1}) - bEphysToMaster(2)) ./ ...
    bEphysToMaster(1) .* samplingRate);
params.spikeClusters = ones(length(frameIndices{1}),1);

wf = getWaveForms(params);

%% Plot all mean traces
meanStd = mean(std(wf.waveFormsMean(1,:,:),0,3),2);
t = ((1:size(wf.waveForms,4)) + params.wfWin(1)) ./ samplingRate;

figure('Position', [1 41 1920 1083])
plot(t,squeeze(wf.waveFormsMean)'./(7*meanStd)+(0:size(wf.waveFormsMean,2)-1), 'k')
axis tight
set(gca, 'box', 'off')

% first take absolute value around mean before averaging across trials
meanAbsResp = NaN(size(wf.waveForms,3), size(wf.waveForms,4));
for chan = 1:size(wf.waveForms,3)
    r = abs(squeeze(wf.waveForms(1,:,chan,:)) - mean(wf.waveFormsMean(1,chan,:),3));
    meanAbsResp(chan,:) = mean(r);
end
meanStdAbs = mean(std(meanAbsResp,0,2));
figure('Position', [1 41 1920 1083])
plot(t, meanAbsResp'./(7*meanStdAbs)+(0:size(meanAbsResp,1)-1), 'k')
axis tight
set(gca, 'box', 'off')

%% Plot all traces of single channel
ch = 50;
reps = 1:500;
chMap = readNPY(fullfile(params.dataDir, 'channel_map.npy'))+1;
chInd = find(chMap >= ch, 1);
meanStd = mean(reshape(std(wf.waveForms(:,:,chInd,:),0,4),[],1));
t = ((1:size(wf.waveForms,4)) + params.wfWin(1)) ./ samplingRate;

figure('Position', [1 41 1920 1083])
plot(t, squeeze(wf.waveForms(1,reps,chInd,:))' + (0:length(reps)-1) .*(7*meanStd), 'k')
axis tight
set(gca, 'box', 'off')

for stim = 1:length(contrastPars)
%     subplot(1, length(contrastPars), stim)
    figure('Position', [1 41 1920 1083])
    hold on
    k = 0;
    for las = 1:length(laserPars)
        stimID = (stim-1)*length(laserPars) + las;
        for rep = 1:reps
            plot(t, squeeze(wf.waveForms(stimID,rep,chInd,:))-7*k*meanStd, 'k')
            k = k + 1;
        end
        k = k + 5;
    end
    axis tight
    set(gca, 'YTick', (-(reps+5)*(length(laserPars)-1)-reps/2 : (reps+5) : -reps/2).*7.*meanStd, ...
        'YTickLabel', flip(laserPars))
    xlabel('Time from stimulus onset (s)')
    ylabel('Laser power (V) (max. 5V)')
    title(sprintf('Contrast: %d%%; Channel: %d', round(contrastPars(stim)*100), ch))
end
h = findobj('type', 'figure');
savefig(h, sprintf('channel_%03d.fig', ch))

%% Plot PSTHs
% sp = loadAllKsDir(subject, date);
probe = 2;
depth = [2370 3240];
% minAmp = 50;

% validClu = find(sp(probe).templateYpos >= depth(1) & ...
%     sp(probe).templateYpos <= depth(2) & sp(probe).tempAmps > minAmp);
validClu = find(sp(probe).templateYpos >= depth(1) & ...
    sp(probe).templateYpos <= depth(2));

goodSpikes = ismember(sp(probe).clu, validClu-1);

relevantStims = length(contrastPars)*length(laserPars) + (-length(laserPars)+1:0);
trialInds = ismember(stimSeq, relevantStims);

% psthViewer(applyCorrection(sp(probe).st(goodSpikes), bEphysToMaster), ...
%     sp(probe).clu(goodSpikes), stimOn(trialInds), [-.1 stimDur+.1], stimSeq(trialInds))
psthViewer(sp(probe).st(goodSpikes), ...
    sp(probe).clu(goodSpikes), stimOn(trialInds), [-.1 stimDur+.1], stimSeq(trialInds))