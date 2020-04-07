%% Defintion of dataset
subject = 'SS090';
% date = '2017-12-19';
% probe = 1;
% exp = 7;
% TLexp = 1;
% date = '2018-01-10';
% probe = 2;
% exp = 5;
% TLexp = 2;
date = '2018-01-31';
probe = 1;
exp = 11;
TLexp = 6;

% date = '2018-01-31';
% tags and hemispheres
% SS090, 2017-12-19: K1 - left, K2 - right
% SS090, 2018-01-10: K1 - right, K2 - left
% SS090, 2018-01-31: K1 - left, K2 - right
ephysTag = 'K1';
ephysMaster = 'K1';
% experiments where visual stimulus and laser come on simultaneously
% SS090, 2018-01-10: left hemisphere: 5, right hemisphere: 6, TL: 2
% SS090, 2018-01-31: left: 11, right: 12, TL: 6
samplingRate = 30000;

%% Folders
ephysFolder = 'J:\Ephys\';
subjectsFolder = '\\zubjects.cortexlab.net\Subjects\';

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

%% Get stimulus data
% Load stimulus parameters
data = load(fullfile(subjectsFolder, subject, date, num2str(exp), ...
    sprintf('%s_%d_%s_Block.mat', date, exp, subject)));
block = data.block;
laser = [block.paramsValues.laserAmp]';
contrast = [block.paramsValues.stimulusContrast]';
contrast = contrast(:,1); % assume that contrast is always the same on right and left side
stimDur = unique([block.paramsValues.flashDur]);
stimDelay = unique([block.paramsValues.visOnsetDelay]);
laserDur = unique([block.paramsValues.laserDuration]);
laserDelay = unique([block.paramsValues.laserOnsetDelay]);

laserPars = unique(laser);
contrastPars = unique(contrast);
stimIDs = reshape(1 : length(laserPars)*length(contrastPars), ...
    length(laserPars), length(contrastPars));
stimSeq = NaN(length(laser),1);
k = 1;
reps = 0;
for c = 1:length(contrastPars)
    for l = 1:length(laserPars)
        ind = laser == laserPars(l) & contrast == contrastPars(c);
        stimSeq(ind) = k;
        reps = max(reps, sum(ind));
        k = k + 1;
    end
end

% -- use actual event times recorded with timeline
stimTimesTL = readNPY(fullfile(alignFolder, ...
    sprintf('block_%d_sw_in_timeline_%d.npy', exp, TLexp)));
% -- or use fit mapping block/mpep times to timeline
% b = readNPY(fullfile(alignFolder, ...
%     sprintf('correct_block_%d_to_timeline_%d.npy', exp, TLexp)));
% fileBlock = dat.expFilePath(subject, date, exp, 'block', 'master');
% data = load(fileBlock);
% sw = data.block.stimWindowUpdateTimes;
% stimTimesTL = applyCorrection(sw,b);
stimOn = applyCorrection(stimTimesTL(1:2:end), bTLtoMaster);
stimOff = applyCorrection(stimTimesTL(2:2:end), bTLtoMaster);

%% Parameters for raw ephys data
% params.dataDir = fullfile(ephysFolder, subject, date, ephysTag);
% params.fileName = sprintf('%s_%s_%s_g1_t0.imec.ap_CAR.bin', subject, date, ephysTag);
% params.dataType = 'int16';
% params.nCh = 385;
% % params.wfWin = [-.1 .6] .* samplingRate;
% params.wfWin = [1.2 1.9] .* samplingRate;
% params.nWf = reps;
% % for spikeTimes: first convert stimulus times (in master time) to time of
% % probe, then multiply by sampling rate
% params.spikeTimes = round((stimOn - bEphysToMaster(2)) ./ ...
%     bEphysToMaster(1) .* samplingRate);
% params.spikeClusters = stimSeq;
% 
% wf = getWaveForms(params);

%% Plot raw data
% % SS090, 2018-01-10, K2 (left): upper SC -> channels 305-324
% ch = 265;
% chMap = readNPY(fullfile(params.dataDir, 'channel_map.npy'))+1;
% chInd = find(chMap >= ch, 1);
% meanStd = mean(reshape(std(wf.waveForms(:,:,chInd,:),0,4),[],1));
% t = ((1:size(wf.waveForms,4)) + params.wfWin(1)) ./ samplingRate;
% indStart = 1;
% indEnd = length(t);
% % indStart = find(t > -.02, 1);
% % indEnd = find(t > 0.12, 1);
% 
% for stim = 1:length(contrastPars)
% %     subplot(1, length(contrastPars), stim)
%     figure('Position', [1 41 1920 1083])
%     hold on
%     k = 0;
%     for las = 1:length(laserPars)
%         stimID = (stim-1)*length(laserPars) + las;
%         for rep = 1:reps
%             plot(t(indStart:indEnd), squeeze(wf.waveForms(stimID,rep,chInd,indStart:indEnd))-7*k*meanStd, 'k')
%             k = k + 1;
%         end
%         k = k + 5;
%     end
%     axis tight
%     set(gca, 'YTick', (-(reps+5)*(length(laserPars)-1)-reps/2 : (reps+5) : -reps/2).*7.*meanStd, ...
%         'YTickLabel', flip(laserPars))
%     xlabel('Time from stimulus onset (s)')
%     ylabel('Laser power (V) (max. 5V)')
%     title(sprintf('Contrast: %d%%; Channel: %d', round(contrastPars(stim)*100), ch))
% end
% h = findobj('type', 'figure');
% % savefig(h, sprintf('channel_%03d.fig', ch))
% savefig(h, sprintf('channel_%03d_laserOff.fig', ch))

%% Plot PSTHs
sp = loadAllKsDir(subject, date);
% depth = [2750 3250]; % sSC in SS090, 2018-01-10, K2 (left)
% depth = [2200 2750]; % dSC in SS090, 2018-01-10, K2 (left)
% minAmp = 50;

% validClu = find(sp(probe).templateYpos >= depth(1) & ...
%     sp(probe).templateYpos <= depth(2));
validClu = sp(probe).cids(sp(probe).cgs == 2);

% Plot PSTH for each unit
% relevantStims = length(contrastPars)*length(laserPars) + (-length(laserPars)+1:0);
% trialInds = ismember(stimSeq, relevantStims);
% goodSpikes = ismember(sp(probe).clu, validClu-1);
% psthViewer(sp(probe).st(goodSpikes), ...
%     sp(probe).clu(goodSpikes), stimOn(trialInds), [-.1 stimDur+.1], stimSeq(trialInds))

% Make population plot
window = [-.5 stimDur+.5];
binSize = 0.01;
relevantStims = length(contrastPars)*length(laserPars) + [-length(laserPars)+1 0]; % [no laser, full power laser]
trials_noLaser = ismember(stimSeq, relevantStims(1));
trials_fullLaser = ismember(stimSeq, relevantStims(2));
[psth, bins] = psthAndBA(sp(probe).st(sp(probe).clu==validClu(1)), ...
    stimOn(trials_noLaser), window, binSize);

psths_noLaser = NaN(length(validClu), length(bins));
psths_fullLaser = NaN(length(validClu), length(bins));
depths = NaN(length(validClu), 1);
for n = 1:length(validClu)
    psths_noLaser(n,:) = psthAndBA(sp(probe).st(sp(probe).clu==validClu(n)), ...
        stimOn(trials_noLaser), window, binSize);
    psths_fullLaser(n,:) = psthAndBA(sp(probe).st(sp(probe).clu==validClu(n)), ...
        stimOn(trials_fullLaser), window, binSize);
    depths(n) = mean(sp(probe).spikeDepths(sp(probe).clu==validClu(n)));
end
[depths, order] = sort(depths, 'descend');
psths_noLaser = psths_noLaser(order,:);
psths_fullLaser = psths_fullLaser(order,:);
maxis = max([psths_noLaser, psths_fullLaser], [], 2);
psths_noLaser = 1 - psths_noLaser ./ maxis;
psths_noLaser(all(isnan(psths_noLaser),2),:) = 1;
psths_fullLaser = 1 - psths_fullLaser ./ maxis;
psths_fullLaser(all(isnan(psths_fullLaser),2),:) = 1;

yTicks = unique([1:5:length(depths) length(depths)]);
figure('Position', [500 42 1400 1074])
subplot(1,2,1)
imagesc(bins([1 end]), [1 length(validClu)], psths_noLaser)
colormap gray
hold on
plot([0 0], [.5 length(validClu)+.5], 'r', 'LineWidth', 2)
plot([1 1].*stimDur, [.5 length(validClu)+.5], 'r', 'LineWidth', 2)
set(gca, 'box', 'off', 'YTick', yTicks, 'YTickLabel', round(depths(yTicks)))
xlabel('Time from stim onset (s)')
ylabel('Probe location (microns)')
title('No laser + full contrast grating')

subplot(1,2,2)
imagesc(bins([1 end]), [1 length(validClu)], psths_fullLaser)
colormap gray
hold on
plot([0 0], [.5 length(validClu)+.5], 'r', 'LineWidth', 2)
plot([1 1].*stimDur, [.5 length(validClu)+.5], 'r', 'LineWidth', 2)
set(gca, 'box', 'off', 'YTick', yTicks, 'YTickLabel', round(depths(yTicks)))
xlabel('Time from stim onset (s)')
ylabel('Probe location (microns)')
title('Full laser + full contrast grating')

annotation('textbox', [0 .95 1 .04], 'String', [subject ' ' date], ...
    'FontWeight', 'bold', 'Interpreter', 'none', 'LineStyle', 'none', ...
    'FontSize', 16, 'HorizontalAlignment', 'center')