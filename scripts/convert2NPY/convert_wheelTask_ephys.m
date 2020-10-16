% Need to run
% C:\dev\workspace\CortexLab\scripts\electrophys\main_preprocess.m first!

%% Folders
% folderTools = 'C:\Users\Flora\Github';
folderTools = 'C:\STORAGE\workspaces';
% folderData = 'Z:';
folderData = '\\zubjects.cortexlab.net\Subjects';
folderEyeData = '\\ZSERVER.cortexlab.net\Data\EyeCamera';
% folderScript = 'C:\Users\Flora\Github\CortexLab\scripts\passiveAnalysis';
folderScript = 'C:\dev\workspace\CortexLab';

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')));
addpath(genpath(fullfile(folderTools, 'spikes')));
addpath(genpath(fullfile(folderTools, 'kilotrodeRig')));
addpath(genpath(fullfile(folderTools, 'alyx-matlab')));
addpath(genpath(fullfile(folderScript)));

%% Define dataset and prepare folders
db = db_ephys_task;
k = 26;

mouseName = db(k).subject; 
date = db(k).date;

root = fullfile(folderData, mouseName, date);

% make alf folder within root
alfDir = fullfile(root, 'alf');
if ~exist(alfDir, 'dir')
   mkdir(alfDir)
end

% make align folder
alignDir = fullfile(root, 'alignments');

%% Basic info
tlToMaster = readNPY(fullfile(alignDir, ...
    sprintf('correct_timeline_%d_to_ephys_%s.npy', db(k).expTL, db(k).probes{1})));
% load timeline
d = load(fullfile(root, num2str(db(k).expTL), ...
    sprintf('%s_%d_%s_Timeline.mat', date, db(k).expTL, mouseName)));
tl = d.Timeline;

%% write ephys data for each probe
ephys.writeEphysToNPY(mouseName, date, db(k).probes);

% write position of probes to csv file
hemisphere = cell(length(db(k).probes),1);
location = cell(length(db(k).probes),1);
probe = string(db(k).probes');
for p = 1:length(probe)
    if contains(db(k).probeLocations{p}, 'L')
        hemisphere{p} = 'left';
    else
        hemisphere{p} = 'right';
    end
    if contains(db(k).probeLocations{p}, 'A')
        location{p} = 'anterior';
    else
        location{p} = 'posterior';
    end
end
hemisphere = string(hemisphere);
location = string(location);
tbl = table(probe, hemisphere, location);
writetable(tbl, fullfile(alfDir, '_ss_probes.location.csv'), ...
    'Delimiter', '\t');

%% write timeline data
smoothWheel = 0.15;
dt = median(diff(tl.rawDAQTimestamps));
% write timeline time
writeNPY(applyCorrection(tl.rawDAQTimestamps, tlToMaster), ...
    fullfile(alfDir, '_ss_signals.timestamps.npy'))
% write wheel signal
wheelPos = wheel.correctCounterDiscont(tl.rawDAQData(:, ...
    strcmp({tl.hw.inputs.name}, 'rotaryEncoder')))';
writeNPY(wheelPos, fullfile(alfDir, '_ss_wheel.position.npy'));
vel = wheel.computeVelocity(wheelPos, round(smoothWheel / dt), 1/dt)';
writeNPY(vel, fullfile(alfDir, '_ss_wheel.velocity.npy'));
% write pieze lick signal
lickPiezo = tl.rawDAQData(:, strcmp({tl.hw.inputs.name}, 'piezoLickDetector'));
writeNPY(lickPiezo, fullfile(alfDir, '_ss_lickPiezo.raw.npy'));
clear lickPiezo vel wheelPos

%% write time of eye video frames
% load eye video frame times
filename1 = fullfile(root, num2str(db(k).expTL), 'eye_timeStamps.mat');
filename2 = fullfile(folderEyeData, mouseName, date, num2str(db(k).expTL), ...
    'eye_timeStamps.mat');

if ~isfile(filename1) && ~isfile(filename2)
    % run \kilotrodeRig\alignVideo! then run this cell again
%     addpath('C:\STORAGE\workspaces\Rigbox')
%     addpath('C:\STORAGE\workspaces\Rigbox\cortexlab')
%     addpath(genpath('C:\STORAGE\workspaces\Rigbox\cb-tools'))
    alignVideo(mouseName, date, db(k).expTL, 'eye');
end

if isfile(filename1)
    d = load(filename1);
    writeNPY(applyCorrection(d.tVid, tlToMaster), ...
        fullfile(alfDir, 'eye.timestamps.npy'));
elseif isfile(filename2)
    d = load(filename2);
    writeNPY(applyCorrection(d.tVid, tlToMaster), ...
        fullfile(alfDir, 'eye.timestamps.npy'));
end
clear d

%% write darkness data
% plot photodiode signal
figure('Position', [34 558 1850 420])
photodiode = tl.rawDAQData(:, strcmp({tl.hw.inputs.name}, 'photoDiode'));
t = applyCorrection(tl.rawDAQTimestamps, tlToMaster);
plot(t, photodiode)
marked = false;
while ~marked
    title('Mark the start of the dark period by placing the crosshair!')
    roiStart = drawcrosshair(gca);
    title('Now mark the end of the dark period by placing the crosshair!')
    roiEnd = drawcrosshair(gca, 'Color', 'r');
    dark = [roiStart.Position(1) roiEnd.Position(1)];
    roiStart.delete;
    roiEnd.delete;
    inds = t > dark(1) & t < dark(2);
    rng = diff(dark);
    m = median(photodiode(inds));
    s = std(photodiode(inds));
    h = drawrectangle('Position', [dark(1) m-7*s rng 14*s]);
    title('Happy? Press enter or y!')
    xlim(dark + [-.1*rng .1*rng])
    ylim(m + [-10 10] .* s)
    str = input('Happy with the marked period? y/n [y]: ', 's');
    if strcmp(str,'y') || isempty(str)
        marked = true;
    else
        h.delete;
    end
end
close gcf
writeNPY(dark, fullfile(alfDir, '_ss_recordings.darkness_intervals.npy'))
clear photodiode t 

%% write choice world data
% Note: for all experiments: interactive delay was NOT extended when animal
% moved wheel during the delay

trialTimes_appr = [];
stimTimes_appr = [];
goCueTimes_appr = [];
feedbackTimes_appr = [];
contrasts = [];
feedbackType = [];
choice = [];
repNum = [];
preStim = [];
delays = [];
reward = [];
delayIntvals = [];
negFb = [];
posFb = [];
preStimIntvals = [];
win = [];
stimAlt = [];
stimAzi = [];
stimOri = [];
stimSig = [];
stimSpatFr = [];
stimThresh = [];
interTrial = [];
goCueFreq = [];
stimOnTimes = [];
stimOffTimes = [];
goCueTimes = [];
extraValveClicks = [];
feedbackTimes = [];

for bl = 1:length(db(k).expTask)
    fprintf('\nTask block %d\n\n', bl)
    
    % load block data
    d = load(fullfile(root, num2str(db(k).expTask(bl)), ...
        sprintf('%s_%d_%s_Block.mat', date, db(k).expTask(bl), mouseName)));
    block = d.block;
    
    if isfield(block, 'trial') % choice world
        % load linear fit from block times to timeline times
        blockToTL = readNPY(fullfile(alignDir, ...
            sprintf('correct_block_%d_to_timeline_%d.npy', db(k).expTask(bl), db(k).expTL)));
        % data from block
        trials = block.trial;
        % excluded unfinished trials
        numTrials = length([trials.trialEndedTime]);
        trials = trials(1:numTrials);
        % trial parameters
        cond = [trials.condition];
        
        % times from block
        trT = reshape(applyCorrection( ...
            [[trials.trialStartedTime]'; [trials.trialEndedTime]'], blockToTL), ...
            [], 2);
        trialTimes_appr = [trialTimes_appr; reshape(applyCorrection(trT, ...
            tlToMaster), [], 2)];
        stT = reshape(applyCorrection( ...
            [[trials.stimulusCueStartedTime]'; [trials.stimulusCueEndedTime]'], blockToTL), ...
            [], 2);
        stimTimes_appr = [stimTimes_appr; reshape(applyCorrection(stT, ...
            tlToMaster), [], 2)];
        goT = applyCorrection([trials.interactiveStartedTime]', blockToTL);
        goCueTimes_appr = [goCueTimes_appr; applyCorrection(goT, tlToMaster)];
        fbT = applyCorrection([trials.feedbackStartedTime]', blockToTL);
        feedbackTimes_appr = [feedbackTimes_appr; applyCorrection(fbT, tlToMaster)];
        % trial parameters/events from block
        contrasts = [contrasts; [cond.visCueContrast]']; % [left right]
        fbType = [trials.feedbackType]' == 1;
        feedbackType = [feedbackType; fbType]; % logical
        ch = [trials.responseMadeID]'; % -1 = left, 0 = still, 1 = right
        ch(ch == 1) = -1;
        ch(ch == 2) = 0;
        ch(ch == 3) = 1;
        choice = [choice; ch];
        repNum = [repNum; [cond.repeatNum]'];
        preStim = [preStim; NaN(numTrials,1)];
        delays = [delays; NaN(numTrials,1)];
        % block parameters
        if isfield(block.parameters, 'rewardVolume')
            reward = [reward; repmat(block.parameters.rewardVolume, numTrials, 1)];
        else
            reward = [reward; [cond.rewardVolume]'];
        end
        delayIntvals = [delayIntvals; repmat(block.parameters.cueInteractiveDelay, numTrials, 1)];
        negFb = [negFb; repmat(block.parameters.negativeFeedbackPeriod, numTrials, 1)];
        posFb = [posFb; repmat(block.parameters.positiveFeedbackPeriod, numTrials, 1)];
        preStimIntvals = [preStimIntvals; repmat(block.parameters.preStimQuiescentPeriod', numTrials, 1)];
        win = [win; repmat(block.parameters.responseWindow, numTrials, 1)];
        stimAlt = [stimAlt; repmat(block.parameters.targetAltitude, numTrials, 1)];
        stimAzi = [stimAzi; repmat(block.parameters.distBetweenTargets, numTrials, 1) ./ 2];
        stimOri = [stimOri; repmat(block.parameters.targetOrientation, numTrials, 1)];
        stimSig = [stimSig; repmat(block.parameters.cueSigma', numTrials, 1)];
        stimSpatFr = [stimSpatFr; repmat(block.parameters.cueSpatialFrequency, numTrials, 1)];
        stimThresh = [stimThresh; repmat(block.parameters.targetThreshold, numTrials, 1)];
        interTrial = [interTrial; repmat(block.parameters.interTrialDelay, numTrials, 1)];
        goCueFreq = [goCueFreq; repmat(block.parameters.onsetToneFreq, numTrials, 1)];
    else % signals
        file = fullfile(alignDir, ...
            sprintf('correct_block_%d_to_timeline_wheel_%d.npy', ...
            db(k).expTask(bl), db(k).expTL));
        if isfile(file)
            blockToTL = readNPY(file);
        else
            blockToTL = preproc.alignBlock2TL(block, tl, true);
            writeNPY(blockToTL, file);
        end
        % data from block
        trials = block.events;
        % excluded unfinished trials
        numTrials = length(trials.endTrialTimes);
        
        % times from block (for SS093 2018-05-24 1st trial, events in order of
        % appearance, time in brackets relative to experiment start time)
        % repeatNumTimes (0.0061) -> contrastLeft/contrastRight (0.0165) ->
        % preStimulusDelay (0.0166) -> interactiveDelay (0.0168) ->
        % contrast (0.0231) -> newTrial (0.034) -> trialNum (0.035) ->
        % stimulusOn (0.404) -> interactiveOn (1.21) -> response (1.718) ->
        % feedback (1.824) -> stimulusOff (2.723) -> endTrial (3.626)
        trT = reshape(applyCorrection( ...
            [trials.repeatNumTimes(1:numTrials)'; ...
            trials.endTrialTimes(1:numTrials)'], blockToTL), [], 2);
        trialTimes_appr = [trialTimes_appr; reshape(applyCorrection(trT, tlToMaster), [], 2)];
        if isfield(trials, 'stimulusOffTimes')
            stT = reshape(applyCorrection( ...
                [trials.stimulusOnTimes(1:numTrials)'; ...
                trials.stimulusOffTimes(1:numTrials)'], blockToTL), [], 2);
        else
            stT = [applyCorrection(trials.stimulusOnTimes(1:numTrials), blockToTL), ...
                NaN(numTrials,1)];
        end
        stimTimes_appr = [stimTimes_appr; reshape(applyCorrection(stT, tlToMaster), [], 2)];
        goT = applyCorrection(trials.interactiveOnTimes(1:numTrials)', blockToTL);
        goCueTimes_appr = [goCueTimes_appr; applyCorrection(goT, tlToMaster)];
        fbT = applyCorrection(trials.feedbackTimes(1:numTrials)', blockToTL);
        feedbackTimes_appr = [feedbackTimes_appr; applyCorrection(fbT, tlToMaster)];
        % trial events from block
        contrasts = [contrasts; [trials.contrastLeftValues(1:numTrials)' ...
            trials.contrastRightValues(1:numTrials)']];
        fbType = trials.feedbackValues(1:numTrials)';
        feedbackType = [feedbackType; fbType];
        choice = [choice; trials.responseValues(1:numTrials)'];
        repNum = [repNum; trials.repeatNumValues(1:numTrials)'];
        valid = trials.preStimulusDelayTimes' > trials.expStartTimes;
        ps = trials.preStimulusDelayValues(valid)';
        preStim = [preStim; ps(1:numTrials)];
        valid = trials.interactiveDelayTimes' > trials.expStartTimes;
        dl = trials.interactiveDelayValues(valid)';
        delays = [delays; dl(1:numTrials)];
        % block parameters
        reward = [reward; [block.paramsValues(1:numTrials).rewardSize]'];
        delayIntvals = [delayIntvals; [block.paramsValues(1:numTrials).interactiveDelay]'];
        negFb = [negFb; [block.paramsValues(1:numTrials).noiseBurstDur]'];
        posFb = [posFb; NaN(numTrials,1)];
        preStimIntvals = [preStimIntvals; [block.paramsValues(1:numTrials).preStimulusDelay]'];
        win = [win; [block.paramsValues(1:numTrials).responseWindow]'];
        stimAlt = [stimAlt; [block.paramsValues(1:numTrials).stimulusAltitude]'];
        stimAzi = [stimAzi; [block.paramsValues(1:numTrials).stimulusAzimuth]'];
        stimOri = [stimOri; [block.paramsValues(1:numTrials).stimulusOrientation]'];
        stimSig = [stimSig; [block.paramsValues(1:numTrials).sigma]'];
        stimSpatFr = [stimSpatFr; [block.paramsValues(1:numTrials).spatialFrequency]'];
        stimThresh = [stimThresh; [block.paramsValues(1:numTrials).stimulusAzimuth]'];
        interTrial = [interTrial; [block.paramsValues(1:numTrials).interTrialDelay]'];
        goCueFreq = [goCueFreq; [block.paramsValues(1:numTrials).onsetToneFrequency]'];
    end
    
    % stimulus times
    tlStimUpdateTimes = readNPY(fullfile(alignDir, ...
        sprintf('block_%d_sw_in_timeline_%d.npy', db(k).expTask(bl), db(k).expTL)));
    % find those photodiode events that are closest to stim onset times
    [indUpdates, indStims, timeDiffs] = ...
        preproc.findMatchingTimes(tlStimUpdateTimes, stT(:,1));
    valid = abs(timeDiffs) < 0.1;
    onT = NaN(numTrials,1);
    onT(indStims(valid)) = tlStimUpdateTimes(indUpdates(valid));
    stimOnTimes = [stimOnTimes; onT];
    % find those photodiode events that are closest to stim offset times
    [indUpdates, indStims, timeDiffs] = ...
        preproc.findMatchingTimes(tlStimUpdateTimes, stT(:,2));
    valid = abs(timeDiffs) < 0.1;
    offT = NaN(numTrials,1);
    offT(indStims(valid)) = tlStimUpdateTimes(indUpdates(valid));
    stimOffTimes = [stimOffTimes; offT];
    if any(any(isnan([onT offT]), 2), 1)
        fprintf('%d stimulus onsets and %d stimulus offsets were NOT detected!\n', ...
            sum(isnan(onT)), sum(isnan(offT)))
    end
    ax = [0 0];
    figure('Position', [30 140 1880 840])
    subplot(2,1,1)
    hold on
    plot(tl.rawDAQTimestamps, tl.rawDAQData(:, strcmp({tl.hw.inputs.name}, ...
        'photoDiode')), 'k')
    plot(onT, ones(length(onT),1), 'r>')
    plot(offT, ones(length(offT),1), 'b<')
    ylabel('Photodiode')
    title(sprintf('Stimulus on and off times (%s, %s, exp %d)', mouseName, date, db(k).expTask(bl)))
    ax(1) = gca;
    subplot(2,1,2)
    hold on
    plot(tlStimUpdateTimes, zeros(length(tlStimUpdateTimes),1), 'ok')
    plot(onT, zeros(length(onT),1), 'rx')
    plot(offT, zeros(length(offT),1), 'bx')
    xlabel('Time (s, tl time)')
    ylabel('Update times of photodiode')
    ax(2) = gca;
    linkaxes(ax, 'x')
    if isnan(offT(end))
        xlim([onT(find(~isnan(onT),1))-3 onT(end)+10])
    else
        xlim([onT(find(~isnan(onT),1))-3 offT(end)+3])
    end
    
    % times from timeline
    tlTaskInds = tl.rawDAQTimestamps > trT(1,1) & ...
        tl.rawDAQTimestamps < trT(end,2);
    t = tl.rawDAQTimestamps(tlTaskInds);
    tlBinSize = median(diff(t));
    
    % go cue/beep times
    audio = tl.rawDAQData(tlTaskInds, strcmp({tl.hw.inputs.name}, 'audioMonitor'));
    smAudio = conv(audio.^2, gausswin(7)./sum(gausswin(7)), 'same');
    [~, beeps, beepsOff] = schmittTimes(t, smAudio, [0.02 0.04]);
    [indBeeps, indGoCues, timeDiffs] = preproc.findMatchingTimes(beeps, goT);
    valid = abs(timeDiffs) < 0.15;
    gcT = NaN(numTrials,1);
    gcT(indGoCues(valid)) = beeps(indBeeps(valid));
    goCueTimes = [goCueTimes; gcT];
    if any(isnan(gcT))
        fprintf('%d go cues were NOT detected!\n', sum(isnan(gcT)))
    end
    figure('Position', [30 560 1880 420])
    hold on
    plot(tl.rawDAQTimestamps(tlTaskInds), smAudio, 'k')
    plot(gcT, ones(length(gcT),1).*0.02, 'mx')
    xlabel('Time (s, tl time)')
    ylabel('Smoothed audio')
    title(sprintf('Go cue times (%s, %s, exp %d)', mouseName, date, db(k).expTask(bl)))
    axis tight
    
    % feedback times
    fdbkT = NaN(numTrials,1);
    posTrials = find(fbType);
    negTrials = find(~fbType);
    % reward times
    rewards = tl.rawDAQData(tlTaskInds, strcmp({tl.hw.inputs.name}, 'rewardEcho'));
    [~, valveClicks] = schmittTimes(t, rewards, [2 3]);
    [indClicks, indFeedback, timeDiffs] = ...
        preproc.findMatchingTimes(valveClicks, fbT(posTrials,1));
    valid = abs(timeDiffs) < 0.12;
    fdbkT(posTrials(indFeedback(valid))) = valveClicks(indClicks(valid));
    extraValveClicks = [extraValveClicks; valveClicks(setdiff(1:length(valveClicks), indClicks(valid)))];
    if length(indFeedback(valid)) < length(posTrials)
        fprintf('%d reward times were NOT detected!\n', length(posTrials)-sum(valid))
    end
    figure('Position', [30 560 1880 420])
    hold on
    plot(tl.rawDAQTimestamps(tlTaskInds), rewards, 'k')
    plot(valveClicks, ones(length(valveClicks),1), 'co')
    xlabel('Time (s, tl time)')
    ylabel('Reward echo')
    title(sprintf('Reward times (%s, %s, exp %d)', mouseName, date, db(k).expTask(bl)))
    axis tight
    % white noise sound times
    % first remove samples corresponding to beeps
    startSamps = round((beeps - t(1)) / tlBinSize) + 1;
    endSamps = round((beepsOff - t(1)) / tlBinSize) + 1;
    for bb = 1:length(startSamps)
%         audio(startSamps(bb)-8:endSamps(bb)+8) = 0;
        audio(startSamps(bb)-30:endSamps(bb)+30) = 0;
    end
    smAudio = conv(audio.^2, gausswin(20)./sum(gausswin(20)), 'same');
    [~, whiteNoise] = schmittTimes(t, smAudio, [1e-4 1e-3]);
    [indNoise, indFeedback, timeDiffs] = ...
        preproc.findMatchingTimes(whiteNoise, fbT(negTrials));
    valid = abs(timeDiffs) < 0.2;
    if sum(valid) < length(negTrials)/2 % if detected white noise burst are very few, noise burst were probably not recorded and detected times are by chance
        valid(:) = false;
    end
    fdbkT(negTrials(indFeedback(valid))) = whiteNoise(indNoise(valid));
    feedbackTimes = [feedbackTimes; fdbkT];
    if length(indFeedback(valid)) < length(negTrials)
        fprintf('%d noise times were NOT detected!\n', length(negTrials)-sum(valid))
    end
    figure('Position', [30 560 1880 420])
    hold on
    plot(tl.rawDAQTimestamps(tlTaskInds), smAudio, 'k')
    plot(whiteNoise(indNoise(valid)), ones(sum(valid),1).*.001, 'go')
    xlabel('Time (s, tl time)')
    ylabel('Smoothed audio w/o beeps')
    title(sprintf('Whitenoise times (%s, %s, exp %d)', mouseName, date, db(k).expTask(bl)))
    axis tight
end

% write variables
writeNPY(trialTimes_appr, fullfile(alfDir, '_ss_trials.intervals.npy'));
writeNPY(stimTimes_appr, fullfile(alfDir, '_ss_trials.approxStimOn_intervals.npy'));
writeNPY(goCueTimes_appr, fullfile(alfDir, '_ss_trials.approxGoCue_times.npy'));
writeNPY(feedbackTimes_appr, fullfile(alfDir, '_ss_trials.approxFeedback_times.npy'));

writeNPY(contrasts(1,:), fullfile(alfDir, '_ss_trials.contrastLeft.npy'));
writeNPY(contrasts(2,:), fullfile(alfDir, '_ss_trials.contrastRight.npy'));
writeNPY(feedbackType, fullfile(alfDir, '_ss_trials.feedbackType.npy'));
writeNPY(choice, fullfile(alfDir, '_ss_trials.choice.npy'));
writeNPY(repNum, fullfile(alfDir, '_ss_trials.repNum.npy'));
writeNPY(preStim, fullfile(alfDir, '_ss_trials.preStimDelays.npy'));
writeNPY(delays, fullfile(alfDir, '_ss_trials.interactiveDelays.npy'));

writeNPY(reward, fullfile(alfDir, '_ss_trials.rewardVolume.npy'));
writeNPY(delayIntvals, fullfile(alfDir, '_ss_trials.interactiveDelay_intervals.npy'));
writeNPY(negFb, fullfile(alfDir, '_ss_trials.negativeFeedbackPeriod.npy'));
writeNPY(posFb, fullfile(alfDir, '_ss_trials.positiveFeedbackPeriod.npy'));
writeNPY(preStimIntvals, fullfile(alfDir, '_ss_trials.preStimDelay_intervals.npy'));
writeNPY(win, fullfile(alfDir, '_ss_trials.responseWindow.npy'));
writeNPY(stimAlt, fullfile(alfDir, '_ss_trials.stimAltitude.npy'));
writeNPY(stimAzi, fullfile(alfDir, '_ss_trials.stimAzimuth.npy'));
writeNPY(stimOri, fullfile(alfDir, '_ss_trials.stimOrientation.npy'));
writeNPY(stimSig, fullfile(alfDir, '_ss_trials.stimSigma.npy'));
writeNPY(stimSpatFr, fullfile(alfDir, '_ss_trials.stimSpatFreq.npy'));
writeNPY(stimThresh, fullfile(alfDir, '_ss_trials.stimTargetThreshold.npy'));
writeNPY(interTrial, fullfile(alfDir, '_ss_trials.interTrialDelay.npy'));
writeNPY(goCueFreq, fullfile(alfDir, '_ss_trials.goCueFreq.npy'));

writeNPY(reshape(applyCorrection([stimOnTimes; stimOffTimes], tlToMaster), [], 2), ...
    fullfile(alfDir, '_ss_trials.stimOn_intervals.npy'));
writeNPY(applyCorrection(goCueTimes, tlToMaster), ...
    fullfile(alfDir, '_ss_trials.goCue_times.npy'));
writeNPY(applyCorrection(extraValveClicks, tlToMaster), ...
    fullfile(alfDir, '_ss_extraReward.times.npy'));
writeNPY(applyCorrection(feedbackTimes, tlToMaster), ...
    fullfile(alfDir, '_ss_trials.feedback_times.npy'));

%% write sparse noise data
% (1) normal (white on black), e.g. SS087 2017-12-12 4
%   10 rows, 36 columns
%   azimuth: [-135 135]
%   elevation: [-37.5 37.5]
%   square size: 1
%   presentation time: 1/6 s
% (2) white and black squares on gray, e.g. SS093 2018-05-24 5
% 
% (3) only left visual field, white on black, e.g. SS093 2018-05-27 13
%
% (4) only right visual field, white on black, e.g. SS093 2018-05-27 14
%
% (5) left monitor covered, normal, e.g. SS093 2018-05-27 15
%
% (6) right monitor covered, normal, e.g. SS093 2018-05-27 16

expNames = {'expNoise', 'expNoiseGray', 'expNoiseLeftOnly', ...
    'expNoiseRightOnly', 'expNoiseLeftCovered', 'expNoiseRightCovered'};
fileNames = {'', 'Gray', 'LeftOnly', 'RightOnly', 'LeftCovered', 'RightCovered'};
positions = [-135 135 -37.5 37.5; -135 135 -37.5 37.5; ...
              135   0 -37.5 37.5;    0 135 -37.5 37.5; ...
             -135 135 -37.5 37.5; -135 135 -37.5 37.5];

for exp = 1:length(expNames)
    if isempty(db(k).(expNames{exp}))
        continue
    end
    % load block file
    d = load(fullfile(root, num2str(db(k).(expNames{exp})), ...
        sprintf('%s_%d_%s_Block.mat', date, db(k).(expNames{exp}), mouseName)));
    block = d.block;
    
    % get times of stimulus updates
    stimTimes = readNPY(fullfile(alignDir, ...
        sprintf('block_%d_sw_in_timeline_%d.npy', ...
        db(k).(expNames{exp}), db(k).expTL)));
    % convert times to master times
    stimTimes = applyCorrection(stimTimes, tlToMaster);
    % get stimulus frames
    stimFrames = block.events.stimuliOnValues;
    numRows = size(stimFrames,1);
    stimFrames = reshape(stimFrames, [], length(stimTimes)); % [px of one frame x time]
    % determine unique frames and indexing into those unique frames
    [uniFrames,~,inds] = unique(stimFrames', 'rows');
    % reshape frames to [rows x columns x time]
    uniFrames = reshape(uniFrames', numRows, [], size(uniFrames,1));            
    
    writeNPY(stimTimes, fullfile(alfDir, ...
        sprintf('_ss_sparseNoise%s.times.npy', fileNames{exp})));
    writeNPY(inds, fullfile(alfDir, ...
        sprintf('_ss_sparseNoise%s._ss_sparseNoise%sID.npy', fileNames{exp}, fileNames{exp})));
    writeNPY(positions(exp,:), fullfile(alfDir, ...
        sprintf('_ss_sparseNoise%sArea.edges.npy', fileNames{exp})));
    writeNPY(uniFrames, fullfile(alfDir, ...
        sprintf('_ss_sparseNoise%sID.map.npy', fileNames{exp})));
end
clear block d fr stimFrames uniFrames

%% write passive viewing data
% Note: in signals, noise bursts don't seem to be recorded -> not saved in
% feedback_times

block = [];

if ~isempty(db(k).expPassive)
    % load block data
    d = load(fullfile(root, num2str(db(k).expPassive), ...
        sprintf('%s_%d_%s_Block.mat', date, db(k).expPassive, mouseName)));
    block = d.block;
end

if isfield(block, 'trial') % choice world
    % load linear fit from block times to timeline times
    blockToTL = readNPY(fullfile(alignDir, ...
        sprintf('correct_block_%d_to_timeline_%d.npy', db(k).expPassive, db(k).expTL)));
    % data from block
    trials = block.trial;
    % excluded unfinished trials
    numTrials = length([trials.trialEndedTime]);
    trials = trials(1:numTrials);
    % trial parameters
    cond = [trials.condition];
    
    % times from block
    trialTimes_appr = reshape(applyCorrection(applyCorrection( ...
        [[trials.trialStartedTime]'; [trials.trialEndedTime]'], blockToTL), ...
        tlToMaster), [], 2);
    stimTimes_appr = reshape(applyCorrection(applyCorrection( ...
        [[trials.stimulusCueStartedTime]'; [trials.stimulusCueEndedTime]'], blockToTL), ...
        tlToMaster), [], 2);
    goCueTimes_appr = applyCorrection(applyCorrection( ...
        [trials.interactiveStartedTime]', blockToTL), tlToMaster);
    feedbackTimes_appr = applyCorrection(applyCorrection( ...
        [trials.feedbackStartedTime]', blockToTL), tlToMaster);
    % trial parameters/events from block
    contrasts = [cond.visCueContrast]'; % [left right]
    fbType = NaN(numTrials,1); % 1: valve click, 0: noise burst, NaN: neither
    fbType([trials.feedbackType]==1) = 1;
    fbType([cond.negFeedbackSoundAmp]>0) = 0;
    feedbackTimes_appr(isnan(fbType)) = NaN;
    stimAlt = [cond.targetAltitude]';
    stimAzi = [cond.distBetweenTargets]' ./ 2;
    % block parameters
    reward = repmat(block.parameters.rewardVolume, numTrials, 1);
    reward([trials.feedbackType] == -1) = 0;
    negFb = repmat(block.parameters.negativeFeedbackPeriod, numTrials, 1);
    posFb = repmat(block.parameters.positiveFeedbackPeriod, numTrials, 1);
    win = repmat(block.parameters.responseWindow, numTrials, 1);
    stimOri = repmat(block.parameters.targetOrientation, numTrials, 1);
    stimSig = repmat(block.parameters.cueSigma', numTrials, 1);
    stimSpatFr = repmat(block.parameters.cueSpatialFrequency, numTrials, 1);
    stimThresh = repmat(block.parameters.targetThreshold, numTrials, 1);
    interTrial = repmat(block.parameters.interTrialDelay, numTrials, 1); % here: 2 columns denoting a period
    goCue = [cond.interactiveOnsetToneRelAmp]' > 0.01;
    goCueTimes_appr(~goCue) = NaN;
    goCueFreq = repmat(block.parameters.onsetToneFreq, numTrials, 1);
elseif ~isempty(block) % signals
    file = fullfile(alignDir, ...
        sprintf('correct_block_%d_to_timeline_wheel_%d.npy', ...
        db(k).expPassive, db(k).expTL));
    if isfile(file)
        blockToTL = readNPY(file);
    else
        blockToTL = preproc.alignBlock2TL(block, tl, true);
        str = input('Happy with alignment? y/n [y]: ', 's');
        if strcmp(str,'y') || isempty(str)
            writeNPY(blockToTL, file);
        end
    end
    % data from block
    trials = block.events;
    % excluded unfinished trials
    numTrials = length(trials.endTrialTimes);
    
    % times from block (for SS093 2018-05-24 1st trial, events in order of 
    % appearance, time in brackets relative to experiment start time)
    % repeatNumTimes (0.0061) -> contrastLeft/contrastRight (0.0165) -> 
    % preStimulusDelay (0.0166) -> interactiveDelay (0.0168) ->
    % contrast (0.0231) -> newTrial (0.034) -> trialNum (0.035) ->
    % stimulusOn (0.404) -> interactiveOn (1.21) -> response (1.718) -> 
    % feedback (1.824) -> stimulusOff (2.723) -> endTrial (3.626)
    trialTimes_appr = reshape(applyCorrection(applyCorrection( ...
        [trials.repeatNumTimes(1:numTrials)'; ...
        trials.endTrialTimes(1:numTrials)'], blockToTL), tlToMaster), [], 2);
    stimTimes_appr = reshape(applyCorrection(applyCorrection( ...
        [trials.stimulusOnTimes(1:numTrials)'; ...
        trials.stimulusOffTimes(1:numTrials)'], blockToTL), tlToMaster), [], 2);
    goCueTimes_appr = applyCorrection(applyCorrection( ...
        trials.interactiveOnTimes(1:numTrials)', blockToTL), tlToMaster);
    feedbackTimes_appr = applyCorrection(applyCorrection( ...
        trials.feedbackTimes(1:numTrials)', blockToTL), tlToMaster);
    % trial events from block
    contrasts = [trials.contrastLeftValues(1:numTrials)' ...
        trials.contrastRightValues(1:numTrials)'];
    % block parameters
    reward = [block.paramsValues(1:numTrials).rewardSize]';
    negFb = [block.paramsValues(1:numTrials).noiseBurstDur]';
    posFb = NaN(numTrials,1);
    win = [block.paramsValues(1:numTrials).responseWindow]';
    stimAlt = [block.paramsValues(1:numTrials).stimulusAltitude]';
    stimAzi = [block.paramsValues(1:numTrials).stimulusAzimuth]';
    stimOri = [block.paramsValues(1:numTrials).stimulusOrientation]';
    stimSig = [block.paramsValues(1:numTrials).sigma]';
    stimSpatFr = [block.paramsValues(1:numTrials).spatialFrequency]';
    interTrial = [block.paramsValues(1:numTrials).interTrialDelay]'; % here: 1 column denoting exact values
    goCue = [block.paramsValues(1:numTrials).onsetToneAmplitude]' > 0.01;
    goCueTimes_appr(~goCue) = NaN;
    goCueFreq = [block.paramsValues(1:numTrials).onsetToneFrequency]';
    fbType = NaN(numTrials,1); % 1: valve click, 0: noise burst, NaN: neither
    fbType(reward > 0) = 1;
    noiseAmp = [block.paramsValues(1:numTrials).noiseBurstAmp]';
    fbType(noiseAmp > 0) = 0;
    feedbackTimes_appr(isnan(fbType)) = NaN;
end

if ~isempty(block)
    % write variables
    writeNPY(trialTimes_appr, fullfile(alfDir, '_ss_passive.intervals.npy'));
    writeNPY(stimTimes_appr, fullfile(alfDir, '_ss_passive.approxStimOn_intervals.npy'));
    writeNPY(goCueTimes_appr, fullfile(alfDir, '_ss_passive.approxGoCue_times.npy'));
    writeNPY(feedbackTimes_appr, fullfile(alfDir, '_ss_passive.approxFeedback_times.npy'));
    
    writeNPY(contrasts(1,:), fullfile(alfDir, '_ss_passive.contrastLeft.npy'));
    writeNPY(contrasts(2,:), fullfile(alfDir, '_ss_passive.contrastRight.npy'));
    writeNPY(fbType, fullfile(alfDir, '_ss_passive.feedbackType.npy'));
    
    writeNPY(reward, fullfile(alfDir, '_ss_passive.rewardVolume.npy'));
    writeNPY(negFb, fullfile(alfDir, '_ss_passive.negativeFeedbackPeriod.npy'));
    writeNPY(posFb, fullfile(alfDir, '_ss_passive.positiveFeedbackPeriod.npy'));
    writeNPY(win, fullfile(alfDir, '_ss_passive.responseWindow.npy'));
    writeNPY(stimAlt, fullfile(alfDir, '_ss_passive.stimAltitude.npy'));
    writeNPY(stimAzi, fullfile(alfDir, '_ss_passive.stimAzimuth.npy'));
    writeNPY(stimOri, fullfile(alfDir, '_ss_passive.stimOrientation.npy'));
    writeNPY(stimSig, fullfile(alfDir, '_ss_passive.stimSigma.npy'));
    writeNPY(stimSpatFr, fullfile(alfDir, '_ss_passive.stimSpatFreq.npy'));
    writeNPY(interTrial, fullfile(alfDir, '_ss_passive.interTrialDelay.npy'));
    writeNPY(goCueFreq, fullfile(alfDir, '_ss_passive.goCueFreq.npy'));
    writeNPY(goCue, fullfile(alfDir, '_ss_passive.goCuePlayed.npy'));
    
    % stimulus times
    tlStimUpdateTimes = readNPY(fullfile(alignDir, ...
        sprintf('block_%d_sw_in_timeline_%d.npy', db(k).expPassive, db(k).expTL)));
    % find those photodiode events that are closest to stim onset times
    [indUpdates, indStims, timeDiffs] = preproc.findMatchingTimes(tlStimUpdateTimes, ...
        (stimTimes_appr(:,1) - tlToMaster(2)) ./ tlToMaster(1));
    valid = abs(timeDiffs) < 0.1;
    stimOnTimes = NaN(numTrials,1);
    stimOnTimes(indStims(valid)) = tlStimUpdateTimes(indUpdates(valid));
    % find those photodiode events that are closest to stim offset times
    [indUpdates, indStims, timeDiffs] = preproc.findMatchingTimes(tlStimUpdateTimes, ...
        (stimTimes_appr(:,2) - tlToMaster(2)) ./ tlToMaster(1));
    valid = abs(timeDiffs) < 0.1;
    stimOffTimes = NaN(numTrials,1);
    stimOffTimes(indStims(valid)) = tlStimUpdateTimes(indUpdates(valid));
    writeNPY(reshape(applyCorrection([stimOnTimes; stimOffTimes], tlToMaster), [], 2), ...
        fullfile(alfDir, '_ss_passive.stimOn_intervals.npy'));
    if any(any(isnan([stimOnTimes stimOffTimes]), 2), 1)
        fprintf('%d stimulus onsets and %d stimulus offsets were NOT detected!\n', ...
            sum(isnan(stimOnTimes)), sum(isnan(stimOffTimes)))
    end
    ax = [0 0];
    figure('Position', [30 140 1880 840])
    subplot(2,1,1)
    hold on
    plot(tl.rawDAQTimestamps, tl.rawDAQData(:, strcmp({tl.hw.inputs.name}, ...
        'photoDiode')), 'k')
    plot(stimOnTimes, ones(length(stimOnTimes),1), 'r>')
    plot(stimOffTimes, ones(length(stimOffTimes),1), 'b<')
    ylabel('Photodiode')
    title(sprintf('Stimulus on and off times (%s, %s)', mouseName, date))
    ax(1) = gca;
    subplot(2,1,2)
    hold on
    plot(tlStimUpdateTimes, zeros(length(tlStimUpdateTimes),1), 'ok')
    plot(stimOnTimes, zeros(length(stimOnTimes),1), 'rx')
    plot(stimOffTimes, zeros(length(stimOffTimes),1), 'bx')
    xlabel('Time (s, tl time)')
    ylabel('Update times of photodiode')
    ax(2) = gca;
    linkaxes(ax, 'x')
    xlim([stimOnTimes(find(~isnan(stimOnTimes),1))-3 stimOffTimes(end)+3])
    
    % times from timeline
    tlTaskInds = tl.rawDAQTimestamps > (trialTimes_appr(1,1) - tlToMaster(2)) ./ tlToMaster(1) & ...
        tl.rawDAQTimestamps < (trialTimes_appr(end,2) - tlToMaster(2)) ./ tlToMaster(1);
    t = tl.rawDAQTimestamps(tlTaskInds);
    tlBinSize = median(diff(t));
    
    % go cue/beep times
    audio = tl.rawDAQData(tlTaskInds, strcmp({tl.hw.inputs.name}, 'audioMonitor'));
    smAudio = conv(audio.^2, gausswin(7)./sum(gausswin(7)), 'same');
    [~, beeps, beepsOff] = schmittTimes(t, smAudio, [0.02 0.04]);
    [indBeeps, indGoCues, timeDiffs] = preproc.findMatchingTimes(beeps, ...
        (goCueTimes_appr(:,1) - tlToMaster(2)) ./ tlToMaster(1));
    valid = abs(timeDiffs) < 0.15;
    goCueTimes = NaN(numTrials,1);
    goCueTimes(indGoCues(valid)) = beeps(indBeeps(valid));
    writeNPY(applyCorrection(goCueTimes, tlToMaster), ...
        fullfile(alfDir, '_ss_passive.goCue_times.npy'));
    if sum(goCueTimes) ~= sum(goCue) && any(isnan(goCueTimes(goCue)))
        fprintf('%d go cues were NOT detected!\n', sum(~valid))
    end
    figure('Position', [30 560 1880 420])
    hold on
    plot(tl.rawDAQTimestamps(tlTaskInds), smAudio, 'k')
    plot(goCueTimes, ones(length(goCueTimes),1).*0.02, 'mx')
    xlabel('Time (s, tl time)')
    ylabel('Smoothed audio')
    title(sprintf('Go cue times (%s, %s)', mouseName, date))
    axis tight
    
    % feedback times
    feedbackTimes = NaN(numTrials,1);
    posTrials = find(fbType == 1);
    negTrials = find(fbType == 0);
    % reward times
    rewards = tl.rawDAQData(tlTaskInds, strcmp({tl.hw.inputs.name}, 'rewardEcho'));
    [~, valveClicks] = schmittTimes(t, rewards, [2 3]);
    [indClicks, indFeedback, timeDiffs] = preproc.findMatchingTimes(valveClicks, ...
        (feedbackTimes_appr(posTrials,1) - tlToMaster(2)) ./ tlToMaster(1));
    valid = abs(timeDiffs) < 0.1;
    feedbackTimes(posTrials(indFeedback(valid))) = valveClicks(indClicks(valid));
    if length(indFeedback(valid)) < length(posTrials)
        fprintf('%d reward times were NOT detected!\n', length(posTrials)-sum(valid))
    end
    figure('Position', [30 560 1880 420])
    hold on
    plot(tl.rawDAQTimestamps(tlTaskInds), rewards, 'k')
    plot(valveClicks, ones(length(valveClicks),1), 'co')
    xlabel('Time (s, tl time)')
    ylabel('Reward echo')
    title(sprintf('Reward times (%s, %s)', mouseName, date))
    axis tight
    % white noise sound times
    % first remove samples corresponding to beeps
    startSamps = round((beeps - t(1)) / tlBinSize) + 1;
    endSamps = round((beepsOff - t(1)) / tlBinSize) + 1;
    for bb = 1:length(startSamps)
        audio(startSamps(bb)-8:endSamps(bb)+8) = 0;
    end
    smAudio = conv(audio.^2, gausswin(20)./sum(gausswin(20)), 'same');
    [~, whiteNoise] = schmittTimes(t, smAudio, [1e-4 1e-3]);
    [indNoise, indFeedback, timeDiffs] = preproc.findMatchingTimes(whiteNoise, ...
        (feedbackTimes_appr(negTrials,1) - tlToMaster(2)) ./ tlToMaster(1));
    valid = abs(timeDiffs) < 0.2;
    feedbackTimes(negTrials(indFeedback(valid))) = whiteNoise(indNoise(valid));
    writeNPY(applyCorrection(feedbackTimes, tlToMaster), ...
        fullfile(alfDir, '_ss_passive.feedback_times.npy'));
    if length(indFeedback(valid)) < length(negTrials)
        fprintf('%d noise times were NOT detected!\n', length(negTrials)-sum(valid))
    end
    figure('Position', [30 560 1880 420])
    hold on
    plot(tl.rawDAQTimestamps(tlTaskInds), smAudio, 'k')
    plot(whiteNoise(indNoise(valid)), ones(sum(valid),1).*.001, 'go')
    xlabel('Time (s, tl time)')
    ylabel('Smoothed audio w/o beeps')
    title(sprintf('Whitenoise times (%s, %s)', mouseName, date))
    axis tight
end