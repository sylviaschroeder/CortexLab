%% Parameters
timeGap = 600; %in s, gap between experiments
smoothWheel = 0.15; % in s

%% Folders
folderTools = 'C:\STORAGE\workspaces';
folderScript = 'C:\dev\workspace\CortexLab';

folderBaseUCL = 'C:\STORAGE\OneDrive - University College London\Lab';
folderBaseSussex = 'C:\STORAGE\OneDrive - University of Sussex\Lab';
folderROIData = fullfile(folderBaseUCL, 'DATA\InfoStructs');
folderSave = fullfile(folderBaseSussex, 'DATA\NPY\task_2p');
folderTimeline = '\\ZSERVER.cortexlab.net\Data\expInfo';
folderConfig = '\\ZSERVER.cortexlab.net\Code\Rigging\config';
folderEyeRaw = '\\ZSERVER.cortexlab.net\Data\EyeCamera';
folderEyeDLC = fullfile(folderBaseUCL, 'DATA\DataToPublish\task2P');

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')));
addpath(genpath(fullfile(folderTools, 'wheelAnalysis')));
addpath('\\ZSERVER.cortexlab.net\Code\2photonPipeline')
addpath(genpath('C:\STORAGE\workspaces\Rigbox\cb-tools'))
addpath('C:\STORAGE\workspaces\Rigbox')
addpath('C:\STORAGE\workspaces\Rigbox\cortexlab')
addpath(genpath(fullfile(folderScript)));

%% Load database
db_wheelTask

%% Collect data
for k = 1:length(db)
    fprintf('Dataset %d of %d\n', k, length(db))
    subject = db(k).subject;
    ind = strfind(subject, 'SS');
    subject = subject(ind:end);
    folderSession = fullfile(folderSave, subject, db(k).date);
    if ~isfolder(folderSession)
        mkdir(folderSession)
    end
    
    valid = [];
    xyzPos = [];
    planePerUnit = [];
    cellIDs = [];
    isGad = [];
    infoDone = false;
    F_all = [];
    time = [];
    
    time_behaviour = [];
    wheelPos = [];
    wheelVelocity = [];
    lickPiezo = [];
    
    eyePos = [];
    eyeDiameter = [];
    eyeTime = [];
    
    trialStart = [];
    stimOnTimes = [];
    stimOffTimes = [];
    beepTimes = [];
    feedbackTimes = [];
    trialEnd = [];
    stimulus = [];
    repeated = [];
    choice = [];
    outcome = [];
    rewardAmount = [];
    preStimQuiescentPeriod = [];
    cueInteractiveDelay = [];
    responseWindow = [];
    posFeedbackPeriod = [];
    negFeedbackPeriod = [];
    stimSigma = [];
    stimAzi = [];
    stimAltitude = [];
    stimTargetThreshold = [];
    stimSpatFreq = [];
    stimOri = [];
    
    t0 = 0;
    for iExp = 1:length(db(k).exp)
        folder = fullfile(folderROIData, db(k).subject, ...
            db(k).date, num2str(iExp));
        files = dir(fullfile(folder, '*_ROI.mat'));
        valid_exp = [];
        F_exp = [];
        for iPlane = 1:length(files)
            % load data
            d = load(fullfile(folder, files(iPlane).name));
            meta = d.meta;
            numCells = size(meta.F_final,2);
            
            if ~infoDone
                xyzPos = [xyzPos; cat(1, meta.ROI.CellXYZMicrons{:})];
                planePerUnit = [planePerUnit; ones(numCells,1) .* iPlane];
                cellIDs = [cellIDs; (1:numCells)'];
                isGad = [isGad; meta.ROI.isGad];
            end
            
            if iPlane == 1
                frameTimes = ppbox.getFrameTimes(meta)';
                time = [time; frameTimes + t0];
                
                if iExp == 1
                    timeShifts = zeros(1, length(files));
                end
            elseif iExp == 1
                t = ppbox.getFrameTimes(meta)';
                l = min(length(t), length(time));
                timeShifts(iPlane) = median(t(1:l) - frameTimes(1:l));
            end
            
            F = meta.F_final;
            if size(F_exp,1) > size(F,1)
                F = padarray(F, size(F_exp,1)-size(F,1), NaN, 'post');
            elseif size(F_exp,1)>0 && size(F_exp,1)<size(F,1)
                F(size(F_exp,1)+1:end,:) = [];
            end
            F_exp = [F_exp, F];
            
            % valid ROIs: unique and no "switch-on" responses
            valid_exp = [valid_exp; meta.ROI.isDuplicate == 0 & ...
                meta.ROI.isSwitchOn == 0 & ~all(isnan(meta.F_final),1)'];
        end
        infoDone = true;
        valid = [valid, valid_exp];
        F_all = [F_all; F_exp];
        
        data = load(fullfile(folder, 'timeAlign.mat'));
        blockAlign = data.alignment;
        
        % load timeline (time stamps, rotary encoder)
        folder = fullfile(folderTimeline, db(k).subject, db(k).date, ...
            num2str(db(k).exp(iExp)));
        file = dir(fullfile(folder, '*_Timeline.mat'));
        if length(file) ~= 1
            disp('Timeline file missing or more than 1 exists!');
            return
        end
        data = load(fullfile(folder, file.name));
        timeline = data.Timeline;
        
        switch db(k).microID
            case 'b'
                data = load(fullfile(folderConfig, 'ZMAZE', ...
                    'hardware_20170519.mat'), 'mouseInput');
            case 'b2'
                data = load(fullfile(folderConfig, 'ZURPRISE', ...
                    'hardware20170131.mat'), 'mouseInput');
        end
        mmFactor = data.mouseInput.MillimetresFactor;
        
        t = timeline.rawDAQTimestamps;
        sr = 1 / median(diff(t));
        behTime_exp = round(t(1)*sr)/sr : 1/sr : round(t(end)*sr)/sr;
        wheel_exp = wheel.correctCounterDiscont(timeline.rawDAQData(:, ...
            strcmp({timeline.hw.inputs.name}, 'rotaryEncoder')));
        wheel_exp = interp1(t, wheel_exp, behTime_exp, 'pchip')' .* mmFactor;
        vel_exp = wheel.computeVelocity(wheel_exp', round(smoothWheel * sr), sr)';
        if any(strcmp({timeline.hw.inputs.name}, 'piezoLickDetector'))
            lickPiezo_exp = interp1(t, timeline.rawDAQData(:, ...
                strcmp({timeline.hw.inputs.name}, 'piezoLickDetector')), ...
                behTime_exp)';
        else
            lickPiezo_exp = NaN(length(behTime_exp),1);
        end
        
        % load eye data
        eyePos_exp = readNPY(fullfile(folderEyeDLC, db(k).subject, ...
            db(k).date, num2str(iExp), 'eye.xyPos.npy'));
        eyeDiam_exp = readNPY(fullfile(folderEyeDLC, db(k).subject, ...
            db(k).date, num2str(iExp), 'eye.diameter.npy'));
        f = fullfile(folderEyeRaw, db(k).subject, db(k).date, ...
            num2str(db(k).exp(iExp)), sprintf('%s_%d_%s_eye_TLtime.mat', ...
            db(k).date, db(k).exp(iExp), db(k).subject));
        if isfile(f)
            data = load(f);
            eyeTime_exp = data.pupilTimes;
        else
            fRaw = fullfile(folderEyeRaw, db(k).subject, db(k).date, ...
                num2str(db(k).exp(iExp)), sprintf('%s_%d_%s_eye.mat', ...
                db(k).date, db(k).exp(iExp), db(k).subject));
            if ~isfile(fRaw)
                fprintf('Pupil data does not exist.')
                eyeTime_exp = NaN(length(eyeDiam_exp), 1);
            else
                d = str2double(strrep(db(k).date, '-', ''));
                pupilTimes = eye.getFrameTimes(db(k).subject, d, db(k).exp(iExp));
                save(f, 'pupilTimes');
                eyeTime_exp = pupilTimes;
            end
        end
        
        % block: stim times, conditions, beep times, choice times
        file = dir(fullfile(folder, '*_Block.mat'));
        if length(file) ~= 1
            disp('Block file missing or more than 1 exists!');
            return
        end
        data = load(fullfile(folder, file.name));
        block_exp = data.block;
        block_exp = psy.stripIncompleteTrials(block_exp);
        start_exp = blockAlign.blockToPdTimeFrame([block_exp.trial.trialStartedTime])';
        stimOnTimes_exp = blockAlign.stimOnTimes';
        stimOffTimes_exp = blockAlign.stimOffTimes';
        beepTimes_exp = blockAlign.blockToPdTimeFrame([block_exp.trial.interactiveStartedTime])';
        feedbackTimes_exp = blockAlign.blockToPdTimeFrame([block_exp.trial.feedbackStartedTime])';
        end_exp = blockAlign.blockToPdTimeFrame([block_exp.trial.trialEndedTime])';
        numTrials = length(block_exp.trial);
        stimulus_exp = zeros(numTrials, 2);
        repeated_exp = zeros(numTrials, 1);
        for i = 1:numTrials
            stimulus_exp(i,:) = block_exp.trial(i).condition.visCueContrast;
            repeated_exp(i) = block_exp.trial(i).condition.repeatNum;
        end
        ch = cat(1, block_exp.trial.responseMadeID);
        ch(ch == 1) = -1;
        ch(ch == 2) = 0;
        ch(ch == 3) = 1;
        outcome_exp = cat(1, block_exp.trial.feedbackType) == 1;
        
        time_behaviour = [time_behaviour; behTime_exp'+t0];
        wheelPos = [wheelPos; wheel_exp];
        wheelVelocity = [wheelVelocity; vel_exp];
        lickPiezo = [lickPiezo; lickPiezo_exp];
        eyePos = [eyePos; eyePos_exp];
        eyeDiameter = [eyeDiameter; eyeDiam_exp];
        eyeTime = [eyeTime; eyeTime_exp+t0];
        trialStart = [trialStart; start_exp+t0];
        stimOnTimes = [stimOnTimes; stimOnTimes_exp+t0];
        stimOffTimes = [stimOffTimes; stimOffTimes_exp+t0];
        beepTimes = [beepTimes; beepTimes_exp+t0];
        feedbackTimes = [feedbackTimes; feedbackTimes_exp+t0];
        trialEnd = [trialEnd; end_exp+t0];
        stimulus = [stimulus; stimulus_exp];
        repeated = [repeated; repeated_exp];
        choice = [choice; ch];
        outcome = [outcome; outcome_exp];
        rewardAmount = [rewardAmount; ...
            ones(numTrials,1) * block_exp.parameters.rewardVolume];
        preStimQuiescentPeriod = [preStimQuiescentPeriod; ...
            repmat(block_exp.parameters.preStimQuiescentPeriod', numTrials, 1)];
        cueInteractiveDelay = [cueInteractiveDelay; ...
            repmat(block_exp.parameters.cueInteractiveDelay', numTrials, 1)];
        responseWindow = [responseWindow; ...
            ones(numTrials,1) * block_exp.parameters.responseWindow];
        posFeedbackPeriod = [posFeedbackPeriod; ...
            ones(numTrials,1) * block_exp.parameters.positiveFeedbackPeriod];
        negFeedbackPeriod = [negFeedbackPeriod; ...
            ones(numTrials,1) * block_exp.parameters.negativeFeedbackPeriod];
        stimSigma = [stimSigma; ...
            repmat(block_exp.parameters.cueSigma', numTrials, 1)];
        stimAzi = [stimAzi; ...
            ones(numTrials,1) * block_exp.parameters.distBetweenTargets ./ 2];
        stimAltitude = [stimAltitude; ...
            ones(numTrials,1) * block_exp.parameters.targetAltitude];
        stimTargetThreshold = [stimTargetThreshold; ...
            ones(numTrials,1) * block_exp.parameters.targetThreshold];
        stimSpatFreq = [stimSpatFreq; ...
            ones(numTrials,1) * block_exp.parameters.cueSpatialFrequency];
        stimOri = [stimOri; ...
            ones(numTrials,1) * block_exp.parameters.targetOrientation];
        
        writeNPY(frameTimes([1 end])' + t0, fullfile(folderSession, ...
            sprintf('_ss_recordings.task%02d_intervals.npy', iExp)));
        
        t0 = time(end) + timeGap;
    end
    valid = all(valid,2);
    
    writeNPY(planePerUnit(valid), ...
        fullfile(folderSession, '_ss_2pRois._ss_2pPlanes.npy'));
    writeNPY(cellIDs(valid), ...
        fullfile(folderSession, '_ss_2pRois.ids.npy'));
    writeNPY(xyzPos(valid,:), ...
        fullfile(folderSession, '_ss_2pRois.xyz.npy'));
    isGad(isGad == 0) = NaN;
    writeNPY(isGad(valid), ...
        fullfile(folderSession, '_ss_2pRois.isGad.npy'));
    writeNPY(timeShifts, fullfile(folderSession, '_ss_2pPlanes.delay.npy'));
    writeNPY(time, fullfile(folderSession, '_ss_2pCalcium.timestamps.npy'));
    writeNPY(F_all(:,valid), ...
        fullfile(folderSession, '_ss_2pCalcium.dff.npy'));
    
    writeNPY(time_behaviour, ...
        fullfile(folderSession, '_ss_wheel.timestamps.npy'));
    writeNPY(wheelPos, ...
        fullfile(folderSession, '_ss_wheel.position.npy'));
    writeNPY(wheelVelocity, ...
        fullfile(folderSession, '_ss_wheel.velocity.npy'));
    writeNPY(lickPiezo, ...
        fullfile(folderSession, '_ss_lickPiezo.raw.npy'));
    
    writeNPY(eyePos, fullfile(folderSession, 'eye.xyPos.npy'));
    writeNPY(eyeDiameter, fullfile(folderSession, 'eye.diameter.npy'));
    writeNPY(eyeTime, fullfile(folderSession, 'eye.timestamps.npy'));
    
    writeNPY([trialStart trialEnd], ...
        fullfile(folderSession, '_ss_trials.intervals.npy'));
    writeNPY([stimOnTimes stimOffTimes], ...
        fullfile(folderSession, '_ss_trials.stimOn_intervals.npy'));
    writeNPY(beepTimes, ...
        fullfile(folderSession, '_ss_trials.goCue_times.npy'));
    writeNPY(feedbackTimes, ...
        fullfile(folderSession, '_ss_trials.feedback_times.npy'));
    writeNPY(stimulus(:,1), ...
        fullfile(folderSession, '_ss_trials.contrastLeft.npy'));
    writeNPY(stimulus(:,2), ...
        fullfile(folderSession, '_ss_trials.contrastRight.npy'));
    writeNPY(repeated, ...
        fullfile(folderSession, '_ss_trials.repNum.npy'));
    writeNPY(choice, ...
        fullfile(folderSession, '_ss_trials.choice.npy'));
    writeNPY(outcome, ...
        fullfile(folderSession, '_ss_trials.feedbackType.npy'));
    writeNPY(rewardAmount, ...
        fullfile(folderSession, '_ss_trials.rewardVolume.npy'));
    writeNPY(preStimQuiescentPeriod, ...
        fullfile(folderSession, '_ss_trials.preStimDelay_intervals.npy'));
    writeNPY(cueInteractiveDelay, ...
        fullfile(folderSession, '_ss_trials.interactiveDelay_intervals.npy'));
    writeNPY(responseWindow, ...
        fullfile(folderSession, '_ss_trials.responseWindow.npy'));
    writeNPY(posFeedbackPeriod, ...
        fullfile(folderSession, '_ss_trials.positiveFeedbackPeriod.npy'));
    writeNPY(negFeedbackPeriod, ...
        fullfile(folderSession, '_ss_trials.negativeFeedbackPeriod.npy'));
    writeNPY(stimSigma, ...
        fullfile(folderSession, '_ss_trials.stimSigma.npy'));
    writeNPY(stimAzi, ...
        fullfile(folderSession, '_ss_trials.stimAzimuth.npy'));
    writeNPY(stimAltitude, ...
        fullfile(folderSession, '_ss_trials.stimAltitude.npy'));
    writeNPY(stimTargetThreshold, ...
        fullfile(folderSession, '_ss_trials.stimTargetThreshold.npy'));
    writeNPY(stimSpatFreq, ...
        fullfile(folderSession, '_ss_trials.stimSpatFreq.npy'));
    writeNPY(stimOri, ...
        fullfile(folderSession, '_ss_trials.stimOrientation.npy'));
end
