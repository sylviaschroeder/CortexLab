%% Parameters
timeGap = 600; %in s, gap between experiments
smoothWheel = 0.15; % in s

%% Folders
% folderTools = 'C:\STORAGE\workspaces';
% folderScript = 'C:\dev\workspace\CortexLab';
% 
% folderBaseUCL = 'C:\STORAGE\OneDrive - University College London\Lab';
% folderBaseSussex = 'C:\STORAGE\OneDrive - University of Sussex\Lab';
% folderROIData = fullfile(folderBaseUCL, 'DATA\InfoStructs');
% folderSave = fullfile(folderBaseSussex, 'DATA\NPY\task_2p');
% folderTimeline = '\\ZSERVER.cortexlab.net\Data\expInfo';
% folderConfig = '\\ZSERVER.cortexlab.net\Code\Rigging\config';
% folderEyeRaw = '\\ZSERVER.cortexlab.net\Data\EyeCamera';
% folderEyeDLC = fullfile(folderBaseUCL, 'DATA\DataToPublish\task2P');

folderTools = 'C:\dev\toolboxes';
folderScript = 'C:\dev\workspace\CortexLab';

folderBaseSussex = 'C:\Users\Sylvia\OneDrive - University of Sussex\Lab';
folderROIData = fullfile('Z:\UCLData\2P_Task');
folderSave = fullfile(folderBaseSussex, 'DATA\NPY\task_2p');
folderTimeline = 'Z:\UCLData\2P_Task';
folderEyeRaw = 'Z:\UCLData\2P_Task';
folderEyeDLC = 'Z:\UCLData\2P_Task';
folderXFiles = fullfile(folderBaseSussex, 'DATA\Raw\UCL data\xfiles');
folderStimulus = fullfile(folderBaseSussex, 'DATA\Raw\UCL data\Stimulus');

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')));
addpath(genpath(fullfile(folderTools, 'wheelAnalysis')));
addpath(fullfile(folderBaseSussex, 'CODE_CortexLab'));
addpath(folderXFiles);
addpath(folderStimulus);
% addpath(genpath('C:\STORAGE\workspaces\Rigbox\cb-tools'))
% addpath('C:\STORAGE\workspaces\Rigbox')
% addpath('C:\STORAGE\workspaces\Rigbox\cortexlab')
% addpath('\\ZSERVER.cortexlab.net\Code\Spikes')
% addpath('\\ZSERVER.cortexlab.net\Data\xfiles')
addpath(genpath(fullfile(folderScript)));

%% Load database
db_wheelTask

%% Collect data
for k = 4:length(db)
    if isempty(db(k).expNoise)
        continue
    end
    
    fprintf('Dataset %d of %d\n', k, length(db))
    subject = db(k).subject;
    ind = strfind(subject, 'SS');
    subject = subject(ind:end);
    folderSession = fullfile(folderSave, subject, db(k).date);
    if ~isfolder(folderSession)
        mkdir(folderSession)
    end
    
%     valid = [];
%     xyzPos = [];
%     planePerUnit = [];
%     cellIDs = [];
%     isGad = [];
%     infoDone = false;
%     F_all = [];
%     time = [];
    
%     time_behaviour = [];
%     wheelPos = [];
%     wheelVelocity = [];
%     moveOn = [];
%     moveOff = [];
%     moveDispl = [];
%     peakVelTime = [];
%     peakAmp = [];
%     lickPiezo = [];
%     
%     eyePos = [];
%     eyeDiameter = [];
%     eyeTime = [];
%     
%     trialStart = [];
%     stimOnTimes = [];
%     stimOffTimes = [];
%     beepTimes = [];
%     feedbackTimes = [];
%     trialEnd = [];
%     stimulus = [];
%     repeated = [];
%     choice = [];
%     outcome = [];
%     rewardAmount = [];
%     preStimQuiescentPeriod = [];
%     cueInteractiveDelay = [];
%     responseWindow = [];
%     posFeedbackPeriod = [];
%     negFeedbackPeriod = [];
%     stimSigma = [];
%     stimAzi = [];
%     stimAltitude = [];
%     stimTargetThreshold = [];
%     stimSpatFreq = [];
%     stimOri = [];
    
%     t0 = 0;
    % task data
%     for iExp = 1:length(db(k).expTask)
%         folder = fullfile(folderROIData, db(k).subject, ...
%             db(k).date, num2str(iExp));
%         files = dir(fullfile(folder, '*_ROI.mat'));
%         valid_exp = [];
%         F_exp = [];
%         for iPlane = 1:length(files)
%             % load data
%             d = load(fullfile(folder, files(iPlane).name));
%             meta = d.meta;
%             numCells = size(meta.F_final,2);
%             
%             if ~infoDone
%                 xyzPos = [xyzPos; cat(1, meta.ROI.CellXYZMicrons{:})];
%                 planePerUnit = [planePerUnit; ones(numCells,1) .* iPlane];
%                 cellIDs = [cellIDs; (1:numCells)'];
%                 isGad = [isGad; meta.ROI.isGad];
%             end
%             
%             if iPlane == 1
%                 frameTimes = ppbox.getFrameTimes(meta)';
%                 time = [time; frameTimes + t0];
%                 
%                 if iExp == 1
%                     timeShifts = zeros(1, length(files));
%                 end
%             elseif iExp == 1
%                 t = ppbox.getFrameTimes(meta)';
%                 l = min(length(t), length(time));
%                 timeShifts(iPlane) = median(t(1:l) - frameTimes(1:l));
%             end
%             
%             F = meta.F_final;
%             if size(F_exp,1) > size(F,1)
%                 F = padarray(F, size(F_exp,1)-size(F,1), NaN, 'post');
%             elseif size(F_exp,1)>0 && size(F_exp,1)<size(F,1)
%                 F(size(F_exp,1)+1:end,:) = [];
%             end
%             F_exp = [F_exp, F];
%             
%             % valid ROIs: unique and no "switch-on" responses
%             valid_exp = [valid_exp; meta.ROI.isDuplicate == 0 & ...
%                 meta.ROI.isSwitchOn == 0 & ~all(isnan(meta.F_final),1)'];
%         end
%         infoDone = true;
%         valid = [valid, valid_exp];
%         F_all = [F_all; F_exp];
        
%         data = load(fullfile(folder, 'timeAlign.mat'));
%         blockAlign = data.alignment;
%         
%         % load timeline (time stamps, rotary encoder)
%         folder = fullfile(folderTimeline, db(k).subject, db(k).date, ...
%             num2str(db(k).expTask(iExp)));
%         file = dir(fullfile(folder, '*_Timeline.mat'));
%         if length(file) ~= 1
%             disp('Timeline file missing or more than 1 exists!');
%             return
%         end
%         data = load(fullfile(folder, file.name));
%         timeline = data.Timeline;
        
%         switch db(k).microID
%             case 'b'
%                 data = load(fullfile(folderConfig, 'ZMAZE', ...
%                     'hardware_20170519.mat'), 'mouseInput');
%             case 'b2'
%                 data = load(fullfile(folderConfig, 'ZURPRISE', ...
%                     'hardware20170131.mat'), 'mouseInput');
%         end
%         mmFactor = data.mouseInput.MillimetresFactor;

%         % NOTE: mmFactor hardcoded on 30.05.2023
%         switch db(k).microID
%             case 'b'
%                 mmFactor = 0.7775441817;
%             case 'b2'
%                 mmFactor = 0.38878905747;
%         end
%         
%         t = timeline.rawDAQTimestamps;
%         sr = 1 / median(diff(t));
%         behTime_exp = round(t(1)*sr)/sr : 1/sr : round(t(end)*sr)/sr;
%         wheel_exp = wheel.correctCounterDiscont(timeline.rawDAQData(:, ...
%             strcmp({timeline.hw.inputs.name}, 'rotaryEncoder')));
%         wheel_exp = interp1(t, wheel_exp, behTime_exp, 'pchip')' .* mmFactor;
%         vel_exp = wheel.computeVelocity(wheel_exp', round(smoothWheel * sr), sr)';
%         % write wheel movement events
%         Fs = round(sr);
%         [moveOn_exp, moveOff_exp, moveDispl_exp, peakVelTime_exp, peakAmp_exp] = ...
%             wheel.findWheelMoves3(wheel_exp, t, Fs, 'posThresh', 0.02, 'tThresh', 0.2, ...
%             'posThreshOnset', 0.001);
%         
%         if any(strcmp({timeline.hw.inputs.name}, 'piezoLickDetector'))
%             lickPiezo_exp = interp1(t, timeline.rawDAQData(:, ...
%                 strcmp({timeline.hw.inputs.name}, 'piezoLickDetector')), ...
%                 behTime_exp)';
%         else
%             lickPiezo_exp = NaN(length(behTime_exp),1);
%         end
%         
%         % load eye data
%         eyePos_exp = readNPY(fullfile(folderEyeDLC, db(k).subject, ...
%             db(k).date, num2str(iExp), 'eye.xyPos.npy'));
%         eyeDiam_exp = readNPY(fullfile(folderEyeDLC, db(k).subject, ...
%             db(k).date, num2str(iExp), 'eye.diameter.npy'));
%         file = dir(fullfile(folderEyeRaw, db(k).subject, db(k).date, ...
%             num2str(db(k).expTask(iExp)), '*_eye_TLtime.mat'));
%         f = fullfile(folderEyeRaw, db(k).subject, db(k).date, ...
%             num2str(db(k).expTask(iExp)), file.name);
%         if isfile(f)
%             data = load(f);
%             eyeTime_exp = data.pupilTimes;
%         else
%             fprintf('Pupil data does not exist.')
%             eyeTime_exp = NaN(length(eyeDiam_exp), 1);
%         end
%         
%         % block: stim times, conditions, beep times, choice times
%         file = dir(fullfile(folder, '*_Block.mat'));
%         if length(file) ~= 1
%             disp('Block file missing or more than 1 exists!');
%             return
%         end
%         warning off
%         data = load(fullfile(folder, file.name));
%         warning on
%         block_exp = data.block;
%         block_exp = beh.stripIncompleteTrials(block_exp);
%         start_exp = blockAlign.blockToPdTimeFrame([block_exp.trial.trialStartedTime])';
%         stimOnTimes_exp = blockAlign.stimOnTimes';
%         stimOffTimes_exp = blockAlign.stimOffTimes';
%         beepTimes_exp = blockAlign.blockToPdTimeFrame([block_exp.trial.interactiveStartedTime])';
%         feedbackTimes_exp = blockAlign.blockToPdTimeFrame([block_exp.trial.feedbackStartedTime])';
%         end_exp = blockAlign.blockToPdTimeFrame([block_exp.trial.trialEndedTime])';
%         numTrials = length(block_exp.trial);
%         stimulus_exp = zeros(numTrials, 2);
%         repeated_exp = zeros(numTrials, 1);
%         for i = 1:numTrials
%             stimulus_exp(i,:) = block_exp.trial(i).condition.visCueContrast;
%             repeated_exp(i) = block_exp.trial(i).condition.repeatNum;
%         end
%         ch = cat(1, block_exp.trial.responseMadeID);
%         ch(ch == 1) = -1;
%         ch(ch == 2) = 1;
%         ch(ch == 3) = 0;
%         outcome_exp = cat(1, block_exp.trial.feedbackType) == 1;
%         
%         time_behaviour = [time_behaviour; behTime_exp'+t0];
%         wheelPos = [wheelPos; wheel_exp];
%         wheelVelocity = [wheelVelocity; vel_exp];
%         moveOn = [moveOn; moveOn_exp];
%         moveOff = [moveOff; moveOff_exp];
%         moveDispl = [moveDispl; moveDispl_exp];
%         peakVelTime = [peakVelTime; peakVelTime_exp];
%         peakAmp = [peakAmp; peakAmp_exp];
%         lickPiezo = [lickPiezo; lickPiezo_exp];
%         eyePos = [eyePos; eyePos_exp];
%         eyeDiameter = [eyeDiameter; eyeDiam_exp];
%         eyeTime = [eyeTime; eyeTime_exp+t0];
%         trialStart = [trialStart; start_exp+t0];
%         stimOnTimes = [stimOnTimes; stimOnTimes_exp+t0];
%         stimOffTimes = [stimOffTimes; stimOffTimes_exp+t0];
%         beepTimes = [beepTimes; beepTimes_exp+t0];
%         feedbackTimes = [feedbackTimes; feedbackTimes_exp+t0];
%         trialEnd = [trialEnd; end_exp+t0];
%         stimulus = [stimulus; stimulus_exp];
%         repeated = [repeated; repeated_exp];
%         choice = [choice; ch];
%         outcome = [outcome; outcome_exp];
%         rewardAmount = [rewardAmount; ...
%             ones(numTrials,1) * block_exp.parameters.rewardVolume];
%         preStimQuiescentPeriod = [preStimQuiescentPeriod; ...
%             repmat(block_exp.parameters.preStimQuiescentPeriod', numTrials, 1)];
%         cueInteractiveDelay = [cueInteractiveDelay; ...
%             repmat(block_exp.parameters.cueInteractiveDelay', numTrials, 1)];
%         responseWindow = [responseWindow; ...
%             ones(numTrials,1) * block_exp.parameters.responseWindow];
%         posFeedbackPeriod = [posFeedbackPeriod; ...
%             ones(numTrials,1) * block_exp.parameters.positiveFeedbackPeriod];
%         negFeedbackPeriod = [negFeedbackPeriod; ...
%             ones(numTrials,1) * block_exp.parameters.negativeFeedbackPeriod];
%         stimSigma = [stimSigma; ...
%             repmat(block_exp.parameters.cueSigma', numTrials, 1)];
%         stimAzi = [stimAzi; ...
%             ones(numTrials,1) * block_exp.parameters.distBetweenTargets ./ 2];
%         stimAltitude = [stimAltitude; ...
%             ones(numTrials,1) * block_exp.parameters.targetAltitude];
%         stimTargetThreshold = [stimTargetThreshold; ...
%             ones(numTrials,1) * block_exp.parameters.targetThreshold];
%         stimSpatFreq = [stimSpatFreq; ...
%             ones(numTrials,1) * block_exp.parameters.cueSpatialFrequency];
%         stimOri = [stimOri; ...
%             ones(numTrials,1) * block_exp.parameters.targetOrientation];
%         
%         % only when 2P frametimes already saved
%         frameTimes = ppbox.getFrameTimes2023(timeline);
%         time = [time; frameTimes + t0];
%         %--------------------------------------
%         writeNPY(frameTimes([1 end])' + t0, fullfile(folderSession, ...
%             sprintf('_ss_recordings.task%02d_intervals.npy', iExp)));
%         
%         t0 = time(end) + timeGap;
%     end

    % for cases where 2P data were converted using prerocess_2023, find t0
    % for visual noise experiments by loading _ss_2pCalcium.timestamps
    % (this assumes that visual noise is always the 2nd experiment)
    frameTimes = readNPY(fullfile(folderSession, '_ss_2pCalcium.timestamps.npy'));
    ind = find(diff(frameTimes) > 300);
    t0 = frameTimes(ind) + timeGap;
    frameTimes(1:ind) = [];
    
    % visual noise data
    if ~isempty(db(k).expNoise)
%         folder = fullfile(folderROIData, db(k).subject, ...
%             db(k).date, num2str(db(k).expNoise));
%         files = dir(fullfile(folder, '*_ROI.mat'));
%         valid_exp = [];
%         F_exp = [];
        for iPlane = 1 %:length(files)
%             % load data
%             d = load(fullfile(folder, files(iPlane).name));
%             meta = d.meta;
%             % find correct expInfo folder
%             stimFolder = dir(fullfile(folderTimeline, ['M*_' subject]));
%             meta.subject = stimFolder.name;
%             meta.folderTL = fullfile(stimFolder.folder, stimFolder.name, ...
%                 db(k).date, num2str(meta.exp));
%             meta.basenameTL = sprintf('%s_%d_%s_Timeline', db(k).date, ...
%                 meta.exp, stimFolder.name);
%             numCells = size(meta.F_final,2);

            meta.folderTL = fullfile(folderTimeline, db(k).subject, db(k).date, ...
                num2str(db(k).expNoise));
            file = dir(fullfile(meta.folderTL, '*_Timeline.mat'));
            meta.basenameTL = file.name;
            
            if iPlane == 1
%                 data = load(fullfile(meta.folderTL, file.name));
%                 timeline = data.Timeline;
%                 frameTimes = ppbox.getFrameTimes2023(timeline);
%                 time = [time; frameTimes + t0];
                
                stimTimes = ppbox.getStimTimes(meta);
%                 [stimFrames, stimPosition] = whiteNoise.getStimulusFrames(meta);
                fileHardware = dir(fullfile(folderTimeline, db(k).subject, db(k).date, ...
                    num2str(db(k).expNoise), '*_hardwareInfo.mat'));
                [stimFrames, stimPosition] = ...
                    whiteNoise.getStimulusFrames2023(...
                    fullfile(folderTimeline, db(k).subject, db(k).date, ...
                    num2str(db(k).expNoise), 'Protocol.mat'), ...
                    fullfile(folderTimeline, db(k).subject, db(k).date, ...
                    num2str(db(k).expNoise), fileHardware.name));
                stimFrames = permute(stimFrames, [3 1 2]);
                stimFrameDur = mean(stimTimes.offset - stimTimes.onset) / size(stimFrames,1);
                stimFrameTimes = ((1:size(stimFrames,1))-1) .* stimFrameDur;
                allStimFrameTimes = reshape((stimTimes.onset + stimFrameTimes)', [], 1);
                
                writeNPY(allStimFrameTimes + t0, ...
                    fullfile(folderSession, '_ss_sparseNoise.times.npy'));
                writeNPY(repmat((1:size(stimFrames,1))', length(stimTimes.onset), 1), ...
                    fullfile(folderSession, '_ss_sparseNoise._ss_sparseNoiseID.npy'));
                writeNPY(stimPosition, ...
                    fullfile(folderSession, '_ss_sparseNoiseArea.edges.npy'));
                writeNPY(stimFrames, ...
                    fullfile(folderSession, '_ss_sparseNoiseID.map.npy'));
            end
            
%             F = meta.F_final;
%             if size(F_exp,1) > size(F,1)
%                 F = padarray(F, size(F_exp,1)-size(F,1), NaN, 'post');
%             elseif size(F_exp,1)>0 && size(F_exp,1)<size(F,1)
%                 F(size(F_exp,1)+1:end,:) = [];
%             end
%             F_exp = [F_exp, F];
%             
%             % valid ROIs: unique and no "switch-on" responses
%             valid_exp = [valid_exp; meta.ROI.isDuplicate == 0 & ...
%                 meta.ROI.isSwitchOn == 0 & ~all(isnan(meta.F_final),1)'];
        end
%         valid = [valid, valid_exp];
%         F_all = [F_all; F_exp];
        
        writeNPY(frameTimes([1 end])', fullfile(folderSession, ...
            '_ss_recordings.sparseNoise_intervals.npy'));
    end
%     valid = all(valid,2);
%     
%     writeNPY(planePerUnit(valid), ...
%         fullfile(folderSession, '_ss_2pRois._ss_2pPlanes.npy'));
%     writeNPY(cellIDs(valid), ...
%         fullfile(folderSession, '_ss_2pRois.ids.npy'));
%     writeNPY(xyzPos(valid,:), ...
%         fullfile(folderSession, '_ss_2pRois.xyz.npy'));
%     isGad(isGad == 0) = NaN;
%     writeNPY(isGad(valid), ...
%         fullfile(folderSession, '_ss_2pRois.isGad.npy'));
%     writeNPY(timeShifts, fullfile(folderSession, '_ss_2pPlanes.delay.npy'));
%     writeNPY(time, fullfile(folderSession, '_ss_2pCalcium.timestamps.npy'));
%     writeNPY(F_all(:,valid), ...
%         fullfile(folderSession, '_ss_2pCalcium.dff.npy'));
    
%     writeNPY(time_behaviour, ...
%         fullfile(folderSession, '_ss_wheel.timestamps.npy'));
%     writeNPY(wheelPos, ...
%         fullfile(folderSession, '_ss_wheel.position.npy'));
%     writeNPY(wheelVelocity, ...
%         fullfile(folderSession, '_ss_wheel.velocity.npy'));
%     writeNPY([moveOn moveOff], ...
%         fullfile(folderSession, '_ss_wheelMoves.intervals.npy'));
%     writeNPY(moveDispl, ...
%         fullfile(folderSession, '_ss_wheelMoves.displacement.npy'));
%     writeNPY(peakVelTime, ...
%         fullfile(folderSession, '_ss_wheelMoves.peak_times.npy'));
%     writeNPY(peakAmp, ...
%         fullfile(folderSession, '_ss_wheelMoves.amplitude.npy'));
%     writeNPY(lickPiezo, ...
%         fullfile(folderSession, '_ss_lickPiezo.raw.npy'));
%     
%     writeNPY(eyePos, fullfile(folderSession, 'eye.xyPos.npy'));
%     writeNPY(eyeDiameter, fullfile(folderSession, 'eye.diameter.npy'));
%     writeNPY(eyeTime, fullfile(folderSession, 'eye.timestamps.npy'));
%     
%     writeNPY([trialStart trialEnd], ...
%         fullfile(folderSession, '_ss_trials.intervals.npy'));
%     writeNPY([stimOnTimes stimOffTimes], ...
%         fullfile(folderSession, '_ss_trials.stimOn_intervals.npy'));
%     writeNPY(beepTimes, ...
%         fullfile(folderSession, '_ss_trials.goCue_times.npy'));
%     writeNPY(feedbackTimes, ...
%         fullfile(folderSession, '_ss_trials.feedback_times.npy'));
%     writeNPY(stimulus(:,1), ...
%         fullfile(folderSession, '_ss_trials.contrastLeft.npy'));
%     writeNPY(stimulus(:,2), ...
%         fullfile(folderSession, '_ss_trials.contrastRight.npy'));
%     writeNPY(repeated, ...
%         fullfile(folderSession, '_ss_trials.repNum.npy'));
%     writeNPY(choice, ...
%         fullfile(folderSession, '_ss_trials.choice.npy'));
%     writeNPY(outcome, ...
%         fullfile(folderSession, '_ss_trials.feedbackType.npy'));
%     writeNPY(rewardAmount, ...
%         fullfile(folderSession, '_ss_trials.rewardVolume.npy'));
%     writeNPY(preStimQuiescentPeriod, ...
%         fullfile(folderSession, '_ss_trials.preStimDelay_intervals.npy'));
%     writeNPY(cueInteractiveDelay, ...
%         fullfile(folderSession, '_ss_trials.interactiveDelay_intervals.npy'));
%     writeNPY(responseWindow, ...
%         fullfile(folderSession, '_ss_trials.responseWindow.npy'));
%     writeNPY(posFeedbackPeriod, ...
%         fullfile(folderSession, '_ss_trials.positiveFeedbackPeriod.npy'));
%     writeNPY(negFeedbackPeriod, ...
%         fullfile(folderSession, '_ss_trials.negativeFeedbackPeriod.npy'));
%     writeNPY(stimSigma, ...
%         fullfile(folderSession, '_ss_trials.stimSigma.npy'));
%     writeNPY(stimAzi, ...
%         fullfile(folderSession, '_ss_trials.stimAzimuth.npy'));
%     writeNPY(stimAltitude, ...
%         fullfile(folderSession, '_ss_trials.stimAltitude.npy'));
%     writeNPY(stimTargetThreshold, ...
%         fullfile(folderSession, '_ss_trials.stimTargetThreshold.npy'));
%     writeNPY(stimSpatFreq, ...
%         fullfile(folderSession, '_ss_trials.stimSpatFreq.npy'));
%     writeNPY(stimOri, ...
%         fullfile(folderSession, '_ss_trials.stimOrientation.npy'));
end
