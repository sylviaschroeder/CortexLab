%% Load meta
subject='M141002_SS026';
expDate='2014-10-29';
% subject='M150323_SS042';
% expDate = '2015-04-14';
exp=4;
info=ppbox.infoPopulate(subject, expDate, exp);

%% Get necessary data
% meta structure with calcium traces
plane = 2;
filePath = fullfile(info.folderProcessed, ...
    sprintf('%s_plane%03d_ROI', info.basename2p, plane));
load(filePath, 'meta')

% calcium traces (+ timing) and stimulus information
[stimTimes, stimSequence, stimMatrix, frameTimes, samplingRate] = ...
    ssLocal.getStimulusResponseInfo(meta);

% pupil data (independent of plane)
[pupilData, pupilTimes] = ssLocal.loadPupilData(info);
pupilTimes(length(pupilData.x)+1:end) = [];

% running data (independent of plane)
ballData = ssLocal.getRunningSpeed(meta);

%% Plot single trial and median responses to all stimuli
ssLocal.plotStimulusResponses(meta.F, samplingRate, stimMatrix, ...
    stimSequence, round(1 * samplingRate), round(4 * samplingRate));

%% Plot position tuning (xy-pos)
% Plot stimulus-triggered responses of cells
ssLocal.mapRetinotopy(traces{plane}, frameTimes, stimMatrix, stimSequence)

neurons = 1:size(deltaTraces,2);
ssLocal.plotPositionTuning(deltaTraces(:,neurons), samplingRate, stimMatrix, stimSequence);

%% Plot orientation tuning
% Plot stimulus-triggered responses of cells
ssLocal.plotOrientationResponses(meta.F, samplingRate, stimMatrix, stimSequence)

% Plot orientation tuning curves
% neurons = [27 25 22 21 20 17 16 15 14 13 12 9 8 6 5 4 2];
% neurons =[45 35 33 29 28 27 26 24 23 22 16 13 12 11 6 5 3 2];
neurons = [47 40 35 32 28 27 26 24 23 22 13 9 6 5 3 1];
ssLocal.plotOrientationTuning(meta.F(:, neurons), samplingRate, stimMatrix, stimSequence)

% Plot stimulus-triggered responses to grating with and without distractor outside of RF
ssLocal.plotOrientationTuningWithCompetitor(meta.F, samplingRate, stimMatrix, stimSequence)

% Plot orientation tuning in different conditions
% (1) Running vs. not running
gratings.plotConditionedOrientationTuning(meta, 'running', 1)
% (2) Pupil dilated vs. constricted
gratings.plotConditionedOrientationTuning(meta, 'pupilSize', 1)
% (3) Eye stationary vs. moved from most common position
gratings.plotConditionedOrientationTuning(meta, 'pupilPos', 1)

%% Plot pupil data, running data, and Ca-traces
nonVis.plotPupilRunningAndCaTraces(meta.F, frameTimes, pupilData, pupilTimes, ...
    ballData);