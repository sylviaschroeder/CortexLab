%% Load meta
subject='M150410_SS045';
expDate='2015-05-04';
exp=1;
planes=1:4;
% subject='M141007_SS029';
% expDate='2014-11-12';
% exp=4;
% planes=2:3;
% subject='M150323_SS042';
% expDate = '2015-04-14';
info=ppbox.infoPopulate(subject, expDate, exp);

%% Get necessary data
% meta structure with calcium traces
plane = 1;
channel = '_channel001';
rect = '_rect01_01';

filePath = fullfile(info.folderProcessed, sprintf('%s_plane%03d%s%s_ROI', ...
    info.basename2p, plane, channel, rect));
[sz, ~, infoROI] = loadArrInfo(strrep(filePath, 'registered', 'ROI'));

% use all data or only neurons
% ind = true(1, size(meta.F,2));
ind = strcmp(infoROI.ROI.CellClasses, 's');

% calcium traces (+ timing) and stimulus information
[stimTimes, stimSequence, stimMatrix, frameTimes, samplingRate] = ...
    ssLocal.getStimulusResponseInfo(infoROI);

% pupil data (independent of plane)
[pupilData, pupilTimes] = ssLocal.loadPupilData(info);
pupilTimes(length(pupilData.x)+1:end) = [];

% running data (independent of plane)
ballData = nonVis.getRunningSpeed(meta);

%% Correct calcium traces
traces = infoROI.F;
traces(infoROI.movingFrames,:) = NaN;

%% Plot orientation tuning
% Plot stimulus-triggered responses of cells
gratings.plotOrientationResponses(traces(:,ind), samplingRate, stimMatrix, stimSequence)

% Plot orientation tuning curves
gratings.plotOrientationTuning(meta.F(:, ind), samplingRate, stimMatrix, stimSequence)

% Plot stimulus-triggered responses to grating with and without distractor outside of RF
gratings.plotOrientationTuningWithCompetitor(meta.F(:,ind), samplingRate, stimMatrix, stimSequence)

% Plot orientation tuning in different conditions
% (1) Running vs. not running
gratings.plotConditionedOrientationTuning(meta, 'running', 1)
% (2) Pupil dilated vs. constricted
gratings.plotConditionedOrientationTuning(meta, 'pupilSize', 1)
% (3) Eye stationary vs. moved from most common position
gratings.plotConditionedOrientationTuning(meta, 'pupilPos', 1)

%% Plot orientation tuning of brain regions
% suffix = 'registered';
% % suffix = 'raw';
% for iPlane = planes
%     filePath = fullfile(info.folderProcessed, ...
%         sprintf('%s_plane%03d_%s', info.basename2p, iPlane, suffix));
%     load(filePath, 'meta')
%     
%     [stimTimes, stimSequence, stimMatrix, frameTimes, samplingRate] = ...
%         ssLocal.getStimulusResponseInfo(meta);
% 
%     filePathMovie = fullfile(info.folderProcessed, [meta.chData(1).basename '_' suffix]);
%     [data, infoReg] = loadArr(filePathMovie);
%     dataPatches = ssLocal.averagePixelsPerPatch(data, 5);
%     gratings.plotTuningOfBrainRegions(dataPatches, samplingRate, stimMatrix, stimSequence, iPlane);
% end

% publish('main_orientation.m','outputDir','\Results_1_tuningOfBrainRegions', 'maxOutputLines', 0, 'figureSnapMethod', 'getframe', 'useNewFigure', false)