%% Specify parameters
numDivs = 4;  % number of times the horizontal and vertical axis of each 
              % frame is divided -> determines number of patches analysed
              % separately

%% Initialize info structure of dataset
subject='M160318_SS059';
expDate = '2016-06-29';
exp=1;
planes=3;
channel = '_channel001';
% channel = '';
suffix = 'raw';
% suffix = 'registered';
% suffix = 'rect01_01_registered';

info=ppbox.infoPopulate(subject, expDate, exp);

for iPlane = planes
    %% Extract single planes
    options.nFramesPerChunk = 4096;
    options.iPlane = iPlane;
    ppbox.extractSinglePlane(info, options);
    
    %% Load movie file
    % fileName = sprintf('%s_%s', info.chData(1).basename, suffix);
    fileName = sprintf('%s_plane%03d%s_%s', info.basename2p, iPlane,channel, suffix);
    filePath = fullfile(info.folderProcessed, fileName);
    [data, infoProc] = loadArr(filePath);
    
    %% Load stimulus and time information
    frameTimes = ppbox.getFrameTimes(infoProc, infoProc.planeFrames);
    stimTimes = ppbox.getStimTimes(infoProc);
    [stimFrames, stimFrameTimes, stimPosition] = whiteNoise.getStimulusFrames(infoProc);
    
    %% Get RFs (only spatial) for all brain patches
    RFs = whiteNoise.getRFsOfBrainRegions(data, frameTimes, stimTimes, ...
        stimFrames, stimFrameTimes, numDivs);
    
    %% Plot contours of RFs for all brain patches
    whiteNoise.plotRFsOfBrainRegions(RFs, stimPosition);
end