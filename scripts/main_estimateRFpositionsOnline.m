%% Specify parameters
numDivs = 4;  % number of times the horizontal and vertical axis of each 
              % frame is divided -> determines number of patches analysed
              % separately

%% Initialize info structure of dataset
subject='M150114_SS035';
expDate = '2015-02-10';
exp=1;
iPlane = 3;

info=ppbox.infoPopulate(subject, expDate, exp);

%% Extract single planes
% Go to local folder first!!! (to save temporary files)
options.nFramesPerChunk = 4096;
options.iPlane = iPlane;
ppbox.extractSinglePlane(info, options);

%% Load unregistered file
fileName = sprintf('%s_plane%03d_raw', info.basename2p, iPlane);
filePath = fullfile(info.folderProcessed, fileName);
[data, infoRaw] = loadArr(filePath);

%% Load stimulus and time information
frameTimes = ppbox.getFrameTimes(infoRaw);
stimTimes = ssLocal.getStimTimes(infoRaw);
[stimFrames, stimFrameTimes, stimPosition] = whiteNoise.getStimulusFrames(infoRaw);

%% Get RFs (only spatial) for all brain patches
RFs = whiteNoise.getRFsOfBrainRegions(data, frameTimes, stimTimes, ...
    stimFrames, stimFrameTimes, numDivs);

%% Plot contours of RFs for all brain patches
whiteNoise.plotRFsOfBrainRegions(RFs, stimPosition);