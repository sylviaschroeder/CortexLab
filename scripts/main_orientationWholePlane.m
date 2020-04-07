%% Specify parameters
numDivs = 4;  % number of times the horizontal and vertical axis of each 
              % frame is divided -> determines number of patches analysed
              % separately
              
%% Initialize info structure of dataset
subject='M141002_SS026';
expDate = '2014-10-29';
exp=4;
iPlane = 3;
% suffix = 'raw';
suffix = 'registered';
% suffix = 'rect01_06_registered';

info=ppbox.infoPopulate(subject, expDate, exp);

%% Extract single planes
options.nFramesPerChunk = 4096;
options.iPlane = iPlane;
% ppbox.extractSinglePlane(info, options);

%% Load unregistered data
basename = sprintf('%s_plane%03d_%s', info.basename2p, iPlane, suffix);
% load meta
load(fullfile(info.folderProcessed, basename), 'meta');
if isfield(meta, 'chData')
    fileName = sprintf('%s_%s', meta.chData(1).basename, suffix);
else
    fileName = basename;
end
filePath = fullfile(info.folderProcessed, fileName);
[data, infoProc] = loadArr(filePath);

%% Load stimulus and time information
[~, stimSequence, stimMatrix, ~, samplingRate] = ...
    ssLocal.getStimulusResponseInfo(infoProc);

%% Plot orientation tuning curves for patches
dataPatches = ssLocal.averagePixelsPerPatch(data, numDivs);

gratings.plotTuningOfBrainRegions(dataPatches, samplingRate, stimMatrix, ...
    stimSequence);