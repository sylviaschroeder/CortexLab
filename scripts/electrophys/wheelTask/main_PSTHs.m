%% Folders
folderData = '\\zubjects.cortexlab.net\Subjects';

folderScript = 'C:\dev\workspace\CortexLab';
folderTools = 'C:\STORAGE\workspaces';

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')));
addpath(genpath(fullfile(folderTools, 'spikes')));
addpath(genpath(fullfile(folderTools, 'wheelAnalysis')));
addpath(genpath(fullfile(folderScript)));

%% Define datasets
subject = 'SS093';
date = '2018-05-24';

folderAlf = fullfile(folderData, subject, date, 'alf');

%% Collect relevant event times etc for data set
% stim times, movement times, go cue, ...
taskData = task.getTaskInfo(folderAlf, {'wheelPos'});
[moveOn, moveOff, moveDispl, peakVel, peakAmp] = wheel.findWheelMoves3();

%% Plot PSTH for each neuron
% For each neuron, rasters for stim onset (separated by stim ID), go cue,
% movement onset (separated by before/after go cue, left/right), feedback 
% (separated by pos/neg), stim offset (separated by stim ID)
