%% Dataset
% subject = 'SS072';
% date = '2016-12-04';
% eyeExp = 1;
subject = 'SS074';
date = '2017-01-05';
eyeExp = 2;

%% Folders
folderEye = '\\ZSERVER.cortexlab.net\Data\EyeCamera';
folderAlign = '\\zubjects.cortexlab.net\Subjects';
folderTools = 'C:\STORAGE\workspaces';

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))

%% Align times
% load time of eye movie frames aligned to "timeline"
d = load(fullfile(folderEye, subject, date, num2str(eyeExp), ...
    'eye_timeStamps.mat'));
t_timeline = d.tVid';

% load alignment from "timeline" to ephys (master)
tlToEphys = readNPY(fullfile(folderAlign, subject, date, 'alignments', ...
    sprintf('correct_timeline_%d_to_ephys_SC.npy', eyeExp)));

% do alignment
t_ephys = t_timeline .* tlToEphys(1) + tlToEphys(2);

%% Save aligned times
save(fullfile(folderEye, subject, date, num2str(eyeExp), 'eyeTime.mat'), ...
    't_ephys')