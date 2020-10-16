%% Folders
foldermyCode = 'C:\dev\workspace\CortexLab';
folderCode = {genpath('C:\dev\toolboxes\matnwb'), myCode, ...
    fullfile(myCode, 'databases'), fullfile(myCode, 'scripts'), ...
    genpath('C:\STORAGE\workspaces\spikes'), ...
    genpath('C:\STORAGE\workspaces\kilotrodeRig')};
folderBase = '\\ZUBJECTS.cortexlab.net\Subjects';
folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS';
folderSave = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\DataToPublish\arousal_NWB\opticTract';

%% Prepare
for f = 1:length(folderCode)
    addpath(folderCode{f});
end

if ~isfolder(folderSave)
    mkdir(folderSave)
end

%% 