units = 'boutons';
% units = 'neurons';

timeGap = 600; %in s, gap between experiments

experiments = {'expGratings', 'expGrayScreen', 'expDark', 'expNoise'};
expNames = {'gratings','grayScreen','darkness','sparseNoise'};

%% Folders
foldermyCode = 'C:\dev\workspace\CortexLab';
folderCode = {genpath('C:\dev\toolboxes\matnwb'), myCode, ...
    fullfile(myCode, 'databases'), fullfile(myCode, 'scripts'), ...
    genpath('C:\STORAGE\workspaces\spikes'), ...
    genpath('C:\STORAGE\workspaces\kilotrodeRig')};
folderBase = 'C:\STORAGE\OneDrive - University College London\Lab';
folderROIData = fullfile(folderBase, 'DATA\InfoStructs');
if strcmp(units, 'boutons')
    folderTuning = fullfile(folderBase, 'RESULTS\boutons\nonVisualEffects');
    folderRF = fullfile(folderBase, 'RESULTS\receptiveFields\boutons');
    d = load(fullfile(folderROIData, 'corrections_boutons.mat'));
    corrections = d.corrections;
    doCorrect = d.doCorrect;
    folderSave = fullfile(folderBase, 'DATA\DataToPublish\arousal_NWB\boutons');
else
    folderTuning = fullfile(folderBase, 'RESULTS\nonVisualEffects\modelGratingResp');
    folderRF = fullfile(folderBase, 'RESULTS\receptiveFields\SC neurons');
    corrections = [];
    folderSave = fullfile(folderBase, 'DATA\DataToPublish\arousal_NWB\sc neurons 2p');
end

%% Prepare
% add paths
for f = 1:length(folderCode)
    addpath(folderCode{f});
end
% make folder to save data
if ~isfolder(folderSave)
    mkdir(folderSave)
end
% load database
if strcmp(units, 'boutons')
    db_boutons_driftingGratings_blanks
else
    db_driftingGratings_blank
end

k = 1;

%% General information
subject = db(k).subject;
ind = strfind(subject, 'SS');
subject = subject(ind:end);
date = db(k).date;
folderSubject = fullfile(folderSave, subject);

nwb = NwbFile( ...
    'session_description', '2P imaging of retinal boutons',...
    'identifier', [subject '_' date], ...
    'session_start_time', datetime(date, 'InputFormat', 'yyyy-MM-dd'), ...
    'general_experimenter', 'Sylvia Schroeder', ...
    'general_institution', 'University College London', ...
    'general_related_publications', 'DOI:10.1016/j.neuron.2020.04.026');

%% Imaging data
