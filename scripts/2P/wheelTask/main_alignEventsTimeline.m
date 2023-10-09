%% Folders
folderTools = 'C:\dev\toolboxes';
folderThisRepo = 'C:\dev\workspaces\CortexLab';

% timelineFolder = '\\ZSERVER.cortexlab.net\Data\expInfo';
% infoStructFolder = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\InfoStructs';

timelineFolder = 'Z:\UCLData\2P_Task';
infoStructFolder = 'Z:\UCLData\2P_Task';

%% Add paths
addpath(genpath(fullfile(folderTools, 'Rigbox')))
addpath(fullfile(folderThisRepo))

%% Build db
db_wheelTask_align;

%% Align times
for iSet = 1:length(db)
    fprintf('Dataset %d of %d: %s %s\n', iSet, length(db), ...
        db(iSet).mouse_name, db(iSet).date)
    for iExp = 1:length(db(iSet).expts)
        fprintf('  Experiment %d (%d)\n', db(iSet).expts(iExp), db(iSet).expNums(iExp))
        data = load(fullfile(timelineFolder, db(iSet).mouse_name, db(iSet).date, ...
            num2str(db(iSet).expts(iExp)), sprintf('%s_%d_%s_Timeline.mat', db(iSet).date, ...
            db(iSet).expts(iExp), db(iSet).mouse_name)));
        timeline = data.Timeline;
        data = load(fullfile(timelineFolder, db(iSet).mouse_name, db(iSet).date, ...
            num2str(db(iSet).expts(iExp)), sprintf('%s_%d_%s_Block.mat', ...
            db(iSet).date, db(iSet).expts(iExp), db(iSet).mouse_name)));
        block = psy.stripIncompleteTrials(data.block);
        
        pdchan = timeline.hw.inputs(strcmp({timeline.hw.inputs.name}, ...
            'photoDiode')).arrayColumn;
        pdchanRaw = timeline.rawDAQData(:,pdchan);
        time = timeline.rawDAQTimestamps;
        
        alignment = preproc.screenSyncAlignment(block, ...
            time, pdchanRaw);
        
        folder = fullfile(infoStructFolder, db(iSet).mouse_name, db(iSet).date, ...
            num2str(db(iSet).expNums(iExp)));
        if ~isdir(folder)
            mkdir(folder)
        end
        save(fullfile(folder, 'timeAlign.mat'), 'alignment');
    end
end