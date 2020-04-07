timelineFolder = '\\ZSERVER.cortexlab.net\Data\expInfo';
infoStructFolder = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\InfoStructs';

db_wheelTask_align;

for iSet = 1:length(db)
    fprintf('Dataset %d of %d: %s %s\n', iSet, length(db), ...
        db(iSet).mouse_name, db(iSet).date)
    for iExp = 1:length(db(iSet).expts)
        fprintf('  Experiment %d (%d)\n', db(iSet).expts(iExp), db(iSet).expNums(iExp))
        data = load(fullfile(timelineFolder, db(iSet).mouse_name, db(iSet).date, ...
            num2str(db(iSet).expts(iExp)), sprintf('%s_%d_%s_Timeline.mat', db(iSet).date, ...
            db(iSet).expts(iExp), db(iSet).mouse_name)));
        timeline = data.Timeline;
        block = psy.stripIncompleteTrials(...
            loadVar(dat.expFilePath(db(iSet).mouse_name, db(iSet).date, ...
            db(iSet).expts(iExp), 'block', 'master'), 'block'));
        
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