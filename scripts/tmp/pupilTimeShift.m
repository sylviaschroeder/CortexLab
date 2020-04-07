folderROIData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\InfoStructs';
% folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\nonVisualEffects\modelGratingResp\';
folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\boutons\nonVisualEffects\';

fields = {'expGratings','expGrayScreen'};
stimuli = {'gratings', 'grayScreen'};

db_driftingGratings_blank

for k = 1:length(db)
    fprintf('Dataset %d: %s %s\n', k, db(k).subject, ...
        db(k).date);
    for st = 1:2 % analyse (1) gratings, (2) blanks
        if isempty(db(k).(fields{st}))
            continue
        end
        folder = fullfile(folderROIData, db(k).subject, ...
            db(k).date, num2str(db(k).(fields{st})));
        fileStart = [db(k).date '_' num2str(db(k).(fields{st})) '_' ...
            db(k).subject];
        file = [fileStart sprintf('_2P_plane%03d_ROI.mat', db(k).planes(1))];
        % load meta
        data=load(fullfile(folder, sprintf(file,db(k).planes(1))));
        meta = data.meta;
        [~, newTime, oldTime] = nonVis.loadPupilData(meta);
        d = median(oldTime(1:length(newTime))'-newTime);
        fprintf('  %s: %.2f s\n', stimuli{st}, d)
    end
end
