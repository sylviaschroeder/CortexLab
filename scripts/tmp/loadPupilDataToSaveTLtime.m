folderROIData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\InfoStructs';
fields = {'expGratings','expGrayScreen'};

db_boutons_driftingGratings_blanks

for k = 1:length(db)
    fprintf('Dataset %d: %s %s\n', k, db(k).subject, db(k).date);
    for st = 1:2 % analyse (1) gratings, (2) blanks
        if isempty(db(k).(fields{st}))
            continue
        end
        folder = fullfile(folderROIData, db(k).subject, ...
            db(k).date, num2str(db(k).(fields{st})));
        file = [db(k).date '_' num2str(db(k).(fields{st})) '_' ...
            db(k).subject  sprintf('_2P_plane%03d_ROI.mat', db(k).planes(1))];
        % load meta
        data=load(fullfile(folder, file));
        meta = data.meta;
        nonVis.loadPupilData(meta);
    end
end