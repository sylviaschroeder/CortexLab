db_driftingGratings_blank;
% db_boutons_driftingGratings_blanks;

folderROIData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\InfoStructs';
fields = {'expGratings','expGrayScreen','expDark'};

badFrames = NaN(1, length(db));
fprintf('Dataset of %d:', length(db))
for k = 1:length(db)
    fprintf(' %d', k)
    bad = zeros(length(fields), length(db(k).planes));
    total = zeros(length(fields), length(db(k).planes));
    for st = 1:length(fields)
        if ~isfield(db, fields{st}) || isempty(db(k).(fields{st}))
            continue
        end
        folder = fullfile(folderROIData, db(k).subject, ...
            db(k).date, num2str(db(k).(fields{st})));
        fileStart = [db(k).date '_' num2str(db(k).(fields{st})) '_' ...
            db(k).subject];
        file = [fileStart '_2P_plane%03d_ROI.mat'];
        for iPlane = 1:length(db(k).planes)
            % load meta
            data=load(fullfile(folder, sprintf(file,db(k).planes(iPlane))));
            meta = data.meta;
            bad(st,iPlane) = bad(st,iPlane) + sum(all(isnan(meta.F_final),2));
            total(st,iPlane) = total(st,iPlane) + size(meta.F_final,1);
        end
    end
    badFrames(k) = sum(bad(:))/sum(total(:));
end
fprintf('\n')