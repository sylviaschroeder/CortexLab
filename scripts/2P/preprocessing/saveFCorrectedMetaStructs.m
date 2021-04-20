%% Load database
db_boutons_driftingGratings_blanks

%% Folders
folderBase = 'C:\STORAGE\OneDrive - University College London\Lab';
folderROIData = fullfile(folderBase, 'DATA\InfoStructs');
d = load(fullfile(folderROIData, 'corrections_boutons.mat'));
corrections = d.corrections;
doCorrect = d.doCorrect;

%% Collect data
experiments = {'expGratings', 'expGrayScreen', 'expDark', 'expNoise'};

for k = 2:length(db)
    fprintf('Dataset %d of %d\n', k, length(db))
    
    for exp = 1:length(experiments)
        if ~isfield(db, experiments{exp}) || isempty(db(k).(experiments{exp}))
            continue
        end
        
        folder = fullfile(folderROIData, db(k).subject, ...
            db(k).date, num2str(db(k).(experiments{exp})));
        file = [sprintf('%s_%d_%s', db(k).date, db(k).(experiments{exp}), ...
            db(k).subject) '_2P_plane%03d_ROI.mat'];
        fileUncorr = [sprintf('%s_%d_%s', db(k).date, db(k).(experiments{exp}), ...
            db(k).subject) '_2P_plane%03d_ROI_uncorrected.mat'];
        for iPlane = 1:length(db(k).planes)
            % load data
            d = load(fullfile(folder, sprintf(file, db(k).planes(iPlane))));
            status = movefile(fullfile(folder, sprintf(file, db(k).planes(iPlane))), ...
                fullfile(folder, sprintf(fileUncorr, db(k).planes(iPlane))));
            if ~status
                return
            end
            meta = d.meta;
            
            F = meta.F_final;
            a = corrections(k).plane(iPlane).a{db(k).(experiments{exp})};
            b = corrections(k).plane(iPlane).b{db(k).(experiments{exp})};
            F = doCorrect(a,b,F);
            meta.F_final = F;
            
            save(fullfile(folder, sprintf(file, db(k).planes(iPlane))), ...
                'meta')
        end
    end
end
