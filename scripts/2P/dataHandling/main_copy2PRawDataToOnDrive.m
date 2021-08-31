%% Folders
folderSource1 = '\\ZSERVER.cortexlab.net\Data\Subjects';
folderSource2 = '\\zarchive.cortexlab.net\Data\Subjects';
folderDest = 'C:\STORAGE\OneDrive - University of Sussex\Projects\2021_Vane_Columns\RawData';

%% Database
db = db_columns_vane;

%% Copy
for k = 2:length(db)
    fprintf('Dataset %s %s\n', db(k).subject, db(k).date)
    source = fullfile(folderSource1, db(k).subject, db(k).date, ...
        'metaStructs');
    if ~isfolder(source)
        source = fullfile(folderSource2, db(k).subject, db(k).date, ...
            'metaStructs');
        if ~isfolder(source)
            fprintf('  Meta structures cannot be found.\n')
            continue
        end
    end
    destination = fullfile(folderDest, db(k).subject, db(k).date, ...
        'metaStructs');
    if ~isfolder(destination)
        mkdir(destination)
    end
    experiments = [db(k).expGratings, db(k).expGrayScreen, db(k).expNoise, ...
        db(k).expDark, db(k).expStatic, db(k).expBars];
    for exp = experiments
        d = fullfile(destination, num2str(exp));
        if ~isfolder(d)
            mkdir(d);
        end
        status = copyfile(fullfile(source, num2str(exp), '*_ROI.mat'), d);
        if ~status
            fprintf('  exp %d: data were not successfully copied.\n', exp)
        end
    end
end