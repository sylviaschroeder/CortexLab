%% Folders
sourceFolder = '\\ZSERVER.cortexlab.net\Data\EyeCamera';
destFolder = 'J:\EyeData';

%% Data base
db_wheelTask;

%% Copy all eye movies
for k = 1:length(db)
    fprintf('%s %s:', db(k).subject, db(k).date)
    for exp = 1:length(db(k).exp)
        source = fullfile(sourceFolder, db(k).subject, db(k).date, ...
            num2str(db(k).exp(exp)));
        dest = fullfile(destFolder, db(k).subject, db(k).date, ...
            num2str(db(k).exp(exp)));
        file = sprintf('%s_%d_%s_eye.mj2', ...
            db(k).date, db(k).exp(exp), db(k).subject);
        if ~isfolder(dest)
            mkdir(dest)
        end
        status = copyfile(fullfile(source, file), fullfile(dest, file));
        fprintf(' exp %d: %d;', db(k).exp(exp), status)
    end
    fprintf('\n')
end