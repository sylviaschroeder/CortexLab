%% Folders
sourceFolder = '\\ZSERVER.cortexlab.net\Data\EyeCamera';
sourceFolder2 = '\\zubjects.cortexlab.net\Subjects';
destFolder = 'J:\EyeData\videos';

folderScript = 'C:\dev\workspace\CortexLab';

addpath(genpath(fullfile(folderScript)));

%% Data base
% db_wheelTask;
db = db_ephys_task;

%% Copy all eye movies
% % A) for 2P data
% for k = 1:length(db)
%     fprintf('%s %s:', db(k).subject, db(k).date)
%     for exp = 1:length(db(k).exp)
%         source = fullfile(sourceFolder, db(k).subject, db(k).date, ...
%             num2str(db(k).exp(exp)));
%         dest = fullfile(destFolder, db(k).subject, db(k).date, ...
%             num2str(db(k).exp(exp)));
%         file = sprintf('%s_%d_%s_eye.mj2', ...
%             db(k).date, db(k).exp(exp), db(k).subject);
%         if ~isfolder(dest)
%             mkdir(dest)
%         end
%         status = copyfile(fullfile(source, file), fullfile(dest, file));
%         fprintf(' exp %d: %d;', db(k).exp(exp), status)
%     end
%     fprintf('\n')
% end


% B) for ephys data
for k = 2:length(db)
    fprintf('%s %s exp %d:', db(k).subject, db(k).date, db(k).expTL)
    source = fullfile(sourceFolder, db(k).subject, db(k).date, ...
        num2str(db(k).expTL));
    if ~isfolder(source)
        source = fullfile(sourceFolder2, db(k).subject, db(k).date, ...
        num2str(db(k).expTL));
    end
    file = sprintf('%s_%d_%s_eye.mj2', ...
        db(k).date, db(k).expTL, db(k).subject);
    status = copyfile(fullfile(source, 'eye.mj2'), fullfile(destFolder, file));
    fprintf(' %d\n', status)
end