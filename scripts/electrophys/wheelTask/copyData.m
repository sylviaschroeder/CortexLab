folderTools = 'C:\STORAGE\workspaces';
folderSource = '\\zubjects.cortexlab.net\Subjects';
folderDestination = 'C:\STORAGE\OneDrive - University of Sussex\Projects\2021_Joanna_competition in SC\Data';
folderScript = 'C:\dev\workspace\CortexLab';

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')));
addpath(genpath(fullfile(folderTools, 'alyx-matlab')));
addpath(genpath(fullfile(folderScript)));

%% Copy files
db = db_ephys_task;

for k = 3:25
    source = fullfile(folderSource, db(k).subject, db(k).date, 'alf');
    files = dir(source);
    dest = fullfile(folderDestination, db(k).subject, db(k).date);
    if ~isfolder(dest)
        mkdir(dest)
    end
    for f = 1:length(files)
        if isfolder(fullfile(source, files(f).name))
            continue
        end
        copyfile(fullfile(source, files(f).name), ...
            fullfile(dest, files(f).name));
    end
end