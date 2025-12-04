folders.data = 'C:\Users\Sylvia\OneDrive - University of Sussex\Projects\2023_OrientationColumns\DataToPublish';
folders.dataProcessedEphys = 'Z:\ProcessedData';
folders.tools = 'C:\dev\toolboxes';
folders.repo = 'C:\dev\workspaces\he_schroeder_columns';

addpath(genpath(fullfile(folders.tools, 'npy-matlab')))
addpath(genpath(fullfile(folders.repo)))

subjDirs = dir(fullfile(folders.data, 'ephys'));
subjDirs = subjDirs(~matches({subjDirs.name}, [".",".."]));
for subj = 1:length(subjDirs) % animals
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folders.data, 'ephys', name, '2*'));
    for dt = 1:length(dateDirs) %dates
        date = dateDirs(dt).name;
        f = fullfile(folders.data, 'ephys', name, date);
        f_nas = fullfile(folders.dataProcessedEphys, name, date, 'exp01');

        copyfile(fullfile(f, '_ss_recordings.scChannels.npy'), ...
            fullfile(f_nas, 'recording.scChannels.npy'))
    end
end