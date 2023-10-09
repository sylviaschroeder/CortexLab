%% Folders
folderBase = 'Z:\UCLData\2P_Task';
folderTools = 'C:\dev\toolboxes';
folderThisRepo = 'C:\dev\workspaces\CortexLab';

%% Add paths
addpath(fullfile(folderThisRepo))

%% Build db
make_db_wheelTask2023;

%% Get time of eye video frames
for iSet = 2:length(db)
    for iExp = 1:length(db(iSet).expts)
        fprintf('%s %s exp %d\n', db(iSet).mouse_name, db(iSet).date, ...
            db(iSet).expts(iExp))
        folder = fullfile(folderBase, db(iSet).mouse_name, db(iSet).date, ...
            num2str(db(iSet).expts(iExp)));
        file = dir(fullfile(folder, sprintf('*_%s_eye.mj2', db(iSet).mouse_name)));
        video = fullfile(folder, file.name(1:end-4));
        file = dir(fullfile(folder, sprintf('*_%s_Timeline.mat', db(iSet).mouse_name)));
        tl = fullfile(folder, file.name);
        pupilTimes = eye.getFrameTimes(video, tl);
        save([video '_TLtime.mat'], "pupilTimes");
    end
end