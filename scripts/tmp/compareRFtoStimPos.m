%% Folders and parameters
timelineFolder = '\\ZSERVER.cortexlab.net\Data\expInfo';
RFfolder = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\wheelTask\receptiveFields';

%% Datasets
make_db_wheelTask_visualNoise

for k = 1:length(db)
    fprintf('Dataset %d: %s %s\n', k, db(k).mouse_name, db(k).date);
    
    % load block structure
    folder = fullfile(timelineFolder, db(k).mouse_name(9:end), db(k).date, ...
        num2str(db(k).taskExp));
    file = dir(fullfile(folder, '*_Block.mat'));
    data = load(fullfile(folder, file.name));
    block = data.block;
    stimPos = [block.parameters.distBetweenTargets/2, ...
        block.parameters.targetAltitude];
    
    % load pixel RF
    folder = fullfile(RFfolder, sprintf('%s_%s_%d', ...
        db(k).mouse_name(9:end), db(k).date, db(k).newExp));
    data = load(fullfile(folder, 'RFinfo.mat'));
    RFinfo = data.RFinfo;
    
    RFdistance = [abs(RFinfo(1).RF_pos(1))-stimPos(1), ...
        -RFinfo(1).RF_pos(2)-stimPos(2)] ./ block.parameters.cueSigma';
    % if x<0, RF is more peripheral than stimulus; if y<0, RF is below
    % stimulus
    
    % save RF distance
    save(fullfile(folder, 'RFdist.mat'), 'RFdistance')
end