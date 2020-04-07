folderROIData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\InfoStructs';

for iExp = 11:length(results)
    if isempty(results(iExp).expGratings)
        continue
    end
    fprintf('Dataset %d: %s %s exp.: %d\n', iExp, results(iExp).subject, ...
        results(iExp).date, results(iExp).expGratings);
    folder = [folderROIData filesep results(iExp).subject filesep ...
        results(iExp).date filesep num2str(results(iExp).expGratings)];
    fileStart = [results(iExp).date '_' num2str(results(iExp).expGratings) '_' ...
        results(iExp).subject];
    file = [fileStart '_2P_plane%03d_ROI.mat'];
    iPlane = 1;
    % load meta
    data=load(fullfile(folder, sprintf(file,results(iExp).planes(iPlane))));
    meta = data.meta;
    meta.folderTL = strrep(meta.folderTL, 'ioo.ucl.ac.uk', ...
        'cortexlab.net');
    nonVis.loadPupilData(meta);
end