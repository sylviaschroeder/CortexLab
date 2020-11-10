%% Load database
db_blanks
% db_driftingGratings

%% Folders
folderROIData = 'C:\DATA\InfoStructs';
folderCorrRunning = ['C:\RESULTS\nonVisualEffects\maps\correlationRunning_' stimType];
folderCorrPupil = ['C:\RESULTS\nonVisualEffects\maps\correlationPupil_' stimType];
if ~exist(folderCorrRunning, 'dir')
    mkdir(folderCorrRunning)
end
if ~exist(folderCorrPupil, 'dir')
    mkdir(folderCorrPupil)
end

%% Parameters
smoothing = 3; %in sec
filtPoly = 3;

%% Loop across datasets
for k = 1:length(db)
    folder = [folderROIData filesep db(k).subject filesep ...
        db(k).date filesep num2str(db(k).exp)];
    fileStart = [db(k).date '_' num2str(db(k).exp) '_' ...
        db(k).subject];
    file = [fileStart '_2P_plane%03d_ROI.mat'];
    for iPlane=1:length(db(k).planes)
        % load meta
        data=load(fullfile(folder, sprintf(file,db(k).planes(iPlane))));
        meta = data.meta;
        ind = find(meta.ROI.isDuplicate == 0);
        traces = bsxfun(@rdivide, meta.Fcorr(:,ind) - ...
            meta.F0(:,ind), bsxfun(@max, 1, mean(meta.F0(:,ind),1)));
        time = ppbox.getFrameTimes(meta);
        ROImaps = meta.ROI.CellMaps(ind);
        
        % load running and pupil data
        if iPlane == 1
            ballData = nonVis.getRunningSpeed(meta);
            [pupilData, pupilTime] = nonVis.loadPupilData(meta);
        end
        
        if ~isempty(ballData)
            filtWindow = ceil(smoothing / median(diff(ballData.t)));
            if mod(filtWindow,2) == 0
                filtWindow = filtWindow-1;
            end
            running = sgolayfilt(ballData.total, filtPoly, filtWindow);
            results = nonVis.getCorrToNonVisData(traces, time, running, ...
                ballData.t, 'Running speed', [], -1, false);
            handle = nonVis.plotCorrelationMap(results.rho, ROImaps, ...
                [length(meta.validY) length(meta.validX)]);
            title(handle, sprintf('%s %s %d plane %d: Corr. with running',db(k).subject,db(k).date, ...
                db(k).exp,db(k).planes(iPlane)), 'Interpreter', 'none')
            savefig(gcf, fullfile(folderCorrRunning, ...
                sprintf('%s_plane%03d', fileStart, db(k).planes(iPlane))), ...
                'compact');
            close(gcf)
        end
        
        if ~isempty(pupilData)
            pupilTime(length(pupilData.x)+1:end) = [];
            pupilDiam = nonVis.getPupilDiam(pupilData);
            results = nonVis.getCorrToNonVisData(traces, time, pupilDiam, ...
                pupilTime, 'Pupil diam.', [], -1, false);
            handle = nonVis.plotCorrelationMap(results.rho, ROImaps, ...
                [length(meta.validY) length(meta.validX)]);
            title(handle, sprintf('%s %s %d plane %d: Corr. with pupil',db(k).subject,db(k).date, ...
                db(k).exp,db(k).planes(iPlane)), 'Interpreter', 'none')
            savefig(gcf, fullfile(folderCorrPupil, ...
                sprintf('%s_plane%03d', fileStart, db(k).planes(iPlane))), ...
                'compact');
            close(gcf)
        end
    end
end