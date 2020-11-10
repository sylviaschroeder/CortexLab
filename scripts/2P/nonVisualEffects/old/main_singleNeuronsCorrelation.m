%% Load database
% db_blanks
% db_driftingGratings
db_sparseNoise

%% Folders
folderROIData = 'C:\Temp2p';
folderCorrRunning = ['C:\Temp2p\Results\nonVisualEffects\correlationRunning_' stimType];
folderCorrPupil = ['C:\Temp2p\Results\nonVisualEffects\correlationPupil_' stimType];

%% Parameters
smoothing = 3; %in sec
filtPoly = 3;

%% Loop across datasets

for k=1:length(db)
    folder = [folderROIData filesep db(k).subject filesep ...
        db(k).date filesep num2str(db(k).exp)];
    fileStart = [db(k).date '_' num2str(db(k).exp) '_' ...
        db(k).subject];
    file = [fileStart '_2P_plane%03d_ROI.mat'];
    neurons = cell(1, length(db(k).planes));
    for iPlane=1:length(db(k).planes)
        % load meta
        data=load(fullfile(folder, sprintf(file,db(k).planes(iPlane))));
        if iPlane == 1;
            info = orderfields(data.meta);
        else
            info(iPlane) = orderfields(data.meta);
        end
        ind = find(strcmp(data.meta.ROI.CellClasses, 's'));
        neurons{iPlane} = ind;
    end
    % interpolate neural responses to get a single series of time points
    [calciumTraces, calciumTime] = ssLocal.matchSampleTimesAcrossPlanes(infos);
    cacliumTraces = cat(2, cacliumTraces{:});
    neuronIDs = [];
    for iPlane = 1:length(db(k).planes)
        neuronIDs = [neuronIDs; [ones(length(neurons{iPlane}),1) * ...
            db(k).planes(iPlane), neurons{iPlane}']];
    end
%     calciumTraces = medfilt1(calciumTraces, nSamples);
    filtWindow = ceil(smoothing / median(diff(calciumTime)));
    if mod(filtWindow,2) == 0
        filtWindow = filtWindow-1;
    end
    calciumTraces = sgolayfilt(calciumTraces, filtPoly, filtWindow);
    
    ballData = nonVis.getRunningSpeed(data.meta);
    if ~isempty(ballData)
        filtWindow = ceil(smoothing / median(diff(ballData.t)));
        if mod(filtWindow,2) == 0
            filtWindow = filtWindow-1;
        end
%         running = medfilt1(ballData.total, nSamples);
        running = sgolayfilt(ballData.total, filtPoly, filtWindow);
        [runResults, figHandles] = nonVis.getCorrToNonVisData(calciumTraces, ...
            calciumTime, running, ballData.t, 'Running Speed', ...
            neuronIDs, .3, 1);
        for iFig = 1:length(figHandles.traces)
            savefig(figHandles.traces(iFig), fullfile(folderCorrRunning, ...
                [fileStart sprintf('_traces%02d',iFig)]), 'compact');
        end
        savefig(figHandles.rho, fullfile(folderCorrRunning, ...
            [fileStart '_rhos']), 'compact');
    end
    
    [pupilData, pupilTime] = nonVis.loadPupilData(data.meta);
    if ~isempty(pupilData)
        pupilTime(length(pupilData.x)+1:end) = [];
        filtWindow = ceil(smoothing / median(diff(pupilTime)));
        if mod(filtWindow,2) == 0
            filtWindow = filtWindow-1;
        end
        pupilDiam = sqrt(4 * pupilData.area / pi);
%         pupilDiam = medfilt1(pupilDiam, nSamples);
        pupilDiam = sgolayfilt(pupilDiam, filtPoly, filtWindow);
        pupilDiam(pupilData.blink | ~pupilData.goodFit) = NaN;
        [pupilResults, figHandles] = nonVis.getCorrToNonVisData(calciumTraces, ...
            calciumTime, pupilDiam, pupilTime, 'Pupil Diam.', ...
            neuronIDs, .3, 1);
        for iFig = 1:length(figHandles.traces)
            savefig(figHandles.traces(iFig), fullfile(folderCorrPupil, ...
                [fileStart sprintf('_traces%02d',iFig)]), 'compact');
        end
        savefig(figHandles.rho, fullfile(folderCorrPupil, ...
            [fileStart '_rhos']), 'compact');
    end
end