%% Load database
% db_blanks
db_driftingGratings
% db_sparseNoise

%% Folders
folderROIData = 'C:\DATA\InfoStructs';
folderCorrRunning = ['C:\RESULTS\nonVisualEffects\correlationRunning_' stimType];
folderCorrPupilDiam = ['C:\RESULTS\nonVisualEffects\correlationPupil_' ...
    stimType '\pupilDiameter'];
folderCorrPupilMov = ['C:\RESULTS\nonVisualEffects\correlationPupil_' ...
    stimType '\pupilMovement'];
folderCorrNonVis = 'C:\RESULTS\nonVisualEffects\correlationNonVisualSignals';
if ~exist(folderCorrRunning, 'dir')
    mkdir(folderCorrRunning)
end
if ~exist(folderCorrPupilDiam, 'dir')
    mkdir(folderCorrPupilDiam)
end
if ~exist(folderCorrPupilMov, 'dir')
    mkdir(folderCorrPupilMov)
end

if ~exist(folderCorrNonVis, 'dir')
    mkdir(folderCorrNonVis)
end

%% Parameters
smoothStd = .25; %in sec

%% Loop across datasets

for k=1:length(db)
    folder = [folderROIData filesep db(k).subject filesep ...
        db(k).date filesep num2str(db(k).exp)];
    fileStart = [db(k).date '_' num2str(db(k).exp) '_' ...
        db(k).subject];
    file = [fileStart '_2P_plane%03d_ROI.mat'];
    neurons = cell(1, length(db(k).planes));
    times = cell(1, length(db(k).planes));
    traces = cell(1, length(db(k).planes));
    isGad = cell(1, length(db(k).planes));
    for iPlane=1:length(db(k).planes)
        % load meta
        data=load(fullfile(folder, sprintf(file,db(k).planes(iPlane))));
        ind = find(data.meta.ROI.isDuplicate == 0 & data.meta.ROI.isSwitchOn == 0);
        neurons{iPlane} = ind;
        times{iPlane} = ppbox.getFrameTimes(data.meta);
        traces{iPlane} = data.meta.F_final(:,ind);
        isGad{iPlane} = data.meta.ROI.isGad(ind);
    end
    % interpolate neural responses to get a single series of time points
    calciumTime = times{1};
    calciumTraces = cell(1, length(db(k).planes));
    calciumTraces{1} = traces{1};
    for iPlane = 2:length(db(k).planes)
        calciumTraces{iPlane} = interp1(times{iPlane}, ...
            traces{iPlane}, calciumTime, 'pchip');
    end
    calciumTraces = cat(2, calciumTraces{:});
    neuronIDs = [];
    for iPlane = 1:length(db(k).planes)
        neuronIDs = [neuronIDs; [ones(length(neurons{iPlane}),1) * ...
            db(k).planes(iPlane), neurons{iPlane}]];
    end
    isGad = cat(1, isGad{:});
%     calciumTraces = medfilt1(calciumTraces, nSamples);
%     filtWindow = ceil(smoothing / median(diff(calciumTime)));
%     if mod(filtWindow,2) == 0
%         filtWindow = filtWindow-1;
%     end
%     calciumTraces = sgolayfilt(calciumTraces, filtPoly, filtWindow);
    
    ballData = nonVis.getRunningSpeed(data.meta);
    if ~isempty(ballData)
        stdSamples = round(smoothStd / median(diff(ballData.t)));
        convWindow = normpdf(-3*stdSamples:3*stdSamples, 0, stdSamples);
        running = conv(ballData.total, convWindow, 'same');
        running = interp1(ballData.t,running, calciumTime);
        [runResults, figHandles] = nonVis.getCorrToNonVisData(calciumTraces, ...
            calciumTime, running, calciumTime, 'running Speed', ...
            neuronIDs, -1, 1, [], isGad);
        for iFig = 1:length(figHandles.traces)
            savefig(figHandles.traces(iFig), fullfile(folderCorrRunning, ...
                [fileStart sprintf('_traces%02d',iFig)]), 'compact');
            close(figHandles.traces(iFig))
        end
%         savefig(figHandles.rho, fullfile(folderCorrRunning, ...
%             [fileStart '_rhos']), 'compact');
%         close(figHandles.rho)
    end
    
    [pupilData, pupilTime] = nonVis.loadPupilData(data.meta);
    if ~isempty(pupilData)
        pupilTime(length(pupilData.x)+1:end) = [];
        [pupilDiam,pupilMov] = nonVis.getPupilDiam(pupilData); % pupil diameter
        pupilMov = abs([0; diff(pupilMov)]);
        
        ind = isnan(pupilDiam);
        indInterp = hist(pupilTime(ind), calciumTime) > 0;
        pupilDiam = interp1(pupilTime, pupilDiam, calciumTime);
        pupilDiam = sgolayfilt(pupilDiam, filtPoly, filtWindow);
        pupilDiam(indInterp) = NaN;
        [pupilDiamResults, diamfigHandles] = nonVis.getCorrToNonVisData(calciumTraces, ...
            calciumTime, pupilDiam, calciumTime, 'pupil diam.', ...
            neuronIDs, -1, 1, [], isGad);
        for iFig = 1:length(diamfigHandles.traces)
            savefig(diamfigHandles.traces(iFig), fullfile(folderCorrPupilDiam, ...
                [fileStart sprintf('_traces%02d',iFig)]), 'compact');
            close(diamfigHandles.traces(iFig))
        end
        savefig(diamfigHandles.rho, fullfile(folderCorrPupilDiam, ...
            [fileStart '_rhos']), 'compact');
        close(diamfigHandles.rho)
        
        ind = isnan(pupilMov);
        indInterp = hist(pupilTime(ind), calciumTime) > 0;
        pupilMov = interp1(pupilTime, pupilMov, calciumTime);
        pupilMov = sgolayfilt(pupilMov, filtPoly, filtWindow);
        pupilDiam(indInterp) = NaN;
        [pupilMovResults, movfigHandles] = nonVis.getCorrToNonVisData(calciumTraces, ...
            calciumTime, pupilMov, calciumTime, 'pupil movem.', ...
            neuronIDs, -1, 1, [], isGad);
        for iFig = 1:length(movfigHandles.traces)
            savefig(movfigHandles.traces(iFig), fullfile(folderCorrPupilMov, ...
                [fileStart sprintf('_traces%02d',iFig)]), 'compact');
            close(movfigHandles.traces(iFig))
        end
        savefig(movfigHandles.rho, fullfile(folderCorrPupilMov, ...
            [fileStart '_rhos']), 'compact');
        close(movfigHandles.rho)
    end
    
    signals = [running(:), pupilDiam(:), pupilMov(:)];
    signalsNorm = bsxfun(@rdivide, bsxfun(@minus, signals, ...
        nanmean(signals,1)), nanstd(signals,0,1));
    signalsNorm = bsxfun(@minus, signalsNorm, [0 3 6]);
    ind = all(~isnan(signals),2);
    r = corr(signals(ind,:));
    figure('Position', [680 42 1060 1074])
    subplot(3,1,1)
    plot(calciumTime, signalsNorm)
    xlim(calciumTime([1 end]))
    xlabel('Time (in s)')
    legend('running', 'pupil diam.', 'pupil movem.')
    subplot(3,2,3)
    plot(running, pupilDiam, 'k.')
    xlabel('Running')
    ylabel('Pupil diameter')
    title(sprintf('Corr.: %.3f',r(1,2)))
    subplot(3,2,4)
    plot(running, pupilMov, 'k.')
    xlabel('Running')
    ylabel('Pupil movements')
    title(sprintf('Corr.: %.3f',r(1,3)))
    subplot(3,2,5)
    plot(pupilDiam, pupilMov, 'k.')
    xlabel('Pupil diameter')
    ylabel('Pupil movements')
    title(sprintf('Corr.: %.3f',r(2,3)))
    savefig(gcf, fullfile(folderCorrNonVis, fileStart), 'compact');
    close(gcf)
end