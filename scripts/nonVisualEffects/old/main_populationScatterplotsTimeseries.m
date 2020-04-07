%% Datasets

% db_driftingGratings;
% db_blanks;
db_sparseNoise;

%% Folders
folderROIData = 'C:\Temp2p';
folderScatterplots = ['C:\Temp2p\Results\nonVisualEffects\scatterplots_' stimType];
folderTimeseries = ['C:\Temp2p\Results\nonVisualEffects\timeseries_' stimType];

%% Loop across datasets

for k=1:length(db)
    folder = [folderROIData filesep db(k).subject filesep ...
        db(k).date filesep num2str(db(k).exp)];
    fileStart = [db(k).date '_' num2str(db(k).exp) '_' ...
        db(k).subject];
    file = [fileStart '_2P_plane%03d_ROI.mat'];
    clear info
    for iPlane=1:length(db(k).planes)
        % load meta
        data=load(fullfile(folder, sprintf(file,db(k).planes(iPlane))));
        m = orderfields(data.meta);
        info(iPlane) = m;
    end
    titleLine = [db(k).subject ', ' db(k).date ', exp. ' ...
        num2str(db(k).exp) ', planes ' num2str(db(k).planes)];
    [popResp, popTimes, running, pupilDiam, pupilDiamDeriv, eyeMove] = ...
        nonVis.exploreNonVisualData(info,titleLine,1);
    savefig(gcf, fullfile(folderScatterplots, fileStart), 'compact');
%     close gcf
    
    
    figure('Position', [15 520 1900 580])
    zPop = (popResp - min(popResp)) ./ (max(popResp)-min(popResp));
    plot(popTimes, zPop, 'k')
    yTicks = 0;
    labels = {'Ca response'};
    
    hold on
    if ~isempty(running)
        zRunning = (running - min(running)) ./ (max(running)-min(running));
        yTicks(end+1) = yTicks(end)-1;
        labels{end+1} = 'Running';
        plot(popTimes, zRunning+yTicks(end), 'r')
    end
    if ~isempty(pupilDiam)
        zDiam = (pupilDiam - min(pupilDiam)) ./ (max(pupilDiam)-min(pupilDiam));
        yTicks(end+1) = yTicks(end)-1;
        labels{end+1} = 'Pupil diam.';
        plot(popTimes, zDiam+yTicks(end), 'b')
    end
    if ~isempty(pupilDiamDeriv)
        zDeriv = (pupilDiamDeriv - min(pupilDiamDeriv)) ./ (max(pupilDiamDeriv)-min(pupilDiamDeriv));
        yTicks(end+1) = yTicks(end)-1;
        labels{end+1} = 'Pupil diam. deriv.';
        plot(popTimes, zDeriv+yTicks(end), 'g')
    end
    if ~isempty(eyeMove)
        zEye = (eyeMove - min(eyeMove)) ./ (max(eyeMove)-min(eyeMove));
        yTicks(end+1) = yTicks(end)-1;
        labels{end+1} = 'Eye movement';
        plot(popTimes, zEye+yTicks(end), 'c')
    end
    xlim(popTimes([1 end]))
    xlabel('Time (in s)')
    title(titleLine, 'Interpreter', 'none')
    set(gca, 'YTick', flip(yTicks), 'YTickLabel', flip(labels))
    savefig(gcf, fullfile(folderTimeseries, fileStart), 'compact');
end