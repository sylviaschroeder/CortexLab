%% Folders
folderROIData = 'C:\DATA\InfoStructs';
folderResults = 'C:\RESULTS\nonVisualEffects\modelGratingResp\';

%% Parameters
nonVisualSignal = 'running'; % 'running' or 'pupil'
% smoothing of signals (calcium traces and nonvisual signal)
smoothing = 3; %in sec
filtPoly = 3;

minR2 = 0.3;

%% Load data
% kernel fitting results
data = load(fullfile(folderResults, nonVisualSignal, ...
    'modelResults_baselineSubtracted_moreModels.mat'));
results = data.results;
% tuning information
data = load(fullfile(folderResults, nonVisualSignal, 'tuning.mat'));
tuning = data.tuning;

%% Correlation with nonvisual signal
for iSet = 1:length(results)
    metaFolder = fullfile(folderROIData, results(iSet).subject, ...
        results(iSet).date, num2str(results(iSet).exp));
    fileStart = [results(iSet).date '_' num2str(results(iSet).exp) '_' ...
        results(iSet).subject];
    file = [fileStart '_2P_plane%03d_ROI.mat'];
    
    for iPlane = 1:length(results(iSet).plane)
        % load meta
        data = load(fullfile(metaFolder, sprintf(file, ...
            results(iSet).planes(iPlane))));
        meta = data.meta;
        
        if iPlane == 1
            switch nonVisualSignal
                case 'running'
                    ballData = nonVis.getRunningSpeed(meta);
                    if isempty(ballData)
                        continue
                    end
                    filtWindow = ceil(smoothing / median(diff(ballData.t)));
                    if mod(filtWindow,2) == 0
                        filtWindow = filtWindow-1;
                    end
                    nonVisData = sgolayfilt(ballData.total, filtPoly, filtWindow);
                    nonVisTime = ballData.t;
                case 'pupil'
                    [pupilData, nonVisTime] = nonVis.loadPupilData(meta);
                    if isempty(pupilData)
                        continue
                    end
                    nonVisTime(length(pupilData.x)+1:end) = [];
                    nonVisData = nonVis.getPupilDiam(pupilData);
            end
        end
        
        IDs = results(iSet).plane(iPlane).cellIDs;
        traces = bsxfun(@rdivide, meta.Fcorr(:,IDs)-meta.F0(:,IDs), ...
            bsxfun(@max, 1, mean(meta.F0(:,IDs),1)));
        time = ppbox.getFrameTimes(meta);
        filtWindow = ceil(smoothing / median(diff(time)));
        if mod(filtWindow,2) == 0
            filtWindow = filtWindow-1;
        end
        traces = sgolayfilt(traces, filtPoly, filtWindow);
        
        nonVisInt = interp1(nonVisTime, nonVisData, time, 'pchip');
        
        [r,p] = corr(traces, nonVisInt(:));
        
        rho(iSet).plane(iPlane).corr = r;
    end
end

%% Compare correlation to tuning
allRho = [];
allSI = [];
for iSet = 1:length(tuning)
    for iPlane = 1:length(tuning(iSet).plane)
        allRho = [allRho; rho(iSet).plane(iPlane).corr(:)];
        allSI = [allSI; max(tuning(iSet).plane(iPlane).stimOnly(:,[2 4]),[],2)];
    end
end

indNaN = isnan(allSI) | allOSI<0.3;
figure
boxplot(allRho, indNaN+1)
set(gca,'XTickLabel',{'not tuned','tuned'})
ylabel(sprintf('Correlation with %s', nonVisualSignal))

indNaN = isnan(allSI);
figure
plot(allSI(~indNaN), allRho(~indNaN), 'k.')
xlabel('Tuning selectivity (max (OSI,DSI))')
ylabel(sprintf('Correlation with %s', nonVisualSignal))