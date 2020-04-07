%% Load database
db_driftingGratings

%% Folders
folderROIData = 'C:\Temp2p';
folderOrientationMaps = 'C:\Temp2p\Results\orientation\maps';
folderTuningModelRunning = ['C:\Temp2p\Results\nonVisualEffects\tuningModelledRunning_' stimType];

%% Parameters
% to determine F_0
window = 400;
percentile = 5;
% to extract stimulus responses
timeLimits = [0 0]; % in sec; (1) before stim onset, (2) after stim offset
% to calculate interneural distances
fieldOfView = [2 500 520; 2.2 450 490; 2.5 400 430; 3 335 370];
% to evaluate fit of tuning curve
minRSquared = 0.2;
% to evaluate orientation/direction vectors (tuning strength)
minVectorLength = 0.3;

%% Loop across datasets
preferredOrientations = [];
preferredDirections = [];
% preferredOrientationsBoth = [];
orientationDifferences = [];
directionDifferences = [];
% orientationBothDifferences = [];
distances = [];

respNeurons = [];       % 1 if ANOVA (on responses: [stimulus x trial]) was significant
tunedNeurons = [];    % 1 if rSquared of fit to tuning curve was >= 0.2
for k=1:length(db)
    
    %% Load data
    folder = [folderROIData filesep db(k).subject filesep ...
        db(k).date filesep num2str(db(k).exp)];
    fileStart = [db(k).date '_' num2str(db(k).exp) '_' ...
        db(k).subject];
    file = [fileStart '_2P_plane%03d_ROI.mat'];

    if ~isdir(fullfile(folderOrientationMaps, [db(k).subject '_' ...
            db(k).date '_' num2str(db(k).exp)]))
        mkdir(fullfile(folderOrientationMaps, [db(k).subject '_' ...
            db(k).date '_' num2str(db(k).exp)]));
    end
    % load results + respondingNeurons from fitting tuning curve and running analysis
    load(fullfile(folderTuningModelRunning, [db(k).subject '_' ...
        db(k).date '_' num2str(db(k).exp)], 'results.mat'), 'results', ...
        'respondingNeurons')
    [~,order] = sort(results.corr.order);
    respNeurons = [respNeurons; respondingNeurons(order)];
    rSquared = results.errors(order,1);
    rSquared(~respondingNeurons(order)) = 0;
    tunedNeurons = [tunedNeurons; rSquared >= minRSquared];
    pars = results.tuningPars(order,1,1);

    prefOrisDataset = [];
    prefDirsDataset = [];
%     prefOrisBothDataset = [];
    numCells = 0;
    for iPlane=1:length(db(k).planes)
        % load meta
        data=load(fullfile(folder, sprintf(file,db(k).planes(iPlane))));
        info = data.meta;
        ind = find(strcmp(data.meta.ROI.CellClasses, 's'));
        neurons = ind;
    
        %% Preprocess calcium traces
        % interpolate neural responses to get a single series of time points
%         calciumTraces = info.F(:,neurons);
%         calciumTime = ppbox.getFrameTimes(info);
%         samplingRate = 1 / median(diff(calciumTime));
%         % caclulate F_0
%         [~,F_0] = ssLocal.removeSlowDrift(calciumTraces, calciumTime, window, percentile);
%         calciumTraces = (calciumTraces - F_0) ./ F_0;
    
        %% Get stimulus information
%         stimTimes = ppbox.getStimTimes(info);
%         stimSequence = ppbox.getStimSequence(info);
%         stimMatrix = ppbox.buildStimMatrix(stimSequence, stimTimes, calciumTime);
%         [orientations, blank] = gratings.getOrientations(stimSequence);
    
        %% Get stimulus pref., maps, pref. differences, and distances
        prefOris = mod(pars((1:length(neurons))+numCells),180);
        prefDirs = mod(pars((1:length(neurons))+numCells),360);
        
%         [prefOris, prefDirs, prefOrisBoth] = ...
%             gratings.getPreferredOrientations(calciumTraces, samplingRate, ...
%             stimMatrix, orientations, blank, timeLimits(1), timeLimits(2));
        prefOrisDataset = [prefOrisDataset; prefOris];
        prefDirsDataset = [prefDirsDataset; prefDirs];
%         prefOrisBothDataset = [prefOrisBothDataset; prefOrisBoth];
%         figure
%         hist(prefOris(prefOris(:,2)>=.3,1), 5:10:180)
%         xlabel('Preferred orientation')
%         ylabel('# Neurons')
%         savefig(gcf, fullfile(folderOrientationMaps, ...
%             [db(k).subject '_' db(k).date '_' num2str(db(k).exp)], ...
%             sprintf('preferredOrientations_plane%02d', db(k).planes(iPlane))), ...
%             'compact');
%         saveas(gcf, fullfile(folderOrientationMaps, ...
%             [db(k).subject '_' db(k).date '_' num2str(db(k).exp)], ...
%             sprintf('preferredOrientations_plane%02d.jpg', ...
%             db(k).planes(iPlane))));
%         close(gcf)
%         figure
%         hist(prefDirs(prefDirs(:,2)>=.3,1), 5:10:360)
%         xlabel('Preferred direction')
%         ylabel('# Neurons')
%         savefig(gcf, fullfile(folderOrientationMaps, ...
%             [db(k).subject '_' db(k).date '_' num2str(db(k).exp)], ...
%             sprintf('preferredDirections_plane%02d', db(k).planes(iPlane))), ...
%             'compact');
%         saveas(gcf, fullfile(folderOrientationMaps, ...
%             [db(k).subject '_' db(k).date '_' num2str(db(k).exp)], ...
%             sprintf('preferredDirections_plane%02d.jpg', ...
%             db(k).planes(iPlane))));
%         close(gcf)
%         figure
%         hist(prefOrisBoth(prefOrisBoth(:,2)>=.3,1), 5:10:180)
%         xlabel('Preferred orientation')
%         ylabel('# Neurons')
%         savefig(gcf, fullfile(folderOrientationMaps, ...
%             [db(k).subject '_' db(k).date '_' num2str(db(k).exp)], ...
%             sprintf('preferredOrientationsCombi_plane%02d', db(k).planes(iPlane))), ...
%             'compact');
%         saveas(gcf, fullfile(folderOrientationMaps, ...
%             [db(k).subject '_' db(k).date '_' num2str(db(k).exp)], ...
%             sprintf('preferredOrientationsCombi_plane%02d.jpg', ...
%             db(k).planes(iPlane))));
%         close(gcf)
        
        figHandles(1) = gratings.plotOrientationMap(...
            [prefOris, rSquared((1:length(neurons))+numCells)], 'ori', ...
            info.ROI.CellMaps(neurons), ...
            [length(info.validY) length(info.validX)], minRSquared);
        figHandles(2) = gratings.plotOrientationMap(...
            [prefDirs, rSquared((1:length(neurons))+numCells)], 'dir', ...
            info.ROI.CellMaps(neurons), ...
            [length(info.validY) length(info.validX)], minRSquared);
        titles = {'orientationMap','directionMap'};
        for iFig = 1:2
            savefig(figHandles(iFig), fullfile(folderOrientationMaps, ...
                [db(k).subject '_' db(k).date '_' num2str(db(k).exp)], ...
                sprintf('%s_plane%02d', titles{iFig}, db(k).planes(iPlane))), ...
                'compact');
            saveas(figHandles(iFig), fullfile(folderOrientationMaps, ...
                [db(k).subject '_' db(k).date '_' num2str(db(k).exp)], ...
                sprintf('%s_plane%02d.jpg',titles{iFig}, ...
                db(k).planes(iPlane))));
            close(figHandles(iFig))
        end
        
        mapSizePixels = [length(info.validY), length(info.validX)];
        ind = fieldOfView(:,1) == info.zoomFactor;
        mapSizeMM = mapSizePixels .* fieldOfView(ind,2:3) ./ ...
            size(info.targetFrame);
        po = prefOris;
        po(rSquared((1:length(neurons))+numCells) < minRSquared) = NaN;
        oriDiffs = gratings.getOrientationDiffs(po,'ori');
        pd = prefDirs;
        pd(rSquared((1:length(neurons))+numCells) < minRSquared) = NaN;
        dirDiffs = gratings.getOrientationDiffs(pd,'dir');
        dist = spatial.getNeuronDistances(info.ROI.CellMaps(neurons), mapSizePixels, ...
            mapSizeMM);
        numNeurons = size(dist,1);
        ind = ones(numNeurons);
        ind = full(spdiags(ind, 1:numNeurons, numNeurons, numNeurons));
        orientationDifferences = [orientationDifferences; oriDiffs(ind==1)];
        directionDifferences = [directionDifferences; dirDiffs(ind==1)];
%         orientationBothDifferences = [orientationBothDifferences; ...
%             combinedDiffs(ind==1)];
        distances = [distances; dist(ind==1)];
        
        numCells = numCells + length(neurons);
    end
    figure
    hist(prefOrisDataset(rSquared>=minRSquared), 0:12.25:180)
    xlabel('Preferred orientation')
    ylabel('# Neurons')
    savefig(gcf, fullfile(folderOrientationMaps, ...
        [db(k).subject '_' db(k).date '_' num2str(db(k).exp)], ...
        'preferredOrientations'), 'compact');
    saveas(gcf, fullfile(folderOrientationMaps, ...
        [db(k).subject '_' db(k).date '_' num2str(db(k).exp)], ...
        'preferredOrientations.jpg'));
%     close(gcf)
    figure
    hist(prefDirsDataset(rSquared>=minRSquared), 0:22.5:360)
    xlabel('Preferred direction')
    ylabel('# Neurons')
    savefig(gcf, fullfile(folderOrientationMaps, ...
        [db(k).subject '_' db(k).date '_' num2str(db(k).exp)], ...
        'preferredDirections'),'compact');
    saveas(gcf, fullfile(folderOrientationMaps, ...
        [db(k).subject '_' db(k).date '_' num2str(db(k).exp)], ...
        'preferredDirections.jpg'));
%     close(gcf)
%     figure
%     hist(prefOrisBothDataset(prefOrisBothDataset(:,2)>=.3,1), 0:12.25:180)
%     xlabel('Preferred orientation')
%     ylabel('# Neurons')
%     savefig(gcf, fullfile(folderOrientationMaps, ...
%         [db(k).subject '_' db(k).date '_' num2str(db(k).exp)], ...
%         'preferredOrientationsCombi'), 'compact');
%     saveas(gcf, fullfile(folderOrientationMaps, ...
%         [db(k).subject '_' db(k).date '_' num2str(db(k).exp)], ...
%         'preferredOrientationsCombi.jpg'));
%     close(gcf)
        
    preferredOrientations = [preferredOrientations; prefOrisDataset];
    preferredDirections = [preferredDirections; prefDirsDataset];
%     preferredOrientationsBoth = [preferredOrientationsBoth; prefOrisBoth];
end

gratings.plotOriDiffsVsDistance(orientationDifferences, distances, 'orientation')
gratings.plotOriDiffsVsDistance(directionDifferences, distances, 'direction')
% gratings.plotOriDiffsVsDistance(orientationBothDifferences, distances, 'orientation')

figure
counts = hist(preferredOrientations(respNeurons & tunedNeurons), 0:11.25:180);
bar(0:11.25:180, counts, 'k')
set(gca, 'box', 'off','XTick',0:45:180)
xlim([-12.25 192.25])
xlabel('Preferred orientation')
ylabel('# Neurons')
text(0,max(counts)*1.1,sprintf('n=%d',sum(respNeurons&tunedNeurons)))
figure
counts = hist(preferredDirections(respNeurons & tunedNeurons), 0:22.5:360);
bar(0:22.5:360, counts, 'k')
set(gca, 'box', 'off','XTick',0:90:360)
xlim([-22.5 382.5])
xlabel('Preferred direction')
ylabel('# Neurons')
text(0,max(counts)*1.1,sprintf('n=%d',sum(respNeurons&tunedNeurons)))
% figure
% hist(preferredOrientationsBoth(preferredOrientationsBoth(:,2)>=.3,1), 0:12.25:180)
% h = findobj(gca, 'Type', 'patch');
% h.FaceColor = 'k';
% h.EdgeColor = 'none';
% set(gca, 'box', 'off')
% xlabel('Preferred orientation')
% ylabel('# Neurons')

figure
bar([sum(~respNeurons), sum(respNeurons&~tunedNeurons), ...
    sum(respNeurons&tunedNeurons)],'k')
set(gca,'box','off','XTickLabel',{'unresponsive','not tuned','tuned'})
ylabel('# Neurons')