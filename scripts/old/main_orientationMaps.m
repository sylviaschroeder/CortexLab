%% Load database
db_driftingGratings

%% Folders
folderROIData = 'C:\Temp2p';
folderOrientationMaps = 'C:\Temp2p\Results\orientation\maps';
folderTuningModelRunning = ['C:\Temp2p\Results\nonVisualEffects\tuningModelledRunning_' stimType];

%% Parameters
fieldOfView = [2 500 520; 2.2 450 490; 2.5 400 430; 3 335 370];
% to evaluate fit of tuning curve
minRSquared = 0.2;

%% Loop across datasets
preferredOrientations = {};
preferredDirections = {};
distances = {};

respNeurons = {};       % 1 if ANOVA (on responses: [stimulus x trial]) was significant
tunedNeurons = {};    % 1 if rSquared of fit to tuning curve was >= 0.2
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
    pars = results.tuningPars(:,1,1);

    prefOrisDataset = [];
    prefDirsDataset = [];
    rSquaredDataset = [];
    for iPlane=1:length(db(k).planes)
        % get results only from neurons of this plane
        indPl = find(results.corr.orderedNeurons(:,1) == db(k).planes(iPlane));
        [~,order] = sort(results.corr.orderedNeurons(indPl,2));
        indPl = indPl(order);
        
        % load meta
        data=load(fullfile(folder, sprintf(file,db(k).planes(iPlane))));
        info = data.meta;
        neurons = find(strcmp(data.meta.ROI.CellClasses, 's'));
    
        %% Get stimulus pref., maps, pref. differences, and distances
        prefOris = mod(pars(indPl),180);
        prefDirs = mod(pars(indPl),360);
        prefOrisDataset = [prefOrisDataset; prefOris];
        prefDirsDataset = [prefDirsDataset; prefDirs];
        preferredOrientations{end+1} = prefOris;
        preferredDirections{end+1} = prefDirs;
        respNeurons{end+1} = respondingNeurons(indPl);
        rSquared = results.errors(indPl,1);
        rSquaredDataset = [rSquaredDataset; rSquared];
        tunedNeurons{end+1} = rSquared >= minRSquared;
        
        figHandles(1) = gratings.plotOrientationMap(...
            [prefOris, rSquared], 'ori', info.ROI.CellMaps(neurons), ...
            [length(info.validY) length(info.validX)], minRSquared);
        figHandles(2) = gratings.plotOrientationMap(...
            [prefDirs, rSquared], 'dir', info.ROI.CellMaps(neurons), ...
            [length(info.validY) length(info.validX)], minRSquared);
        titles = {'orientationMap','directionMap'};
%         for iFig = 1:2
%             savefig(figHandles(iFig), fullfile(folderOrientationMaps, ...
%                 [db(k).subject '_' db(k).date '_' num2str(db(k).exp)], ...
%                 sprintf('%s_plane%02d', titles{iFig}, db(k).planes(iPlane))), ...
%                 'compact');
%             saveas(figHandles(iFig), fullfile(folderOrientationMaps, ...
%                 [db(k).subject '_' db(k).date '_' num2str(db(k).exp)], ...
%                 sprintf('%s_plane%02d.jpg',titles{iFig}, ...
%                 db(k).planes(iPlane))));
%             close(figHandles(iFig))
%         end
        
        mapSizePixels = [length(info.validY), length(info.validX)];
        ind = fieldOfView(:,1) == info.zoomFactor;
        mapSizeMM = mapSizePixels .* fieldOfView(ind,2:3) ./ ...
            size(info.targetFrame);
        dist = spatial.getNeuronDistances(info.ROI.CellMaps(neurons), mapSizePixels, ...
            mapSizeMM);
        distances{end+1} = dist;
    end
    figure
    hist(prefOrisDataset(rSquared>=minRSquared), 0:12.25:180)
    xlabel('Preferred orientation')
    ylabel('# Neurons')
%     savefig(gcf, fullfile(folderOrientationMaps, ...
%         [db(k).subject '_' db(k).date '_' num2str(db(k).exp)], ...
%         'preferredOrientations'), 'compact');
%     saveas(gcf, fullfile(folderOrientationMaps, ...
%         [db(k).subject '_' db(k).date '_' num2str(db(k).exp)], ...
%         'preferredOrientations.jpg'));
%     close(gcf)
    figure
    hist(prefDirsDataset(rSquared>=minRSquared), 0:22.5:360)
    xlabel('Preferred direction')
    ylabel('# Neurons')
%     savefig(gcf, fullfile(folderOrientationMaps, ...
%         [db(k).subject '_' db(k).date '_' num2str(db(k).exp)], ...
%         'preferredDirections'),'compact');
%     saveas(gcf, fullfile(folderOrientationMaps, ...
%         [db(k).subject '_' db(k).date '_' num2str(db(k).exp)], ...
%         'preferredDirections.jpg'));
%     close(gcf)
end

% Plot measured distances versus differences in pref. orientation/direction
orientationDifferences = cell(size(preferredOrientations));
directionDifferences = cell(size(preferredOrientations));
dist = cell(size(preferredOrientations));
permOriDiffs = cell(size(preferredOrientations));
permDirDiffs = cell(size(preferredOrientations));
numPerm = 1000;
for k = 1:length(preferredOrientations)
    prefOris = preferredOrientations{k};
    prefOris(respNeurons{k}==0 & tunedNeurons{k}==0) = [];
    oriDiffs = gratings.getOrientationDiffs(prefOris,'ori');
    prefDirs = preferredDirections{k};
    prefDirs(respNeurons{k}==0 & tunedNeurons{k}==0) = [];
    dirDiffs = gratings.getOrientationDiffs(prefDirs,'dir');
    
    numNeurons = size(oriDiffs,1);
    ind = ones(numNeurons);
    ind = full(spdiags(ind, 1:numNeurons, numNeurons, numNeurons));
    orientationDifferences{k} = oriDiffs(ind==1);
    directionDifferences{k} = dirDiffs(ind==1);
    d = distances{k};
    d(respNeurons{k}==0 & tunedNeurons{k}==0,:) = [];
    d(:,respNeurons{k}==0 & tunedNeurons{k}==0) = [];
    dist{k} = d(ind==1);
    
    permODs = NaN(length(orientationDifferences{k}),numPerm);
    permDDs = NaN(length(orientationDifferences{k}),numPerm);
    for p = 1:numPerm
        s = RandStream('mt19937ar','Seed',p);
        order = randperm(s,numNeurons);
        po = prefOris(order);
        pod = gratings.getOrientationDiffs(po,'ori');
        permODs(:,p) = pod(ind==1);
        pd = prefDirs(order);
        pdd = gratings.getOrientationDiffs(pd,'dir');
        permDDs(:,p) = pdd(ind==1);
    end
    permOriDiffs{k} = permODs;
    permDirDiffs{k} = permDDs;
    gratings.plotOriDiffsVsDistance(orientationDifferences{k}, dist{k}, ...
        'orientation', permODs);
    gratings.plotOriDiffsVsDistance(directionDifferences{k}, dist{k}, ...
        'direction', permDDs);
end

oriDiffs = cat(1, orientationDifferences{:});
permODs = cat(1, permOriDiffs{:});
dirDiffs = cat(1, directionDifferences{:});
permDDs = cat(1, permDirDiffs{:});
allDist = cat(1, dist{:});

gratings.plotOriDiffsVsDistance(oriDiffs, allDist, 'orientation', permODs)
gratings.plotOriDiffsVsDistance(dirDiffs, allDist, 'direction', permDDs)

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