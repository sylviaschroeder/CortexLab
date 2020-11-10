%% Load database
db_blanks
% db_driftingGratings

%% Folders
folderROIData = 'C:\DATA\InfoStructs';
folderCorrelationMaps = 'C:\RESULTS\nonVisualEffects\maps';
folderTuningModelRunning = ['C:\RESULTS\nonVisualEffects\tuningModelledRunning_' stimType];

%% Parameters
fieldOfView = [2 500 520; 2.2 450 490; 2.5 400 430; 3 335 370];

%% Loop across datasets
corrValues = {};
distances = {};

for k = 1:length(db)
    if ~isdir(fullfile(folderCorrelationMaps, [db(k).subject '_' ...
            db(k).date '_' num2str(db(k).exp)]))
        mkdir(fullfile(folderCorrelationMaps, [db(k).subject '_' ...
            db(k).date '_' num2str(db(k).exp)]));
    end
    
    % load correlation results
    data = load(fullfile(folderTuningModelRunning, [db(k).subject '_' ...
        db(k).date '_' num2str(db(k).exp)], 'results.mat'), 'results');
    results = data.results;
    folder = fullfile(folderROIData, db(k).subject, db(k).date, ...
        num2str(db(k).exp));
    file = [db(k).date '_' num2str(db(k).exp) '_' ...
        db(k).subject '_2P_plane%03d_correctedF.mat'];
    allRhos = results.corr.rho;
    allNeurons = results.corr.orderedNeurons;
    for iPlane = 1:length(db(k).planes)
        % load ROIs
        data = load(fullfile(folder, sprintf(file, db(k).planes(iPlane))));
        ROImaps = data.meta.ROI.CellMaps;
        % find ROIs that have rho value and that are neurons
        neuronInds = find(allNeurons(:,1)==db(k).planes(iPlane));
        [sortedInds, order] = sort(allNeurons(neuronInds,2));
        ROImaps = ROImaps(sortedInds);
        rhos = allRhos(neuronInds);
        rhos = rhos(order);
        invalidInds = ~strcmp('s', data.meta.ROI.CellClasses(sortedInds));
        rhos(invalidInds) = [];
        ROImaps(invalidInds) = [];
        
        corrValues{end+1} = rhos;
        
        handle = nonVis.plotCorrelationMap(rhos,ROImaps,data.arrSize([1 2]));
        title(handle, sprintf('%s %s %d plane %d',db(k).subject,db(k).date, ...
            db(k).exp,db(k).planes(iPlane)), 'Interpreter', 'none')
        savefig(gcf, fullfile(folderCorrelationMaps, ...
            [db(k).subject '_' db(k).date '_' num2str(db(k).exp)], ...
            sprintf('corrMap_plane%02d', db(k).planes(iPlane))), ...
            'compact');
        saveas(gcf, fullfile(folderCorrelationMaps, ...
            [db(k).subject '_' db(k).date '_' num2str(db(k).exp)], ...
            sprintf('corrMap_plane%02d.jpg', db(k).planes(iPlane))));
        close(gcf)
        
        mapSizePixels = [length(data.meta.validY), length(data.meta.validX)];
        ind = fieldOfView(:,1) == data.meta.zoomFactor;
        mapSizeMM = mapSizePixels .* fieldOfView(ind,2:3) ./ ...
            size(data.meta.targetFrame);
        dist = spatial.getNeuronDistances(ROImaps, mapSizePixels, ...
            mapSizeMM);
        distances{end+1} = dist;
    end
end

corrDifferences = cell(size(corrValues));
dist = cell(size(corrValues));
permCorrDiffs = cell(size(corrValues));
numPerm = 1000;
for k = 1:length(corrValues)
    cVals = corrValues{k};
    cDiffs = abs(bsxfun(@minus, cVals, cVals'));
    
    numNeurons = size(cDiffs,1);
    ind = ones(numNeurons);
    ind = full(spdiags(ind, 1:numNeurons, numNeurons, numNeurons));
    corrDifferences{k} = cDiffs(ind==1);
    dist{k} = distances{k}(ind==1);
    
    permCorrs = NaN(length(corrDifferences{k}),numPerm);
    for p = 1:numPerm
        s = RandStream('mt19937ar','Seed',p);
        order = randperm(s,numNeurons);
        cv = cVals(order);
        cvd = abs(bsxfun(@minus, cv, cv'));
        permCorrs(:,p) = cvd(ind==1);
    end
    permCorrDiffs{k} = permCorrs;
    
    nonVis.plotCorrDiffsVsDistance(corrDifferences{k}, dist{k}, ...
        permCorrs)
end

cDiffs = cat(1, corrDifferences{:});
permCorrs = cat(1, permCorrDiffs{:});
allDist = cat(1, dist{:});

nonVis.plotCorrDiffsVsDistance(cDiffs, allDist, permCorrs);