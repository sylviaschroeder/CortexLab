%% Folders:
folders.data = 'C:\Users\Sylvia\OneDrive - University of Sussex\Projects\2023_OrientationColumns\DataToPublish';
folders.tools = 'C:\dev\toolboxes';
folders.repo = 'C:\dev\workspaces\he_schroeder_columns';
folderScript = 'C:\dev\workspaces\CortexLab\scripts\convert2NPY';

%% Add paths
addpath(genpath(fullfile(folders.tools, 'npy-matlab')))
addpath(genpath(fullfile(folders.repo)))
addpath(genpath(fullfile(folderScript)));

%% Parameters
minCorr = 0.4; % neurons: 0.4;
maxDistXY = 5; % neurons: 5;
filtWindow = 5; % in samples
stimTypes = {'gratingsDrifting', 'gratingsStatic', 'bars'};

%% Load database
k = 1;
db(k).subject = 'SS041'; % Vane
db(k).date = '2015-04-11';
db(k).planes = 1:4; % set to be 8 um apart
k=k+1;
db(k).subject = 'SS044'; % Vane
db(k).date = '2015-04-30';
db(k).planes = 1:3; % set to be 10 um apart
k=k+1;
db(k).subject = 'SS044'; % Vane
db(k).date = '2015-05-15';
db(k).planes = 1:3; % set to be 10 um apart

%% Convert data
for k = 2:length(db)
    f = fullfile(folders.data, 'neurons', db(k).subject, db(k).date);

    % load data
    calc = io.getCalciumData(f);
    recData = io.getRecordingInfo(f);
    pos = recData.roiPositions;

    % pairwise distances in XY
    dist = sqrt( (pos(:,1)-pos(:,1)').^2 + (pos(:,2)-pos(:,2)').^2 );
    % disregard pairs in the same plane or farther apart than 2 planes
    planeDist = abs(calc.planes - calc.planes');
    dist(planeDist == 0) = NaN;
    dist(planeDist > 2) = NaN;

    % disregard neurons too far apart in XY
    dist(dist > maxDistXY) = NaN;

    % find cells that have close neighbour
    hasNeighbour = find(any(~isnan(dist),2));

    % reduce distance matrix to those neurons with neighbours
    dist = dist(hasNeighbour, hasNeighbour);

    % determine correlation between traces
    corrMat = corr(medfilt1(calc.traces(:,hasNeighbour), ...
        filtWindow, [], 1, 'omitnan'), 'rows', 'pairwise');

    % disregard pairs with small correlations
    dist(corrMat < minCorr) = NaN;
    ind = find(~all(isnan(dist),1));
    dist = dist(ind,ind);
    hasNeighbour = hasNeighbour(ind);

    duplicates = {};
    dist2 = dist;
    for j = 1:size(dist,2)
        % find duplicates for item j
        d = find(~isnan(dist2(:,j)));
        if isempty(d)
            continue
        end
        % put them in a queue
        queue =  d;
        % add j itself to the list of duplicates
        d = [j; d];
        while ~isempty(queue)
            % for each member in the queue check whether it has further
            % duplicates and add those to the queue
            d_add = find(~isnan(dist2(:,queue(1))));
            for i = d_add'
                if ~ismember(i, d)
                    d = [d; i];
                    queue(end+1,1) = i;
                end
            end
            % do not check found duplicates again
            dist2(:,queue(1)) = NaN;
            queue(1) = [];
        end
        duplicates{end + 1} = d;
    end

    % go through list of duplicates, choose the best representative and
    % disregard others
    isUnique = true(size(calc.traces,2),1);
    for j = 1:length(duplicates)
        neighbours = hasNeighbour(duplicates{j});
        SNR = diff( prctile( medfilt1(calc.traces(:, neighbours), ...
            filtWindow, [], 1, "omitnan"), [50 98], 1), 1, 1) ./ ...
            mad(calc.traces(:, neighbours), 1, 1);
        dupl = neighbours(SNR < max(SNR));
        isUnique(dupl) = false;
    end
    
    % update data --> delete duplicate neurons
    % 2pCalcium, 2pRois
    data = io.getCalciumData(f);
    writeNPY(data.traces(:, isUnique), fullfile(f, '_ss_2pCalcium.dff.npy'))
    writeNPY(data.planes(isUnique), fullfile(f, '_ss_2pRois._ss_2pPlanes.npy'))
    writeNPY(data.ids(isUnique), fullfile(f, '_ss_2pRois.ids.npy'))
    data = io.getRecordingInfo(f);
    writeNPY(data.roiMasks(isUnique,:), fullfile(f, '_ss_2pRois.masks.npy'))
    writeNPY(data.roiPositions(isUnique,:), fullfile(f, '_ss_2pRois.xyz.npy'))
    writeNPY(data.isInhibitory(isUnique), fullfile(f, '_ss_2pRois.isGad.npy'))
    % _ss_bars..., _ss_gratingsDrifting...,  _ss_gratingsStatic...
    for st = 1:length(stimTypes)
        type = stimTypes{st};
        % ignore session if stimulus was not presented
        if ~isfile(fullfile(f, sprintf('_ss_%s.intervals.npy', type)))
            continue
        end
        data = io.getStimResponseFits(f, type);
        data.kernel(:,~isUnique) = [];
        data.amplitudes(:,:,~isUnique) = [];
        data.prediction(:,~isUnique) = [];
        data.pValue(~isUnique) = [];
        data.R2(~isUnique) = [];
        if st == 3 % bars
            data.lags(:,:,~isUnique) = [];
        end
        time.kernel = data.time_kernel;
        time.prediction = data.time_prediction;
        io.writeKernelFitResults(data, time, f, type);

        [dirTuning, oriTuning] = io.getTuningResults(f, type);
        oriTuning.preference(~isUnique) = [];
        oriTuning.selectivity(~isUnique) = [];
        oriTuning.pValue(~isUnique) = [];
        oriTuning.responseSign(~isUnique) = [];
        if st ~= 2
            dirTuning.preference(~isUnique) = [];
            dirTuning.selectivity(~isUnique) = [];
            dirTuning.pValue(~isUnique) = [];
            dirTuning.responseSign(~isUnique) = [];
        end
        io.writeTuningResults(dirTuning, oriTuning, f, type);
    end
    % _ss_rf...
    data=io.getRFFits(f);
    writeNPY(data.maps(isUnique,:,:,:,:), fullfile(f, '_ss_rf.maps.npy'));
    writeNPY(data.bestSubFields(isUnique), fullfile(f, '_ss_rf.bestSubField.npy'));
    writeNPY(data.subFieldSigns(isUnique,:), fullfile(f, '_ss_rf.subFieldSigns.npy'));
    writeNPY(data.fitParameters(isUnique,:), fullfile(f, '_ss_rf.gaussFitPars.npy'));
    writeNPY(data.peaks(isUnique), fullfile(f, '_ss_rf.peak.npy'));
    writeNPY(data.gaussMasks(isUnique,:,:,:), fullfile(f, '_ss_rf.gaussMask.npy'));
    writeNPY(data.timeWeights(isUnique,:), fullfile(f, '_ss_rf.gaussTimeWeights.npy'));
    writeNPY(data.EV(isUnique), fullfile(f, '_ss_rf.explVars.npy'));
    writeNPY(data.predictions(:,isUnique), fullfile(f, '_ss_rfPrediction.traces.npy'));

    outliers = readNPY(fullfile(f, '_ss_rf.outliers.npy'));
    writeNPY(outliers(isUnique), fullfile(f, '_ss_rf.outliers.npy'));
    pred = readNPY(fullfile(f, '_ss_rf.posRetinotopy.npy'));
    writeNPY(pred(isUnique,:), fullfile(f, '_ss_rf.posRetinotopy.npy'));
end