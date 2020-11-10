%% Parameters

stimuli = {'gratings', 'grayScreen', 'dark'};
% label = 'neurons';
label = 'boutons';

% for evaluation of receptive fields (significance/goodness)
minExplainedVarianceStim = 0.01;
minPVal = 0.05;
maxLambda = 1;

%% Folders
% file location depending on PCs
folderPC = 'C:\STORAGE\OneDrive - University College London';
% folderPC = 'C:\Users\skgtchr\OneDrive - University College London'; % ZEN notebook
% folderPC = 'C:\Storage\OneDrive - University College London'; % Lenovo laptop
% data folders
folderROIData = fullfile(folderPC, 'Lab\DATA\InfoStructs');
if strcmp(label, 'neurons')
    folderResults = fullfile(folderPC, 'Lab\RESULTS\nonvisualEffects\modelGratingResp');
    folderRFResults = fullfile(folderPC, 'Lab\RESULTS\receptiveFields\SC neurons');
    corrections = [];
else
    folderResults = fullfile(folderPC, 'Lab\RESULTS\boutons\nonvisualEffects');
    folderRFResults = fullfile(folderPC, 'Lab\RESULTS\receptiveFields\boutons');
    data = load(fullfile(folderROIData, 'corrections_boutons.mat'));
    corrections = data.corrections;
    doCorrect = data.doCorrect;
    stimuliCompare = [3 1];
    nonVis = {'running','pupil'};
end

%% Collect correlation values from analyses of continuous traces
data = load(fullfile(folderResults, 'pupil', ...
    'corrsDuringGratingsAndGrayScreen_sigma1.00.mat'));
corrsPupil = data.corrs;
data = load(fullfile(folderResults, 'running', ...
    'corrsDuringGratingsAndGrayScreen_sigma1.00.mat'));
corrsRun = data.corrs;

dataset = [];
plane = [];
neuron = [];
rhosPupil = [];
nullRhosPupil = [];
rhosRun = [];
nullRhosRun = [];
resp = [];
isGad = [];
n = 0;
draws = size(corrsRun(1).plane(1).gratings.nullRhos,2);
for k = 1:length(corrsRun)
    for iPlane = 1:length(corrsRun(k).plane)
        fields = fieldnames(corrsRun(k).plane(iPlane));
        stRun = find(~structfun(@isempty, corrsRun(k).plane(iPlane)),1);
        numCellsRun = length(corrsRun(k).plane(iPlane).(fields{stRun}).rhos);
        stPupil = find(~structfun(@isempty, corrsPupil(k).plane(iPlane)),1);
        numCellsPupil = length(corrsPupil(k).plane(iPlane).(fields{stPupil}).rhos);
        numCells = max(numCellsRun, numCellsPupil);
        dataset = [dataset; ones(numCells,1).*k];
        plane = [plane; ones(numCells,1).*iPlane];
        neuron = [neuron; (1:numCells)'];
        rhosPupil = [rhosPupil; NaN(numCells, length(stimuli))];
        nullRhosPupil = [nullRhosPupil; NaN(numCells, length(stimuli), draws)];
        rhosRun = [rhosRun; NaN(numCells, length(stimuli))];
        nullRhosRun = [nullRhosRun; NaN(numCells, length(stimuli), draws)];
        for st = 1:length(stimuli)
            if isfield(corrsPupil(k).plane(iPlane), stimuli{st}) && ...
                    ~isempty(corrsPupil(k).plane(iPlane).(stimuli{st}).rhos)
                rhosPupil((1:numCells) + n, st) = ...
                    corrsPupil(k).plane(iPlane).(stimuli{st}).rhos;
                nullRhosPupil((1:numCells) + n, st, :) = ...
                    corrsPupil(k).plane(iPlane).(stimuli{st}).nullRhos;
            end
            if isfield(corrsRun(k).plane(iPlane), stimuli{st}) && ...
                    ~isempty(corrsRun(k).plane(iPlane).(stimuli{st}).rhos)
                rhosRun((1:numCells) + n, st) = ...
                    corrsRun(k).plane(iPlane).(stimuli{st}).rhos;
                nullRhosRun((1:numCells) + n, st, :) = ...
                    corrsRun(k).plane(iPlane).(stimuli{st}).nullRhos;
            end
        end
        if isfield(corrsRun(k).plane(iPlane), 'responsive')
            resp = [resp; corrsRun(k).plane(iPlane).responsive];
        elseif isfield(corrsPupil(k).plane(iPlane), 'responsive')
            resp = [resp; corrsPupil(k).plane(iPlane).responsive];
        else resp = [resp; NaN(numCells,1)];
        end
        if isfield(corrsRun(k).plane, 'isGad')
            isGad = [isGad; corrsRun(k).plane(iPlane).isGad];
        elseif isfield(corrsPupil(k).plane, 'isGad')
            isGad = [isGad; corrsPupil(k).plane(iPlane).isGad];
        else isGad = [isGad; NaN(numCells,1)];
        end
        
        n = n + numCells;
    end 
end

%% Collect relevant variables from visual noise data
file = 'receptiveFields.mat';
data = load(fullfile(folderRFResults, file));
RFs = data.RFs;

evStim = NaN(0,1);
pValues = NaN(0,1);
lambdasStim = NaN(0,1);
onOffPeaks = NaN(0,2);
for k = 1:length(corrsRun)
    krf = find(strcmp(corrsRun(k).subject, {RFs.subject}) & ...
        strcmp(corrsRun(k).date, {RFs.date}));
    if isempty(krf) || isempty(RFs(krf).RFTimes) || sum(dataset == k)==0
        numCells = sum(dataset == k);
        evStim(end+(1:numCells),1) = NaN;
        pValues(end+(1:numCells),1) = NaN;
        lambdasStim(end+(1:numCells),1) = NaN;
        onOffPeaks(end+(1:numCells),:) = NaN;
        continue
    end
    for iPlane = 1:length(RFs(krf).planes)
        evStim = [evStim; RFs(krf).plane(iPlane).explainedVariances_stimOnly];
        pValues = [pValues; RFs(krf).plane(iPlane).pVal_RFonly];
        lambdasStim = [lambdasStim; RFs(krf).plane(iPlane).lambdasStim'];
        
        rfs = RFs(krf).plane(iPlane).receptiveFields;
        peaks = NaN(size(rfs,5),2);
        [~,t] = max(max(reshape(permute(abs(rfs),[1 2 4 3 5]), [], ...
            size(rfs,3), size(rfs,5)), [], 1), [], 2);
        t = squeeze(t);
        for iCell = 1:size(rfs,5)
            rf = squeeze(rfs(:,:,t(iCell),:,iCell));
            if all(isnan(rf(:)))
                continue
            end
            [mx,row] = max(max(abs(rf),[],3),[],1);
            [~,col] = max(mx);
            row = row(col);
            rf = squeeze(rf(row,col,:));
            rf(2) = -rf(2); % set values for drive by black squares to positive values
            peaks(iCell,:) = rf;
        end
        onOffPeaks = [onOffPeaks; peaks];
    end
end

validRF = pValues < minPVal & evStim > minExplainedVarianceStim & ...
    lambdasStim < maxLambda;

% get ON-OFF-ratio
[~,type] = max(abs(onOffPeaks),[],2); % determine whether On or Off is stronger
inds = sub2ind(size(onOffPeaks), (1:size(onOffPeaks,1))', type);
signs = sign(onOffPeaks(inds)); % determine whether neuron is driven or suppressed
OnOffRatios = onOffPeaks;
OnOffRatios(signs<0,:) = -OnOffRatios(signs<0,:); % if suppressed, reverse values
OnOffRatios(OnOffRatios<0) = 0; % set values that have opposite drive to dominant RF type to zero
OnOffRatios = (OnOffRatios(:,1)-OnOffRatios(:,2)) ./ sum(OnOffRatios,2);

%% Collect tuning data
data = load(fullfile(folderResults, 'pupil', ...
    'tuning_prefDirSigmaDIFixed_lineFit.mat'));
tuning = data.tuning;
degrees = data.x;
data = load(fullfile(folderResults, 'pupil', ...
    'nullTuning_prefDirSigmaDIFixed.mat'));
null = data.null;

prefDirs = [];
responses = cell(1, length(corrsRun));
isSuppr = [];

fprintf('Dataset (of %d): ', length(corrsRun))
for k = 1:length(corrsRun)
    fprintf('%d ', k)
    ktun = find(strcmp(corrsRun(k).subject, {tuning.subject}) & ...
        strcmp(corrsRun(k).date, {tuning.date}));
    if isempty(ktun) || isempty(tuning(ktun).plane)
        numCells = sum(dataset==k);
        prefDirs = [prefDirs; NaN(numCells,1)];
        responses{k} = NaN(numCells,12,1,2);
        isSuppr = [isSuppr; NaN(numCells,1)];
        continue
    end
    numResp = 1;
    for iPlane = 1:length(tuning(ktun).plane)
        if isempty(tuning(ktun).plane(iPlane).cellIDs)
            continue
        end
        n = find(~cellfun(@isempty, {tuning(ktun).plane(iPlane).cond(1).cell.responses}), 1);
        if isempty(n)
            continue
        end
        numResp = size(tuning(ktun).plane(iPlane).cond(1).cell(n).responses,2);
        break
    end
    for iPlane = 1:length(tuning(ktun).plane)
        numCells = sum(dataset==k & plane==iPlane);
        data = tuning(ktun).plane(iPlane);
        pd = NaN(numCells,1);
        resp = NaN(numCells, 12, numResp, 2);
        spp = NaN(numCells,1);
        if isempty(data.cellIDs)
            prefDirs = [prefDirs; pd];
            responses{k} = [responses{k}; resp];
            isSuppr = [isSuppr; spp];
            continue
        end
        
        spp(data.cellIDs) = data.isSuppressed;
        neurons = find(~cellfun(@isempty, {data.cond(1).cell.parameters}));
        for j = 1:length(neurons)
            iCell = neurons(j);
            ID = data.cellIDs(iCell);
            for c = 1:2
                resp(ID,:,:,c) = data.cond(c).cell(iCell).responses;
                pars = data.cond(c).cell(iCell).parameters;
                if length(pars) > 1 && c == 1
                    pd(ID) = pars(1);
                end
            end
            
            if ~isempty(corrections)
                a = corrections(ktun).plane(iPlane).a{tuning(ktun).exp}(ID);
                b = corrections(ktun).plane(iPlane).b{tuning(ktun).exp}(ID);
                resp(ID,:,:,:) = doCorrect(a,b,resp(ID,:,:,:));
            end
        end
        prefDirs = [prefDirs; pd];
        responses{k} = [responses{k}; resp];
        isSuppr = [isSuppr; spp];
    end
end
fprintf('\n')

%% Determine DSIs and OSIs
cond = 1; % small or large pupil
directions = 0:30:330;
shuffles = 1000;
dirVectors = exp(directions./180.*pi .* 1i);
oriVectors = exp(directions./180.*2.*pi .* 1i);

meanResp = [];
shuffleResp = [];
for k = 1:length(responses)
    r = responses{k}(:,:,:,cond);
    meanResp = [meanResp; nanmean(r, 3)];
    r1 = reshape(r, size(r,1), []);
    shR = NaN(size(r,1),length(dirVectors),shuffles);
    for sh = 1:shuffles
        p = randperm(size(r1,2));
        r2 = reshape(r1(:,p), size(r));
        shR(:,:,sh) = nanmean(r2, 3);
    end
    shuffleResp = [shuffleResp; shR];
end

meanResp(isSuppr==1,:) = -meanResp(isSuppr==1,:);
shuffleResp(isSuppr==1,:,:) = -shuffleResp(isSuppr==1,:,:);
meanResp(meanResp<0) = 0;
shuffleResp(shuffleResp<0) = 0;
meanResp = meanResp ./ max(meanResp,[],2);
shuffleResp = shuffleResp ./ max(shuffleResp,[],2);

% Determine DSIs
vects = mean(dirVectors .* meanResp, 2);
shuffleVects = squeeze(mean(dirVectors .* shuffleResp, 2));
DSIs = abs(vects);
nullDSIs = abs(shuffleVects);
p_DSI = sum(nullDSIs > DSIs,2) ./ shuffles;
p_DSI(isnan(DSIs)) = NaN;

% Determine OSIs
vects = mean(oriVectors .* meanResp, 2);
shuffleVects = squeeze(mean(oriVectors .* shuffleResp, 2));
OSIs = abs(vects);
nullOSIs = abs(shuffleVects);
p_OSI = sum(nullOSIs > OSIs,2) ./ shuffles;
p_OSI(isnan(OSIs)) = NaN;

%% Plot correlations against each other
onThr = 0.6;
offThr = -0.6;
axLimits = [-.65 .85];
binSize = 0.2;
edges = -0.6:binSize:0.8;
bins = edges(1:end-1)+binSize/2;

groupNames = {{'all'},{'ON', 'ON/OFF', 'OFF'},{'enhanced','suppressed'}, ...
    {'DS','OS','not selective'}};
groupInds = {true(size(rhosRun,1),1), ...
    [OnOffRatios>onThr, OnOffRatios>=offThr&OnOffRatios<=onThr, OnOffRatios<offThr], ...
    [isSuppr==-1, isSuppr==1], ...
    [p_DSI<.05, p_OSI<.05, p_DSI>=.05&p_OSI>=.05]};
groupColors = {[0 0 0], [1 0 0;.5 0 .5;0 0 1], lines(2), [1 0 0;.5 0 .5;0 0 1]};

rhos = cell(1,2);
for c = 1:2
    if strcmp(nonVis{c}, 'running')
        rhos{c} = rhosRun(:,stimuliCompare(c));
    else
        rhos{c} = rhosPupil(:,stimuliCompare(c));
    end
end
valid = all(~isnan(cat(2, rhos{:})),2);

for g = 1:length(groupNames)
    figure
    hold
    plot(axLimits, [0 0], 'k')
    plot([0 0], axLimits, 'k')
    h = zeros(1, size(groupInds{g},2));
    for type = 1:size(groupInds{g},2)
        r1 = rhos{1}(groupInds{g}(:,type));
        r2 = rhos{2}(groupInds{g}(:,type));
        h(type) = scatter(r1, r2, ...
            [], 'filled', 'MarkerFaceColor', groupColors{g}(type,:), ...
            'MarkerFaceAlpha', 0.1);
    end
    for type = 1:size(groupInds{g},2)
        r1 = rhos{1}(groupInds{g}(:,type));
        r2 = rhos{2}(groupInds{g}(:,type));
        [N,~,b] = histcounts(r1, edges);
        yMeans = zeros(1, length(bins));
        ySEMs = zeros(1, length(bins));
        for bin = 1:length(bins)
            vals = r2(b==bin);
            yMeans(bin) = nanmean(vals);
            ySEMs(bin) = nanstd(vals)./sqrt(sum(~isnan(vals)));
        end
        errorbar(bins, yMeans, ySEMs, 'o', 'MarkerSize', 8, ...
            'Color', groupColors{g}(type,:), ...
            'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', groupColors{g}(type,:), 'LineWidth', 2, 'CapSize', 0)
    end
    legend(h, groupNames{g}, 'Location', 'NorthWest')
    axis([axLimits axLimits])
    axis square
    set(gca, 'box', 'off')
    xlabel(sprintf('Correlation with %s during %s', nonVis{1}, stimuli{stimuliCompare(1)}))
    ylabel(sprintf('Correlation with %s during %s', nonVis{2}, stimuli{stimuliCompare(2)}))
    title(sprintf('%s (n=%d)', label, sum(any(groupInds{g},2) & valid)))
end