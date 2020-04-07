%% Define data

% label = 'neurons';
label = 'boutons';

if strcmp(label, 'neurons')
    exampleSets(1).subject = 'M150410_SS044'; % corr gratings: -0.2212, p = 0.024
    exampleSets(1).date = '2015-05-29';       % no gray screen
    exampleSets(1).plane = 2;                 % DI pref resp: -0.0776, p = 0.02
    exampleSets(1).neuron = 72;
    exampleSets(1).RFtype = 'ON';
    exampleSets(1).corrDirection = -1; 
    exampleSets(2).subject = 'M150410_SS044'; % corr gratings: 0.2021, p = 0.036
    exampleSets(2).date = '2015-04-28';       % corr gray screen: -0.0144, p = 0.436
    exampleSets(2).plane = 2;                 % DI pref resp: 0.1324, p = 0.01
    exampleSets(2).neuron = 204;
    exampleSets(2).RFtype = 'ON';
    exampleSets(2).corrDirection = 1;
    exampleSets(3).subject = 'M150305_SS041'; % corr gratings: -0.2171, p = 0.034
    exampleSets(3).date = '2015-04-23';       % corr gray screen: -0.2397, p = 0.002
    exampleSets(3).plane = 4;                 % DI pref resp: -0.1449, p = 0.01
    exampleSets(3).neuron = 58;
    exampleSets(3).RFtype = 'OFF';
    exampleSets(3).corrDirection = -1;
    exampleSets(4).subject = 'M150611_SS048'; % corr gratings: 0.3802, p = 0.002
    exampleSets(4).date = '2015-12-02';       % corr gray screen: 0.5246, p = 0.01
    exampleSets(4).plane = 3;                 % DI pref resp: 0.1956, p = 0.02
    exampleSets(4).neuron = 62;               % alternatives: set 10, pl 3, neuron 15;
    exampleSets(4).RFtype = 'OFF';            % set 3, pl 3, neuron 159
    exampleSets(4).corrDirection = 1;         % set 3, pl 4, neuron 60
    exampleSets(5).subject = 'M150305_SS041'; % corr gratings: -0.3017, p = 0.0000
    exampleSets(5).date = '2015-04-23';       % corr gray screen: -0.5497, p = 0.000
    exampleSets(5).plane = 4;                 % DI pref resp: -0.2956, p = 0.0000
    exampleSets(5).neuron = 43;
    exampleSets(5).RFtype = 'ON+OFF';
    exampleSets(5).corrDirection = -1;
    exampleSets(6).subject = 'M150410_SS044';
    exampleSets(6).date = '2015-04-28';
    exampleSets(6).plane = 3;
    exampleSets(6).neuron = 188;
    exampleSets(6).RFtype = 'ON+OFF';
    exampleSets(6).corrDirection = 1;
else % boutons in sSC
    exampleSets(1).subject = 'SS078';
    exampleSets(1).date = '2017-10-05';
    exampleSets(1).plane = 1;
    exampleSets(1).neuron = 156;
    exampleSets(1).RFtype = 'ON';
    exampleSets(1).corrDirection = -1;
    exampleSets(2).subject = 'SS078';
    exampleSets(2).date = '2017-10-05';
    exampleSets(2).plane = 1;
    exampleSets(2).neuron = 225;
    exampleSets(2).RFtype = 'ON';
    exampleSets(2).corrDirection = 1;
    exampleSets(3).subject = 'SS078';
    exampleSets(3).date = '2017-10-05';
    exampleSets(3).plane = 1;
    exampleSets(3).neuron = 22;
    exampleSets(3).RFtype = 'OFF';
    exampleSets(3).corrDirection = -1;
    exampleSets(4).subject = 'SS078';
    exampleSets(4).date = '2017-09-28';
    exampleSets(4).plane = 1;
    exampleSets(4).neuron = 243;
    exampleSets(4).RFtype = 'OFF';
    exampleSets(4).corrDirection = 1;
    exampleSets(5).subject = 'M160923_SS069';
    exampleSets(5).date = '2016-10-21';
    exampleSets(5).plane = 1;
    exampleSets(5).neuron = 223;
    exampleSets(5).RFtype = 'ON+OFF';
    exampleSets(5).corrDirection = -1;
    exampleSets(6).subject = 'SS078';
    exampleSets(6).date = '2017-09-28';
    exampleSets(6).plane = 1;
    exampleSets(6).neuron = 263;
    exampleSets(6).RFtype = 'ON+OFF';
    exampleSets(6).corrDirection = 1;
end

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
end

%% Parameters
RFtypes = {'ON','OFF','ON+OFF'};

% for evaluation of receptive fields (significance/goodness)
minExplainedVarianceStim = 0.01;
minPVal = 0.05;
maxLambda = 1;

% to group into ON, OFF and ON+OFF
onThr = 0.5; % ON field is at least three times stronger than OFF field
offThr = -0.5; % OFF field is at least three times stronger than ON field

stimuli = {'gratings', 'grayScreen', 'dark'};

% for tests across datasets
minSamples = 10;

% colors
red = [1 0 .5];
blue = [0 .5 1];

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

pValues_rhoRun = sum(nullRhosRun > rhosRun,3) ./ size(nullRhosRun,3);
pValues_rhoRun(pValues_rhoRun > 0.5) = 1 - pValues_rhoRun(pValues_rhoRun > 0.5);
pValues_rhoRun(isnan(rhosRun)) = NaN;
pValues_rhoPupil = sum(nullRhosPupil > rhosPupil,3) ./ size(nullRhosPupil,3);
pValues_rhoPupil(pValues_rhoPupil > 0.5) = 1 - pValues_rhoPupil(pValues_rhoPupil > 0.5);
pValues_rhoPupil(isnan(rhosPupil)) = NaN;

%% Collect relevant variables from visual noise data
file = 'receptiveFields.mat';
data = load(fullfile(folderRFResults, file));
RFs = data.RFs;

evStim = NaN(0,1);
pValues_RF = NaN(0,1);
lambdasStim = NaN(0,1);
onOffPeaks = NaN(0,2);
for k = 1:length(corrsRun)
    krf = find(strcmp(corrsRun(k).subject, {RFs.subject}) & ...
        strcmp(corrsRun(k).date, {RFs.date}));
    if isempty(krf) || isempty(RFs(krf).RFTimes) || sum(dataset == k)==0
        numCells = sum(dataset == k);
        evStim(end+(1:numCells),1) = NaN;
        pValues_RF(end+(1:numCells),1) = NaN;
        lambdasStim(end+(1:numCells),1) = NaN;
        onOffPeaks(end+(1:numCells),:) = NaN;
        continue
    end
    for iPlane = 1:length(RFs(krf).planes)
        evStim = [evStim; RFs(krf).plane(iPlane).explainedVariances_stimOnly];
        pValues_RF = [pValues_RF; RFs(krf).plane(iPlane).pVal_RFonly];
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

validRF = pValues_RF < minPVal & evStim > minExplainedVarianceStim & ...
    lambdasStim < maxLambda;

% get ON-OFF-ratio
[~,type] = max(abs(onOffPeaks),[],2); % determine whether On or Off is stronger
inds = sub2ind(size(onOffPeaks), (1:size(onOffPeaks,1))', type);
signs = sign(onOffPeaks(inds)); % determine whether neuron is driven or suppressed
OnOffRatios = onOffPeaks;
OnOffRatios(signs<0,:) = -OnOffRatios(signs<0,:); % if suppressed, reverse values
OnOffRatios(OnOffRatios<0) = 0; % set values that have opposite drive to dominant RF type to zero
OnOffRatios = (OnOffRatios(:,1)-OnOffRatios(:,2)) ./ sum(OnOffRatios,2);

%% Plot cumulative distributions of correlations with pupil/running for cell types
rhos = {rhosRun, rhosPupil};
nullRhos = {nullRhosRun, nullRhosPupil};
nonvisNames = {'running','pupil'};
neuronSets = {{OnOffRatios>onThr & validRF; ...
    OnOffRatios>=offThr & OnOffRatios<=onThr & validRF; ...
    OnOffRatios<offThr & validRF} ... % ON, Absolute, OFF
    };
setNames = {{'ON','ON+OFF','OFF'}};
colors = lines(3);
exColors = [red; (red+blue)./2; blue];
if strcmp(label, 'neurons')
    combinations = [2 2; 2 1]; % [nonvis stimulus]
else
    combinations = [1 3; 2 1]; % [nonvis stimulus]
end

% Plot cumulative distribtutions of correlation values
xLimits = [-0.5 0.6];
indRhos = ~any(isnan([rhos{combinations(1,1)}(:,combinations(1,2)), ...
    rhos{combinations(2,1)}(:,combinations(2,2))]),2);
for combi = 1:size(combinations,1)
    nonvis = combinations(combi,1);
    st = combinations(combi,2);
        
    if all(isnan(rhos{nonvis}(:,st)))
        continue
    end
    for s = 1:length(neuronSets)
        figure
        hold on
        plot([0 0], [0 1.2], 'k')
        h = zeros(1, length(neuronSets{s}));
        lbls = cell(1, length(neuronSets{s}));
        for type = 1:length(neuronSets{s})
            ind = neuronSets{s}{type} & indRhos;
            hs = general.plotCumHist(gca, rhos{nonvis}(ind,st), ...
                squeeze(nullRhos{nonvis}(ind,st,:)), [-1 1], colors(type,:));
            plot(mean(rhos{nonvis}(ind,st)), 1.1, 'v', 'MarkerEdgeColor', 'none', ...
                'MarkerFaceColor', colors(type,:))
            h(type) = hs(1);
            lbls{type} = sprintf('%s (n=%d)', setNames{s}{type}, sum(ind));
            
            for ex = 1:length(exampleSets)
                if ~strcmp(exampleSets(ex).RFtype, setNames{s}{type})
                    continue
                end
                exSet = find(strcmp(exampleSets(ex).subject, {corrsRun.subject}) & ...
                    strcmp(exampleSets(ex).date, {corrsRun.date}));
                exInd = find(dataset==exSet & plane==exampleSets(ex).plane & ...
                    neuron==exampleSets(ex).neuron);
                exX = rhos{nonvis}(exInd,st);
                exY = sum(rhos{nonvis}(ind,st) <= exX) / sum(ind);
                if exampleSets(ex).corrDirection > 0
                    c = exColors(type,:) .* 0.7 + [1 1 1].*0.3;
                else
                    c = exColors(type,:) .* 0.7 + [0 0 0].*0.3;
                end
                plot(exX, exY, 'o', 'MarkerFaceColor', c, ...
                    'MarkerEdgeColor', 'none')
            end
        end
        xlim(xLimits)
        ylim([0 1.2])
        legend(h, lbls, 'Location', 'SouthEast')
        xlabel(sprintf('Correlation with %s', nonvisNames{nonvis}))
        ylabel(sprintf('Proportion of %s', label))
        title(sprintf('Correlation with %s during %s', ...
            nonvisNames{nonvis}, stimuli{st}))
    end
end

%% Test across datasets whether median correlation is different between RF types
% (ON, OFF, ON+OFF)
groups = [validRF & OnOffRatios>onThr, ... % ON
    validRF & OnOffRatios<offThr, ...      % OFF
    validRF & OnOffRatios<=onThr & OnOffRatios>=offThr, ... % ON+OFF
    ~validRF]; % no RF fitted

dSets = unique(dataset);
corrMedians = NaN(length(dSets), size(groups,2));
pPerDataset = NaN(length(dSets), size(groups,2), size(groups,2)); % ON - OFF, ON - ON+OFF, OFF - ON+OFF
for d = 1:length(dSets)
    for g = 1:size(groups,2)
        ind = groups(:,g) & dataset==dSets(d);
        if sum(ind) < minSamples
            continue
        end
        corrMedians(d,g) = nanmedian(rhosPupil(ind,1));
        for g2 = g+1:size(groups,2)
            ind2 = groups(:,g2) & dataset==dSets(d);
            if sum(ind) < minSamples
                continue
            end
            p = ranksum(rhosPupil(ind,1), rhosPupil(ind2,1));
            pPerDataset(d,g,g2) = p;
        end
    end
end
pPerDataset = reshape(pPerDataset, length(dSets), []);
ind = ~any(isnan(corrMedians(:,1:3)),2);
[p,~,stats] = friedman(corrMedians(ind,1:3),1,'off');
fprintf('%d datasets, p = %.4f, medians: ', sum(ind), p)
fprintf('%.4f ', median(corrMedians(ind,1:3)))
fprintf('\n')