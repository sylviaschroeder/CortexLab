%% Define data

% label = 'neurons';
label = 'boutons';

%% Parameter
% for evaluation of receptive fields (significance/goodness)
minExplainedVarianceStim = 0.01;
minPVal = 0.05;
maxLambda = 1;

onThresh = 0.5;
offThresh = -0.5;

% for tests across datasets
minSamples = 10;

%% Folders
% file location depending on PCs
folderPC = 'C:\STORAGE\OneDrive - University College London';
% folderPC = 'C:\Users\skgtchr\OneDrive - University College London'; % ZEN notebook
% folderPC = 'C:\Storage\OneDrive - University College London'; % Lenovo laptop
% data folders
folderROIData = fullfile(folderPC, 'Lab\DATA\InfoStructs');
if strcmp(label, 'neurons')
    folderResults = fullfile(folderPC, 'Lab\RESULTS\nonvisualEffects\modelGratingResp');
    folderRF = fullfile(folderPC, 'Lab\RESULTS\receptiveFields\SC neurons');
    
    corrections = [];
else
    folderResults = fullfile(folderPC, 'Lab\RESULTS\boutons\nonvisualEffects');
    folderRF = fullfile(folderPC, 'Lab\RESULTS\receptiveFields\boutons');
    
    data = load(fullfile(folderROIData, 'corrections_boutons.mat'));
    corrections = data.corrections;
    doCorrect = data.doCorrect;
end

data = load(fullfile(folderResults, 'pupil', 'tuning_prefDirSigmaDIFixed_lineFit.mat'));
tuning = data.tuning;
degrees = data.x;

data = load(fullfile(folderResults, 'kernelFit', 'results.mat'));
results = data.results;

expField = 'exp';
if ~isfield(results,expField)
    expField = 'expGratings';
end

%% Collect variables from visual noise
file = 'receptiveFields.mat';
% file = 'receptiveFields_separateFitEVs.mat';
data = load(fullfile(folderRF, file));
RFs = data.RFs;

datasets = [];
planes = [];
neurons = [];
evTotal = [];
evStim = [];
evRun = [];
lambdasStim = [];
pValues = [];
OnOffValues = [];

for k = 1:length(RFs)
    for iPlane = 1:length(RFs(k).plane)
        n = length(RFs(k).plane(iPlane).explainedVariances);
        datasets = [datasets; ones(n,1).*k];
        planes = [planes; ones(n,1).*iPlane];
        neurons = [neurons; (1:n)'];
        evTotal = [evTotal; RFs(k).plane(iPlane).explainedVariances];
        evStim = [evStim; RFs(k).plane(iPlane).explainedVariances_stimOnly];
        evRun = [evRun; RFs(k).plane(iPlane).explainedVariances_runOnly];
        lambdasStim = [lambdasStim; RFs(k).plane(iPlane).lambdasStim'];
        pValues = [pValues; RFs(k).plane(iPlane).pVal_RFonly];
        
        oov = NaN(n,2);
        for iCell = 1:n
            rf = RFs(k).plane(iPlane).receptiveFields(:,:,:,:,iCell);
            [~,t] = max(max(reshape(permute(abs(rf),[1 2 4 3]), [], ...
                size(rf,3)), [], 1));
            rf = squeeze(rf(:,:,t,:));
            [mx,row] = max(max(abs(rf),[],3),[],1);
            [~,col] = max(mx);
            row = row(col);
            rf = squeeze(rf(row,col,:));
            rf(2) = -rf(2);
            oov(iCell,:) = rf;
        end
        OnOffValues = [OnOffValues; oov];
    end
end

[~,type] = max(abs(OnOffValues),[],2);
ind = sub2ind(size(OnOffValues), (1:size(OnOffValues,1))', type);
signs = sign(OnOffValues(ind));
OnOffRatios = OnOffValues;
OnOffRatios(signs<0,:) = -OnOffRatios(signs<0,:);
OnOffRatios(OnOffRatios<0) = 0;
OnOffRatios = (OnOffRatios(:,1)-OnOffRatios(:,2))./sum(OnOffRatios,2);

validRF = pValues < minPVal & evStim > minExplainedVarianceStim & lambdasStim < maxLambda;

%% Determine identity of mouse for each unit
mice = NaN(size(datasets));
for k = 1:length(RFs)
    ind = datasets == k;
    if sum(ind)==0
        continue
    end
    mouse = str2double(RFs(k).subject(end-2:end));
    mice(ind) = mouse;
end

%% Collect variables of tuning curves
prefDirs = [];
responses = cell(1, length(RFs));
maxis = []; % resp to pref stimulus
stimMeans = [];
isSuppr = [];
isGad = [];

for k = 1:length(RFs)
    ktun = find(strcmp(RFs(k).subject, {tuning.subject}) & ...
        strcmp(RFs(k).date, {tuning.date}));
    if isempty(ktun) || isempty(tuning(ktun).plane)
        numCells = sum(datasets==k);
        prefDirs = [prefDirs; NaN(numCells,1)];
        responses{k} = NaN(numCells,12,1,2);
        maxis = [maxis; NaN(numCells,2)];
        stimMeans = [stimMeans; NaN(numCells,2)];
        isSuppr = [isSuppr; NaN(numCells,1)];
        isGad = [isGad; NaN(numCells,1)];
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
        numCells = sum(datasets==k & planes==iPlane);
        data = tuning(ktun).plane(iPlane);
        prefDir = NaN(numCells,1);
        resp = NaN(numCells, 12, numResp, 2);
        prefs = NaN(numCells, 2);
        means = NaN(numCells, 2);
        spp = NaN(numCells, 1);
        gd = NaN(numCells, 1);
        if isempty(data.cellIDs)
            prefDirs = [prefDirs; prefDir];
            responses{k} = [responses{k}; resp];
            maxis = [maxis; prefs];
            stimMeans = [stimMeans; means];
            isSuppr = [isSuppr; spp];
            isGad = [isGad; gd];
            continue
        end
        units = find(~cellfun(@isempty, {data.cond(1).cell.parameters}));
        for j = 1:length(units)
            iCell = units(j);
            for c = 1:2
                % NOTE: indexing with iCell is WRONG! It should be
                % data.cellIDs(iCell)!
                resp(iCell,:,:,c) = data.cond(c).cell(iCell).responses;
                
                pars = data.cond(c).cell(iCell).parameters;
                curve = data.cond(c).cell(iCell).curve;
                if length(pars) == 1 % not tuned
                    means(iCell,c) = pars;
                else
                    prefDir(iCell) = pars(1);
                    means(iCell,c) = mean(curve);
                    oris = mod(pars(1) + [0 90 180], 360);
                    sr = gratings.orituneWrappedConditions(pars, oris);
                    prefs(iCell,c) = sr(1);
                end
            end
            
            if ~isempty(corrections)
                a = corrections(ktun).plane(iPlane).a{tuning(ktun).exp} ...
                    (data.cellIDs(iCell));
                b = corrections(ktun).plane(iPlane).b{tuning(ktun).exp} ...
                    (data.cellIDs(iCell));
                resp(iCell,:,:,:) = doCorrect(a,b,resp(iCell,:,:,:));
                prefs(iCell,:) = doCorrect(a,b,prefs(iCell,:));
                means(iCell,:) = doCorrect(a,b,means(iCell,:));
            end
            
            spp(iCell) = data.isSuppressed(iCell);
            if isfield(data, 'isGad')
                gd(iCell) = data.isGad(iCell);
            end
        end
        prefDirs = [prefDirs; prefDir];
        responses{k} = [responses{k}; resp];
        maxis = [maxis; prefs];
        stimMeans = [stimMeans; means];
        isSuppr = [isSuppr; spp];
        isGad = [isGad; gd];
    end
end

%% Determine DSIs and OSIs
cond = 1; % small or large pupil
shuffles = 1000;
dirVectors = exp((0:30:330)./180.*pi .* 1i);
oriVectors = exp((0:30:330)./180.*2.*pi .* 1i);

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
meanResp = meanResp ./ sum(meanResp,2);
shuffleResp = shuffleResp ./ sum(shuffleResp,2);

% Determine DSIs
vects = sum(dirVectors .* meanResp, 2);
shuffleVects = squeeze(sum(dirVectors .* shuffleResp, 2));
DSIs = abs(vects);
nullDSIs = abs(shuffleVects);
p_DSI = sum(nullDSIs > DSIs,2) ./ shuffles;
p_DSI(isnan(DSIs)) = NaN;

% Determine OSIs
vects = sum(oriVectors .* meanResp, 2);
shuffleVects = squeeze(sum(oriVectors .* shuffleResp, 2));
OSIs = abs(vects);
nullOSIs = abs(shuffleVects);
p_OSI = sum(nullOSIs > OSIs,2) ./ shuffles;
p_OSI(isnan(OSIs)) = NaN;

%% Plot RF type (ON/OFF-ratio) for enhanced and suppressed
% check whether driven and suppressed units have different mean
% ON/OFF-ratios
ind = validRF;
tbl = table(OnOffRatios(ind), nominal(isSuppr(ind)), datasets(ind), mice(ind), ...
    'VariableNames', {'OnOff','isSuppr','session','mouse'});
lme = fitlme(tbl, 'OnOff ~ isSuppr + (isSuppr|session) + (isSuppr|mouse)','DummyVarCoding','effects');
lme = fitlme(tbl, 'OnOff ~ isSuppr + (1|session) + (isSuppr-1|session) + (1|mouse) + (isSuppr-1|mouse)','DummyVarCoding','effects');
lme = fitlme(tbl, 'OnOff ~ isSuppr + (1|session) + (1|mouse) + (isSuppr-1|mouse)','DummyVarCoding','effects');
lme = fitlme(tbl, 'OnOff ~ isSuppr + (isSuppr|session)','DummyVarCoding','effects');
lme = fitlme(tbl, 'OnOff ~ isSuppr + (1|session) + (isSuppr-1|session)','DummyVarCoding','effects')
lme = fitlme(tbl, 'OnOff ~ isSuppr + (isSuppr|mouse)','DummyVarCoding','effects')

binSize = 0.05;
edges = -1:binSize:1;
bins = edges(1:end-1)+binSize/2;
N1 = histcounts(OnOffRatios(isSuppr == -1 & validRF), edges);
N2 = histcounts(OnOffRatios(isSuppr == 1 & validRF), edges);
figure
b = bar(bins, [N2; N1], 'stacked');
b(1).FaceColor = 'k';
b(2).FaceColor = 'w';
legend('suppr-by-contrast','enhanced-by-contrast')
xlim([-1 1])
set(gca, 'box', 'off')
xlabel('(ON-OFF)/(ON+OFF)')
ylabel(sprintf('#%s',label))
title(sprintf('%s (n=%d)', label, sum(N1)+sum(N2)))

groups = validRF & [isSuppr == -1, isSuppr == 1];
groupNames = {'driven-by-contrast','suppressed-by-contrast'};
figure
hold
h = zeros(1, size(groups,2));
colors = lines(size(groups,2));
b = fixedEffects(lme);
means = [sum(b) b(1)-b(2)];
for g = 1:size(groups,2)
    x = sort(OnOffRatios(groups(:,g) & ~isnan(OnOffRatios)), 'ascend');
    y = (1:sum(groups(:,g) & ~isnan(OnOffRatios))) ./ ...
        sum(groups(:,g) & ~isnan(OnOffRatios));
    x = [-1; x; 1];
    y = [0 y 1];
    h(g) = plot(x, y, 'Color', colors(g,:), 'LineWidth', 2);
    plot(means(g), 1.1, 'v', 'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', colors(g,:))
end
xlim([min(OnOffRatios) max(OnOffRatios)])
ylim([0 1.2])
legend(h, groupNames)
xlabel('ON/OFF-ratio')
ylabel(['Proportion of ' label])
title(sprintf('%d driven, %d suppressed', ...
    sum(groups(:,1)&~isnan(OnOffRatios)), sum(groups(:,2)&~isnan(OnOffRatios))))

% % Test across datasets whether median ON-OFF-ratio is different between
% % enhanced and suppressed neurons
% dSets = unique(datasets);
% ratioMedians = NaN(length(dSets), size(groups,2));
% for d = 1:length(dSets)
%     for g = 1:size(groups,2)
%         ind = groups(:,g) & datasets==dSets(d);
%         if sum(ind) < minSamples
%             continue
%         end
%         ratioMedians(d,g) = nanmedian(OnOffRatios(ind));
%     end
% end
% ind = ~any(isnan(ratioMedians),2);
% p = signrank(ratioMedians(ind,1), ratioMedians(ind,2));
% fprintf('%d datasets, p = %.4f, medians: ', sum(ind), p)
% fprintf('%.4f ', median(ratioMedians(ind,:)))
% fprintf('\n')

%% Plot RF type (ON/OFF-ratio) for excitatory and inhibitory
% check whether excitatory and inhibitory units have different mean
% ON/OFF-ratios
ind = validRF & (isGad==1 | isGad==-1);
tbl = table(OnOffRatios(ind), nominal(isGad(ind)), datasets(ind), mice(ind), 'VariableNames', ...
    {'OnOff','inhibitory','session','mouse'});
lme = fitlme(tbl, 'OnOff ~ inhibitory + (inhibitory|session) + (inhibitory|mouse)','DummyVarCoding','effects');
lme = fitlme(tbl, 'OnOff ~ inhibitory + (1|session) + (inhibitory-1|session) + (1|mouse) + (inhibitory-1|mouse)','DummyVarCoding','effects');
lme = fitlme(tbl, 'OnOff ~ inhibitory + (1|session) + (1|mouse) + (inhibitory-1|mouse)','DummyVarCoding','effects');
lme = fitlme(tbl, 'OnOff ~ inhibitory + (1|session) + (1|mouse)','DummyVarCoding','effects');
lme = fitlme(tbl, 'OnOff ~ inhibitory + (inhibitory|session)','DummyVarCoding','effects');
lme = fitlme(tbl, 'OnOff ~ inhibitory + (1|session) + (inhibitory-1|session)','DummyVarCoding','effects');
lme = fitlme(tbl, 'OnOff ~ inhibitory + (1|session)','DummyVarCoding','effects');
lme = fitlme(tbl, 'OnOff ~ inhibitory + (1|mouse) + (inhibitory-1|mouse)','DummyVarCoding','effects');
lme = fitlme(tbl, 'OnOff ~ inhibitory + (1|mouse)','DummyVarCoding','effects');
lme = fitlme(tbl, 'OnOff ~ inhibitory','DummyVarCoding','effects')

groups = validRF & [isGad == -1, isGad == 1];
groupNames = {'excitatory','inhibitory'};
figure
hold
h = zeros(1, size(groups,2));
colors = lines(size(groups,2));
b = fixedEffects(lme);
means = [sum(b) b(1)-b(2)];
for g = 1:size(groups,2)
    x = sort(OnOffRatios(groups(:,g) & ~isnan(OnOffRatios)), 'ascend');
    y = (1:sum(groups(:,g) & ~isnan(OnOffRatios))) ./ ...
        sum(groups(:,g) & ~isnan(OnOffRatios));
    x = [-1; x; 1];
    y = [0 y 1];
    h(g) = plot(x, y, 'Color', colors(g,:), 'LineWidth', 2);
    plot(means(g), 1.1, 'v', 'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', colors(g,:))
end
xlim([min(OnOffRatios) max(OnOffRatios)])
ylim([0 1.2])
legend(h, groupNames)
xlabel('ON/OFF-ratio')
ylabel(['Proportion of ' label])
title(sprintf('%d %s, %d %s', sum(groups(:,1)&~isnan(OnOffRatios)), ...
    groupNames{1}, sum(groups(:,2)&~isnan(OnOffRatios)), groupNames{2}))

%% Plot RF type versus preferred direction
figure
scatter(prefDirs(validRF), OnOffRatios(validRF), [], 'k', 'filled', 'MarkerFaceAlpha', 0.2)
axis([0 360 -1 1])
set(gca, 'box', 'off', 'XTick', 0:90:360, 'YTick', -1:1)
xlabel('Preferred direction')
ylabel('ON/OFF-ratio')
title(sprintf('%s (n=%d)', label, sum(~any(isnan([prefDirs(validRF) OnOffRatios(validRF)]),2))))

binSize = 90;
edges = -binSize/2:binSize:360;
bins = edges(1:end-1)+binSize/2;
pd = prefDirs;
pd(pd>=360-binSize/2) = pd(pd>=360-binSize/2) - 360;
[N1,~,b1] = histcounts(pd(validRF & OnOffRatios > onThresh & ~isnan(pd)), edges);
[N2,~,b2] = histcounts(pd(validRF & OnOffRatios < offThresh & ~isnan(pd)), edges);
[N3,~,b3] = histcounts(pd(validRF & OnOffRatios >= offThresh & ...
    OnOffRatios <= onThresh & ~isnan(pd)), edges);
[N4,~,b4] = histcounts(pd(~validRF & ~isnan(pd)), edges);
figure
hold on
plot([bins 360], N1([1:end 1]) ./ sum(N1), 'LineWidth', 1)
plot([bins 360], N2([1:end 1]) ./ sum(N2), 'LineWidth', 1)
plot([bins 360], N3([1:end 1]) ./ sum(N3), 'LineWidth', 1)
plot([bins 360], N4([1:end 1]) ./ sum(N4), 'LineWidth', 1)
legend(sprintf('ON (n=%d)',sum(N1)), sprintf('OFF (n=%d)',sum(N2)), ...
    sprintf('ON+OFF (n=%d)',sum(N3)), sprintf('no RF (n=%d)',sum(N4)))
xlim([0 360])
set(gca, 'XTick', 0:90:360, 'box', 'off')
xlabel('Preferred direction')
ylabel('Proportion of boutons')

% [p,tbl,stats] = anova1([b1;b2;b3], [ones(length(b1),1); ones(length(b2),1).*2; ...
%     ones(length(b3),1).*3]);

%% Plot DSI/OSI for each RF type (ON, OFF, ON+OFF)
groups = [validRF & OnOffRatios>onThresh, ...
    validRF & OnOffRatios<offThresh, ...
    validRF & OnOffRatios<=onThresh & OnOffRatios>=offThresh];
groupNames = {'ON', 'OFF', 'ON+OFF'};

RFtypes = cell(length(validRF),1);
for g = 1:3
    [RFtypes{groups(:,g)}] = deal(groupNames{g});
end
% check whether ON, OFF, ON+OFF units have different mean DSI
tbl = table(DSIs(validRF), RFtypes(validRF), datasets(validRF), mice(validRF), 'VariableNames', ...
    {'DSI','RF','session','mouse'});
lme = fitlme(tbl, 'DSI ~ RF + (RF|session) + (RF|mouse)', 'DummyVarCoding', 'effects');
lme = fitlme(tbl, 'DSI ~ RF + (RF|session)', 'DummyVarCoding', 'effects')
lme = fitlme(tbl, 'DSI ~ RF + (1|session) + (RF-1|session)', 'DummyVarCoding', 'effects')
lme = fitlme(tbl, 'DSI ~ RF + (RF|session)', 'DummyVarCoding', 'effects')
lme = fitlme(tbl, 'DSI ~ RF + (RF|mouse)', 'DummyVarCoding', 'effects')
lme = fitlme(tbl, 'DSI ~ RF + (1|mouse) + (RF-1|mouse)', 'DummyVarCoding', 'effects')
lme = fitlme(tbl, 'DSI ~ RF + (RF-1|mouse)', 'DummyVarCoding', 'effects')
anova(lme)

[b, bnames] = fixedEffects(lme);
bnames = bnames.Variables;
means = NaN(1, size(groups,2));
figure;
hold
h = zeros(1, size(groups,2));
lbls = cell(1, size(groups,2));
colors = lines(size(groups,2));
for g = 1:size(groups,2)
    ind = groups(:,g) & ~isnan(DSIs);
    x = sort(DSIs(ind), 'ascend');
    y = (1:sum(ind)) ./ sum(ind);
    x = [0; x; 1];
    y = [0 y 1];
    h(g) = plot(x, y, 'Color', colors(g,:), 'LineWidth', 2);
    v = find(strcmp(bnames, ['RF_' groupNames{g}]));
    if ~isempty(v)
        means(g) = b(1) + b(v);
    else
        means(g) = b(1) - sum(b(2:end));
    end
    plot(means(g), 1.1, 'v', 'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', colors(g,:))
    lbls{g} = sprintf('%s (n=%d)', groupNames{g}, sum(ind));
end
xlim([0 max(DSIs)])
legend(h, lbls)
xlabel('DSI')
ylabel(['Proportion of ' label])


% check whether ON, OFF, ON+OFF units have different mean OSI
tbl = table(OSIs(validRF), RFtypes(validRF), datasets(validRF), mice(validRF), 'VariableNames', ...
    {'OSI','RF','session','mouse'});
lme = fitlme(tbl, 'OSI ~ RF + (RF|session) + (RF|mouse)', 'DummyVarCoding', 'effects');
lme = fitlme(tbl, 'OSI ~ RF + (RF|session)', 'DummyVarCoding', 'effects')
lme = fitlme(tbl, 'OSI ~ RF + (1|session) + (RF-1|session)', 'DummyVarCoding', 'effects')
lme = fitlme(tbl, 'OSI ~ RF + (1|session)', 'DummyVarCoding', 'effects')
lme = fitlme(tbl, 'OSI ~ RF + (RF|mouse)', 'DummyVarCoding', 'effects')
lme = fitlme(tbl, 'OSI ~ RF + (1|mouse) + (RF-1|mouse)', 'DummyVarCoding', 'effects')
lme = fitlme(tbl, 'OSI ~ RF', 'DummyVarCoding', 'effects')
anova(lme)

[b, bnames] = fixedEffects(lme);
bnames = bnames.Variables;
means = NaN(1, size(groups,2));
figure;
hold
h = zeros(1, size(groups,2));
lbls = cell(1, size(groups,2));
colors = lines(size(groups,2));
for g = 1:size(groups,2)
    ind = groups(:,g) & ~isnan(OSIs);
    x = sort(OSIs(ind), 'ascend');
    y = (1:sum(ind)) ./ sum(ind);
    x = [0; x; 1];
    y = [0 y 1];
    h(g) = plot(x, y, 'Color', colors(g,:), 'LineWidth', 2);
    v = find(strcmp(bnames, ['RF_' groupNames{g}]));
    if ~isempty(v)
        means(g) = b(1) + b(v);
    else
        means(g) = b(1) - sum(b(2:end));
    end
    plot(means(g), 1.1, 'v', 'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', colors(g,:))
    lbls{g} = sprintf('%s (n=%d)', groupNames{g}, sum(ind));
end
xlim([0 max(OSIs)])
legend(h, lbls)
xlabel('OSI')
ylabel(['Proportion of ' label])

%% Plot DSI/OSI for excitatory and inhibitory units
groups = [isGad==-1, isGad==1];
groupNames = {'excitatory','inhibitory'};

% check whether excitatory and inhibitory units have different mean DSI
ind = groups(:,1) | groups(:,2);
tbl = table(DSIs(ind), nominal(isGad(ind)), datasets(ind), mice(ind), 'VariableNames', ...
    {'DSI','inhibitory','session','mouse'});
lme = fitlme(tbl, 'DSI ~ inhibitory + (inhibitory|session) + (inhibitory|mouse)', 'DummyVarCoding', 'effects');
lme = fitlme(tbl, 'DSI ~ inhibitory + (inhibitory|session)', 'DummyVarCoding', 'effects')
lme = fitlme(tbl, 'DSI ~ inhibitory + (1|session) + (inhibitory-1|session)', 'DummyVarCoding', 'effects')

b = fixedEffects(lme);
means = [sum(b) b(1)-b(2)];
figure;
hold
h = zeros(1, size(groups,2));
lbls = cell(1, size(groups,2));
colors = lines(size(groups,2));
for g = 1:size(groups,2)
    ind = groups(:,g) & ~isnan(DSIs);
    x = sort(DSIs(ind), 'ascend');
    y = (1:sum(ind)) ./ sum(ind);
    x = [0; x; 1];
    y = [0 y 1];
    h(g) = plot(x, y, 'Color', colors(g,:), 'LineWidth', 2);
    plot(means(g), 1.1, 'v', 'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', colors(g,:))
    lbls{g} = sprintf('%s (n=%d)', groupNames{g}, sum(ind));
end
xlim([0 max(DSIs)])
legend(h, lbls)
xlabel('DSI')
ylabel(['Proportion of ' label])


% check whether excitatory and inhibitory units have different mean OSI
ind = groups(:,1) | groups(:,2);
tbl = table(OSIs(ind), nominal(isGad(ind)), datasets(ind), mice(ind), 'VariableNames', ...
    {'OSI','inhibitory','session','mouse'});
lme = fitlme(tbl, 'OSI ~ inhibitory + (inhibitory|session) + (inhibitory|mouse)', 'DummyVarCoding', 'effects');
lme = fitlme(tbl, 'OSI ~ inhibitory + (inhibitory|session)', 'DummyVarCoding', 'effects')
lme = fitlme(tbl, 'OSI ~ inhibitory + (1|session) + (inhibitory-1|session)', 'DummyVarCoding', 'effects')
lme = fitlme(tbl, 'OSI ~ inhibitory + (1|session) + (inhibitory-1|session) + (1|mouse)', 'DummyVarCoding', 'effects')

b = fixedEffects(lme);
means = [sum(b) b(1)-b(2)];
figure;
hold
h = zeros(1, size(groups,2));
lbls = cell(1, size(groups,2));
colors = lines(size(groups,2));
for g = 1:size(groups,2)
    ind = groups(:,g) & ~isnan(DSIs);
    x = sort(OSIs(ind), 'ascend');
    y = (1:sum(ind)) ./ sum(ind);
    x = [0; x; 1];
    y = [0 y 1];
    h(g) = plot(x, y, 'Color', colors(g,:), 'LineWidth', 2);
    plot(means(g), 1.1, 'v', 'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', colors(g,:))
    lbls{g} = sprintf('%s (n=%d)', groupNames{g}, sum(ind));
end
xlim([0 max(OSIs)])
legend(h, lbls)
xlabel('OSI')
ylabel(['Proportion of ' label])

%% Plot preferred direction versus DSI/OSI
% only consider neurons that have a significant DSI or OSI
figure
hold on
scatter(prefDirs(p_DSI>=.05&p_OSI<.05), DSIs(p_DSI>=.05&p_OSI<.05), [], ...
    'k', 'filled', 'MarkerFaceAlpha', 0.2)
scatter(prefDirs(p_DSI<.05), DSIs(p_DSI<.05), [], ...
    'r', 'filled', 'MarkerFaceAlpha', 0.2)
legend('DSI: p>0.05','DSI: p<0.05')
axis([0 360 0 .32])
set(gca, 'box', 'off', 'XTick', 0:90:360, 'YTick', 0:.1:.4)
xlabel('Preferred direction')
ylabel('DSI')
title(sprintf('%s (n=%d)', label, sum(~any(isnan([DSIs prefDirs]),2) & ...
    (p_DSI<.05 | p_OSI<.05))))

figure
hold on
scatter(prefDirs(p_OSI>=.05&p_DSI<.05), OSIs(p_OSI>=.05&p_DSI<.05), [], ...
    'k', 'filled', 'MarkerFaceAlpha', 0.2)
scatter(prefDirs(p_OSI<.05), OSIs(p_OSI<.05), [], ...
    'r', 'filled', 'MarkerFaceAlpha', 0.2)
legend('OSI: p>0.05','OSI: p<0.05')
axis([0 360 0 .32])
set(gca, 'box', 'off', 'XTick', 0:90:360, 'YTick', 0:.1:.4)
xlabel('Preferred direction')
ylabel('OSI')
title(sprintf('%s (n=%d)', label, sum(~any(isnan([OSIs prefDirs]),2) & ...
    (p_DSI<.05 | p_OSI<.05))))

binSize = 10;
edges = 0:binSize:360;
bins = edges(1:end-1)+binSize/2;
figure
N1 = histcounts(prefDirs(p_DSI<.05&p_OSI>=.05), edges);
N2 = histcounts(prefDirs(p_DSI>=.05&p_OSI<.05), edges);
N3 = histcounts(prefDirs(p_DSI<.05&p_OSI<.05), edges);
hold on
plot(bins, N1./sum(N1), 'r', 'LineWidth', 2)
plot(bins, N2./sum(N2), 'b', 'LineWidth', 2)
plot(bins, N3./sum(N3), 'Color', [.5 0 .5], 'LineWidth', 2)
legend('DSI: p<0.05','OSI: p<0.05','DSI & OSI: p<0.05')
xlim([0 360])
set(gca, 'box', 'off', 'XTick', 0:90:360)
xlabel('Preferred direction')
ylabel('Proportion of cell type')
title(sprintf('%s (n=%d)', label, sum([N1,N2,N3])))

%% Plot preferred direction versus enhanced/suppressed
binSize = 10;
edges = 0:binSize:360;
bins = edges(1:end-1)+binSize/2;
N1 = histcounts(prefDirs(isSuppr == -1), edges);
N2 = histcounts(prefDirs(isSuppr == 1), edges);
figure
b = bar(bins, [N2; N1], 'stacked');
b(1).FaceColor = 'k';
b(2).FaceColor = 'w';
legend('suppr-by-contrast','enhanced-by-contrat')
set(gca, 'box', 'off')
xlabel('Preferred direction (deg)')
ylabel(sprintf('#%s',label))
title(sprintf('%s (n=%d)', label, sum(N1)+sum(N2)))

%% Compare enhanced/driven by gratings to sign of RF
ind = validRF & ~isnan(isSuppr);
N = sum([isSuppr(ind)==-1 & signs(ind)==1, isSuppr(ind)==-1 & signs(ind)==-1, ...
    isSuppr(ind)==1 & signs(ind)==1, isSuppr(ind)==1 & signs(ind)==-1]);
figure
bar(N)
set(gca, 'XTickLabel', {'grat+ RF+','grat+ RF-','grat- RF+','grat- RF-'})
fprintf('%.2f%% %s driven by gratings were driven by visual noise\n', N(1)/sum(N(1:2)).*100, label)
fprintf('%.2f%% %s suppressed by gratings were suppressed by visual noise\n', N(4)/sum(N(3:4)).*100, label)

%% Plot preferred direction versus excitatory/inhibitory
binSize = 10;
edges = 0:binSize:360;
bins = edges(1:end-1)+binSize/2;
N1 = histcounts(prefDirs(isGad == -1), edges);
N2 = histcounts(prefDirs(isGad == 1), edges);
figure
b = bar(bins, [N1; N2], 'stacked');
b(1).FaceColor = 'k';
b(2).FaceColor = 'w';
legend(sprintf('excitatory (%d)',sum(N1)),sprintf('inhibitory (%d)',sum(N2)))
set(gca, 'box', 'off')
xlabel('Preferred direction (deg)')
ylabel(sprintf('#%s',label))
title(sprintf('%s (n=%d)', label, sum(N1)+sum(N2)))
