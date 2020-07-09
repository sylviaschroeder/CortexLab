%% Parameters
stimuli = {'gratings', 'grayScreen', 'dark'};
% label = 'neurons';
label = 'boutons';

% for evaluation of receptive fields (significance/goodness)
minExplainedVarianceStim = 0.01;
minPVal = 0.05;
maxLambda = 1;
onThr = 0.5; % threshold for ON cells
offThr = -0.5; % threshold for OFF cells

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
        
        n = n + numCells;
    end 
end

pValues_rhoRun = sum(nullRhosRun > rhosRun,3) ./ size(nullRhosRun,3);
pValues_rhoRun(pValues_rhoRun > 0.5) = 1 - pValues_rhoRun(pValues_rhoRun > 0.5);
pValues_rhoRun = pValues_rhoRun .* 2;
pValues_rhoRun(isnan(rhosRun)) = NaN;
pValues_rhoPupil = sum(nullRhosPupil > rhosPupil,3) ./ size(nullRhosPupil,3);
pValues_rhoPupil(pValues_rhoPupil > 0.5) = 1 - pValues_rhoPupil(pValues_rhoPupil > 0.5);
pValues_rhoPupil = pValues_rhoPupil .* 2;
pValues_rhoPupil(isnan(rhosPupil)) = NaN;

%% Determine identity of mouse for each unit
mice = NaN(size(dataset));
for k = 1:length(corrsRun)
    ind = dataset == k;
    if sum(ind)==0
        continue
    end
    mouse = str2double(corrsRun(k).subject(end-2:end));
    mice(ind) = mouse;
end

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
minis = []; % resp to non-pref stimulus
maxis = []; % resp to pref stimulus
stimMeans = [];
nullMinis = [];
nullMaxis = [];
nullStimMeans = [];
isSuppr = [];
isGad = [];
isTuned = [];
explVar = [];
draws = cellfun(@size,squeeze(struct2cell(null(1).plane(1).cond(1).cell)),'UniformOutput',false);
draws = max(cell2mat(draws),[],1);
draws = draws(1);

behCondStim = [];
directionStim = [];

fprintf('Dataset (of %d): ', length(corrsRun))
for k = 1:length(corrsRun)
    fprintf('%d ', k)
    ktun = find(strcmp(corrsRun(k).subject, {tuning.subject}) & ...
        strcmp(corrsRun(k).date, {tuning.date}));
    if isempty(ktun) || isempty(tuning(ktun).plane)
        numCells = sum(dataset==k);
        prefDirs = [prefDirs; NaN(numCells,1)];
        responses{k} = NaN(numCells,12,1,2);
        minis = [minis; NaN(numCells,2)];
        maxis = [maxis; NaN(numCells,2)];
        stimMeans = [stimMeans; NaN(numCells,2)];
        nullMinis = [nullMinis; NaN(numCells,2,draws)];
        nullMaxis = [nullMaxis; NaN(numCells,2,draws)];
        nullStimMeans = [nullStimMeans; NaN(numCells,2,draws)];
        isSuppr = [isSuppr; NaN(numCells,1)];
        isGad = [isGad; NaN(numCells,1)];
        isTuned = [isTuned; NaN(numCells,1)];
        explVar = [explVar; NaN(numCells,1)];
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
        mn = NaN(numCells,2);
        mx = NaN(numCells,2);
        sm = NaN(numCells,2);
        nmn = NaN(numCells,2,draws);
        nmx = NaN(numCells,2,draws);
        nsm = NaN(numCells,2,draws);
        spp = NaN(numCells,1);
        gad = NaN(numCells,1);
        tun = NaN(numCells,1);
        eV = NaN(numCells,1);
        if isempty(data.cellIDs)
            prefDirs = [prefDirs; pd];
            responses{k} = [responses{k}; resp];
            minis = [minis; mn];
            maxis = [maxis; mx];
            stimMeans = [stimMeans; sm];
            nullMinis = [nullMinis; nmn];
            nullMaxis = [nullMaxis; nmx];
            nullStimMeans = [nullStimMeans; nsm];
            isSuppr = [isSuppr; spp];
            isGad = [isGad; gad];
            isTuned = [isTuned; tun];
            explVar = [explVar; eV];
            continue
        end
        
        if iPlane == 1
            % get condition of trials according to pupil size (large or
            % small)
            validPlanes = find(~cellfun(@isempty, {tuning(ktun).plane.cellIDs}));
            validROI = [];
            for p = 1:length(validPlanes)
                validROI = find(~isnan(tuning(ktun).plane(validPlanes(p)).isTuned), 1);
                if ~isempty(validROI)
                    break
                end
            end
            if ~isempty(validROI)
                conditions = NaN(size(tuning(ktun).plane(validPlanes(p)) ...
                    .cond(1).cell(validROI).responses));
                for c = 1:2
                    conditions(~isnan(tuning(ktun).plane(validPlanes(p)).cond(c) ...
                        .cell(validROI).responses)) = c;
                end
                behCondStim = [behCondStim; conditions(:)];
                directionStim = [directionStim; ...
                    repmat(tuning(ktun).plane(validPlanes(p)).cond(c) ...
                    .cell(validROI).directions, size(conditions,2), 1)];
            end
        end
        
        spp(data.cellIDs) = data.isSuppressed;
        if isfield(data, 'isGad')
            gad(data.cellIDs) = data.isGad;
        end
        tun(data.cellIDs) = data.isTuned;
        eV(data.cellIDs) = data.crossValExplVar;
        neurons = find(~cellfun(@isempty, {data.cond(1).cell.parameters}));
        for j = 1:length(neurons)
            iCell = neurons(j);
            ID = data.cellIDs(iCell);
            for c = 1:2
                resp(ID,:,:,c) = data.cond(c).cell(iCell).responses;
                pars = data.cond(c).cell(iCell).parameters;
                curve = data.cond(c).cell(iCell).curve;
                if length(pars) == 1 % not tuned
                    sm(ID,c) = pars;
                else
                    if c ==1
                        pd(ID) = pars(1);
                    end
                    sm(ID,c) = mean(curve);
                    oris = mod(pars(1) + [0 90 180], 360);
                    sr = gratings.orituneWrappedConditions(pars, oris);
                    mx(ID,c) = sr(1);
                    if sr(1)-sr(2)>0
                        [mn(ID,c), ind] = min(sr(2:3));
                    else % suppressed neurons
                        [mn(ID,c), ind] = max(sr(2:3));
                    end
                    ind = ind+1;
                end
                
                pars = null(ktun).plane(iPlane).cond(c).cell(iCell).parameters;
                if size(pars,2) == 1 % not tuned
                    nsm(ID,c,:) = pars;
                else
                    sr = NaN(size(pars,1), 3);
                    curves = NaN(size(pars,1), length(degrees));
                    for p = 1:size(pars,1)
                        oris = mod(pars(p,1) + [0 90 180], 360);
                        sr(p,:) = gratings.orituneWrappedConditions(pars(p,:), oris);
                        curves(p,:) = gratings.orituneWrappedConditions(pars(p,:), degrees);
                    end
                    nsm(ID,c,:) = mean(curves,2);
                    nmx(ID,c,:) = sr(:,1);
                    nmn(ID,c,:) = sr(:,ind);
                end
            end
            
            if ~isempty(corrections)
                a = corrections(ktun).plane(iPlane).a{tuning(ktun).exp} ...
                    (data.cellIDs(iCell));
                b = corrections(ktun).plane(iPlane).b{tuning(ktun).exp} ...
                    (data.cellIDs(iCell));
                resp(ID,:,:,:) = doCorrect(a,b,resp(ID,:,:,:));
                mn(ID,:) = doCorrect(a,b,mn(ID,:));
                mx(ID,:) = doCorrect(a,b,mx(ID,:));
                sm(ID,:) = doCorrect(a,b,sm(ID,:));
                nmn(ID,:,:) = doCorrect(a,b,nmn(ID,:,:));
                nmx(ID,:,:) = doCorrect(a,b,nmx(ID,:,:));
                nsm(ID,:,:) = doCorrect(a,b,nsm(ID,:,:));
            end
        end
        prefDirs = [prefDirs; pd];
        responses{k} = [responses{k}; resp];
        minis = [minis; mn];
        maxis = [maxis; mx];
        stimMeans = [stimMeans; sm];
        nullMinis = [nullMinis; nmn];
        nullMaxis = [nullMaxis; nmx];
        nullStimMeans = [nullStimMeans; nsm];
        isSuppr = [isSuppr; spp];
        isGad = [isGad; gad];
        isTuned = [isTuned; tun];
        explVar = [explVar; eV];
    end
end
fprintf('\n')

mods = abs(maxis - minis);
nullMods = abs(nullMaxis - nullMinis);

ind = all(isnan(maxis),2); % untuned neurons
mx = maxis;
mx(ind,:) = stimMeans(ind,:);
nmx = nullMaxis;
nmx(ind,:,:) = nullStimMeans(ind,:,:);
modFun = @(a,b) (b-a)./(abs(a)+abs(b));
DImaxis = modFun(maxis(:,1), maxis(:,2));
nullDImaxis = modFun(squeeze(nullMaxis(:,1,:)), squeeze(nullMaxis(:,2,:)));
pValues_DImaxis = sum(DImaxis > nullDImaxis, 2) ./ size(nullDImaxis,2);
pValues_DImaxis(pValues_DImaxis > 0.5) = 1 - pValues_DImaxis(pValues_DImaxis > 0.5);
pValues_DImaxis = pValues_DImaxis .* 2;
DIminis = modFun(minis(:,1), minis(:,2));
nullDIminis = modFun(squeeze(nullMinis(:,1,:)), squeeze(nullMinis(:,2,:)));
pValues_DIminis = sum(DIminis > nullDIminis, 2) ./ size(nullDIminis,2);
pValues_DIminis(pValues_DIminis > 0.5) = 1 - pValues_DIminis(pValues_DIminis > 0.5);
pValues_DIminis = pValues_DIminis .* 2;
DImaxisAll = modFun(mx(:,1), mx(:,2));
nullDImaxisAll = modFun(squeeze(nmx(:,1,:)), squeeze(nmx(:,2,:)));
pValues_DImaxisAll = sum(DImaxisAll > nullDImaxisAll, 2) ./ size(nullDImaxisAll,2);
pValues_DImaxisAll(pValues_DImaxisAll > 0.5) = 1 - pValues_DImaxisAll(pValues_DImaxisAll > 0.5);
pValues_DImaxisAll = pValues_DImaxisAll .* 2;
DImods = modFun(mods(:,1), mods(:,2));
nullDImods = modFun(squeeze(nullMods(:,1,:)), squeeze(nullMods(:,2,:)));
pValues_DImods = sum(DImods > nullDImods, 2) ./ size(nullDImods,2);
pValues_DImods(pValues_DImods > 0.5) = 1 - pValues_DImods(pValues_DImods > 0.5);
pValues_DImods = pValues_DImods .* 2;

% for suppressed units, change sign of DIs of preferred responses so that
% a more negative response with arousal results in a positive DI
DImaxis(isSuppr==1) = -DImaxis(isSuppr==1);
nullDImaxis(isSuppr==1,:) = -nullDImaxis(isSuppr==1,:);
DIminis(isSuppr==1) = -DIminis(isSuppr==1);
nullDIminis(isSuppr==1,:) = -nullDIminis(isSuppr==1,:);
DImaxisAll(isSuppr==1) = -DImaxisAll(isSuppr==1);
nullDImaxisAll(isSuppr==1,:) = -nullDImaxisAll(isSuppr==1,:);

isResponsive = ~isnan(explVar);

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

%% Determine DSIs and OSIs
directions = 0:30:330;
shuffles = 1000;
dirVectors = exp(directions./180.*pi .* 1i);
oriVectors = exp(directions./180.*2.*pi .* 1i);

meanResp = cell(1,2);
shuffleResp = cell(1,2);
for cond = 1:2 % small or large pupil
    for k = 1:length(responses)
        r = responses{k}(:,:,:,cond);
        meanResp{cond} = [meanResp{cond}; nanmean(r, 3)];
        r1 = reshape(r, size(r,1), []);
        shR = NaN(size(r,1),length(dirVectors),shuffles);
        for sh = 1:shuffles
            p = randperm(size(r1,2));
            r2 = reshape(r1(:,p), size(r));
            shR(:,:,sh) = nanmean(r2, 3);
        end
        shuffleResp{cond} = [shuffleResp{cond}; shR];
    end
end
meanResp = cat(3, meanResp{:}); % [neuron x directions x small/large pupil]
shuffleResp = cat(4, shuffleResp{:}); % [neuron x directions x draws x small/large pupil]

meanResp(isSuppr==1,:,:) = -meanResp(isSuppr==1,:,:);
shuffleResp(isSuppr==1,:,:,:) = -shuffleResp(isSuppr==1,:,:,:);
meanResp(meanResp<0) = 0;
shuffleResp(shuffleResp<0) = 0;
meanResp = meanResp ./ sum(meanResp,2);
shuffleResp = shuffleResp ./ sum(shuffleResp,2);

% Determine DSIs
vects = squeeze(sum(dirVectors .* meanResp, 2)); % [neuron x small/large pupil]
shuffleVects = squeeze(sum(dirVectors .* shuffleResp, 2)); % [neuron x draws x small/large pupil]
DSIs = abs(vects); % [neuron x small/large pupil]
nullDSIs = abs(shuffleVects);
p_DSI = squeeze(sum(nullDSIs > permute(DSIs,[1 3 2]),2) ./ shuffles);
p_DSI(isnan(DSIs)) = NaN; % [neuron x small/large pupil]

% Determine OSIs
vects = squeeze(sum(oriVectors .* meanResp, 2)); % [neuron x small/large pupil]
shuffleVects = squeeze(sum(oriVectors .* shuffleResp, 2)); % [neuron x draws x small/large pupil]
OSIs = abs(vects); % [neuron x small/large pupil]
nullOSIs = abs(shuffleVects);
p_OSI = squeeze(sum(nullOSIs > permute(OSIs,[1 3 2]),2) ./ shuffles);
p_OSI(isnan(OSIs)) = NaN; % [neuron x small/large pupil]

% Determine pref directions
vects = squeeze(sum(dirVectors .* meanResp, 2)); % [neuron x small/large pupil]
prefDirs_vectAvg = mod(angle(vects) ./ pi .* 180, 360);

% Determine pref orientations
vects = squeeze(sum(oriVectors .* meanResp, 2)); % [neuron x small/large pupil]
prefOris_vectAvg = mod(angle(vects) ./ pi .* 180, 360) ./ 2;

%% Get percentages of units: tuned to dir./ori., ...
% direction or orientation selective:
% (1) based on tuning curve fits
sum(isTuned==1) / sum(isResponsive)
% (2) based on DSI/OSI
sum(p_DSI(:,1)<0.05|p_OSI(:,1)<0.05) / sum(isResponsive)

%% Compare DSI/OSI and preferred directions/orientations between small and large pupil
edges = 0 : 1/100 : 1;
bins = edges(1:end-1)+diff(edges(1:2));
m = hot(200);
m = m(1:180,:);

% Compare DSIs (all neurons)
% use linear mixed-effects model to estimate slope and its significance
tbl = table(DSIs(:,1), DSIs(:,2), mice, dataset, 'VariableNames', ...
    {'DSI_small', 'DSI_large', 'mouse', 'session'});
lme = fitlme(tbl, 'DSI_large ~ -1 + DSI_small + (-1 + DSI_small | session) + (-1 + DSI_small | mouse)');

N = histcounts2(DSIs(:,1), DSIs(:,2), edges, edges);
[xout, yout, zout] = prepareSurfaceData(bins, bins, N);
f = fit([xout, yout], zout, 'linearinterp');
densities = f(DSIs);
figure
hold on
scatter(DSIs(:,1), DSIs(:,2), [], densities, 'filled')
plot([0 1], [0 1], 'k', 'LineWidth', 1)
plot([0 1], [0 1] .* fixedEffects(lme), 'r', 'LineWidth', 1)
colormap(m)
colorbar
axis square
axis([0 1 0 1])
xlabel('DSI small pupil')
ylabel('DSI large pupil')
title(sprintf('%s (n=%d, slope: %.3f, p = %.2e)',label, ...
    sum(~any(isnan(DSIs),2)), fixedEffects(lme), coefTest(lme)))

% Compare OSIs (all neurons)
% use linear mixed-effects model to estimate slope and its significance
tbl = table(OSIs(:,1), OSIs(:,2), mice, dataset, 'VariableNames', ...
    {'OSI_small', 'OSI_large', 'mouse', 'session'});
lme = fitlme(tbl, 'OSI_large ~ -1 + OSI_small + (-1 + OSI_small | session) + (-1 + OSI_small | mouse)');

N = histcounts2(OSIs(:,1), OSIs(:,2), edges, edges);
[xout, yout, zout] = prepareSurfaceData(bins, bins, N);
f = fit([xout, yout], zout, 'linearinterp');
densities = f(OSIs);
figure
hold on
scatter(OSIs(:,1), OSIs(:,2), [], densities, 'filled')
plot([0 1], [0 1], 'k', 'LineWidth', 1)
plot([0 1], [0 1] .* fixedEffects(lme), 'r', 'LineWidth', 1)
colormap(m)
colorbar
axis square
axis([0 1 0 1])
xlabel('OSI small pupil')
ylabel('OSI large pupil')
title(sprintf('%s (n=%d, slope: %.3f, p = %.2e)',label, ...
    sum(~any(isnan(DSIs),2)), fixedEffects(lme), coefTest(lme)))

% Compare preferred directions for DS neurons (p < 0.05 during small OR
% large pupil)
ind = all(p_DSI < 0.05,2);
a = prefDirs_vectAvg(ind,1)./180.*pi - prefDirs_vectAvg(ind,2)./180.*pi;
ind = a > pi;
a(ind) = 2*pi - a(ind);
ind = a < - pi;
a(ind) = 2*pi + a(ind);
circ_mtest(a, 0, 1e-100)

ind = all(p_DSI < 0.05,2);
figure
hold on
scatter(prefDirs_vectAvg(ind,1), prefDirs_vectAvg(ind,2), 'filled')
plot([0 360], [0 360], 'k', 'LineWidth', 1)
axis square
axis([0 360 0 360])
set(gca, 'XTick', 0:90:360, 'YTick', 0:90:360)
xlabel('Preferred direction, small pupil')
ylabel('Preferred direction, large pupil')
title(sprintf('%s (n=%d)',label, sum(ind)))

% Compare preferred orientations for OS neurons (p < 0.05 during small OR
% large pupil)
ind = all(p_OSI < 0.05,2);
a = prefOris_vectAvg(ind,1)./180.*2.*pi - prefOris_vectAvg(ind,2)./180.*2.*pi;
ind = a > pi;
a(ind) = 2*pi - a(ind);
ind = a < - pi;
a(ind) = 2*pi + a(ind);
circ_mtest(a, 0, 1e-10)

figure
hold on
ind = all(p_OSI < 0.05,2);
scatter(prefOris_vectAvg(ind,1), prefOris_vectAvg(ind,2), 'filled')
plot([0 180], [0 180], 'k', 'LineWidth', 1)
axis square
axis([0 180 0 180])
set(gca, 'XTick', 0:45:180, 'YTick', 0:45:180)
xlabel('Preferred orientation, small pupil')
ylabel('Preferred orientation, large pupil')
title(sprintf('%s (n=%d)',label, sum(ind)))

%% Plot cumulative distributions of correlations with pupil/running for cell types
rhos = {rhosRun, rhosPupil};
nullRhos = {nullRhosRun, nullRhosPupil};
nonvisNames = {'running','pupil'};
if strcmp(label, 'neurons')
    neuronSets = {{isGad==-1; isGad==1&isSuppr==-1; isGad==1&isSuppr==1}, ... % exc, inh&enh, inh&suppr
        {OnOffRatios>onThr & validRF; ...
        OnOffRatios>=offThr & OnOffRatios<=onThr & validRF; ...
        OnOffRatios<offThr & validRF}, ... % ON, Absolute, OFF
        {any(p_DSI<.05,2) & all(p_OSI>=0.05,2); ...
        any(p_OSI<.05,2) & all(p_DSI>=0.05,2); ...
        any(p_DSI<.05,2) & any(p_OSI<.05,2)}}; % DS, OS, DS & OS
    setNames = {{'exc', 'inh&enh', 'inh&sup'}, {'ON','ON+OFF','OFF'}, ...
        {'DS','OS','DS & OS'}};
else % boutons
    neuronSets = {{isSuppr==-1; isSuppr==1}, ... %enh, suppr
        {OnOffRatios>onThr & validRF; ...
        OnOffRatios>=offThr & OnOffRatios<=onThr & validRF; ...
        OnOffRatios<offThr & validRF}, ... % ON, Absolute, OFF
        {any(p_DSI<.05,2) & all(p_OSI>=0.05,2); ...
        any(p_OSI<.05,2) & all(p_DSI>=0.05,2); ...
        any(p_DSI<.05,2) & any(p_OSI<.05,2)}}; % DS, OS, DS & OS
%         {any(p_DSI<.05,2) & all(p_OSI>=0.05,2) & ~isnan(DImaxis); ...
%         any(p_OSI<.05,2) & all(p_DSI>=0.05,2) & ~isnan(DImaxis); ...
%         any(p_DSI<.05,2) & any(p_OSI<.05,2) & ~isnan(DImaxis)}}; % DS, OS, DS & OS
    setNames = {{'enh', 'sup'}, {'ON','ON+OFF','OFF'}, {'DS','OS','DS & OS'}};
end
colors = lines(3);
if strcmp(label, 'boutons')
    combinations = [1 3; 2 1]; % 1. running in darkness; 2. pupil during gratings
else
    combinations = [2 2; 2 1]; % 1. pupil during gray screens; 2. pupil during gratings
end
indRhos = ~any(isnan([rhos{combinations(1,1)}(:,combinations(1,2)), ...
    rhos{combinations(2,1)}(:,combinations(2,2))]),2);


% Plot cumulative distribtutions of correlation values
xLimits = [-0.5 0.5];
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
        test = [];
        for type = 1:length(neuronSets{s})
            ind = neuronSets{s}{type} & indRhos;
            hs = general.plotCumHist(gca, rhos{nonvis}(ind,st), ...
                squeeze(nullRhos{nonvis}(ind,st,:)), [-1 1], colors(type,:));
            plot(mean(rhos{nonvis}(ind,st)), 1.1, 'v', 'MarkerEdgeColor', 'none', ...
                'MarkerFaceColor', colors(type,:))
            h(type) = hs(1);
            lbls{type} = sprintf('%s (n=%d)', setNames{s}{type}, sum(ind));
            test = [test; [rhos{nonvis}(ind,st), ones(sum(ind),1).*type]];
        end
        xlim(xLimits)
        ylim([0 1.2])
        legend(h, lbls, 'Location', 'SouthEast')
        xlabel(sprintf('Correlation with %s', nonvisNames{nonvis}))
        ylabel(sprintf('Proportion of %s', label))
        
        [p,~,stats] = anova1(test(:,1), test(:,2), 'off');
        title(sprintf('Correlation with %s during %s\n(p=%.1e, ANOVA across groups)', ...
            nonvisNames{nonvis}, stimuli{st}, p))
        figure
        stats.gnames = setNames{s}(unique(test(:,2)))';
        multcompare(stats);
        title(sprintf('Correlation with %s during %s\n(p=%.1e, ANOVA across groups)', ...
            nonvisNames{nonvis}, stimuli{st}, p))
    end
end

% Test for significant differences between cell types
xLimits = [-0.5 0.5];
for combi = 1:size(combinations,1)
    nonvis = combinations(combi,1);
    st = combinations(combi,2);
        
    if all(isnan(rhos{nonvis}(:,st)))
        continue
    end
    for s = 1:length(neuronSets)
        types = cell(size(rhos{nonvis},1),1);
        for type = 1:length(neuronSets{s})
            [types{neuronSets{s}{type}}] = deal(setNames{s}{type});
        end
        
        ind = indRhos & ~cellfun(@isempty,types);
        tbl = table(rhos{nonvis}(ind,st), nominal(types(ind)), ...
            dataset(ind), mice(ind), ...
            'VariableNames', {'corr','type','session','mouse'});
        lme = fitlme(tbl, 'corr ~ type + (type|session) + (type|mouse)', 'DummyVarCoding', 'effects')
        lme = fitlme(tbl, 'corr ~ type + (type|session)', 'DummyVarCoding', 'effects')
        lme = fitlme(tbl, 'corr ~ type + (1|session) + (-1 + type|session)', 'DummyVarCoding', 'effects')
        lme = fitlme(tbl, 'corr ~ type + (1|session)', 'DummyVarCoding', 'effects')
        lme = fitlme(tbl, 'corr ~ type + (1|session) + (-1 + type|session) + (type|mouse)', 'DummyVarCoding', 'effects')
        lme = fitlme(tbl, 'corr ~ type + (1|session) + (-1 + type|session) + (1|mouse) + (-1 + type|mouse)', 'DummyVarCoding', 'effects')
        lme = fitlme(tbl, 'corr ~ type + (1|session) + (-1 + type|session) + (1|mouse)', 'DummyVarCoding', 'effects')
        lme = fitlme(tbl, 'corr ~ type + (1|session) + (1|mouse)', 'DummyVarCoding', 'effects')
        lme = fitlme(tbl, 'corr ~ type + (1|session) + (-1 + type|mouse)', 'DummyVarCoding', 'effects')
        anova(lme)
        coefTest(lme, [0 1 -1])
        
        [b, bnames] = fixedEffects(lme);
        bnames = bnames.Variables;
        means = NaN(1, length(neuronSets{s}));
        figure
        hold on
        plot([0 0], [0 1.2], 'k')
        h = zeros(1, length(neuronSets{s}));
        lbls = cell(1, length(neuronSets{s}));
        for type = 1:length(neuronSets{s})
            ind = neuronSets{s}{type} & indRhos;
            hs = general.plotCumHist(gca, rhos{nonvis}(ind,st), ...
                squeeze(nullRhos{nonvis}(ind,st,:)), [-1 1], colors(type,:));
            v = find(strcmp(bnames, ['type_' setNames{s}{type}]));
            if ~isempty(v)
                means(type) = b(1) + b(v);
            else
                means(type) = b(1) - sum(b(2:end));
            end
            plot(means(type), 1.1, 'v', 'MarkerEdgeColor', 'none', ...
                'MarkerFaceColor', colors(type,:))
            h(type) = hs(1);
            lbls{type} = sprintf('%s (n=%d)', setNames{s}{type}, sum(ind));
        end
        xlim(xLimits)
        ylim([0 1.2])
        legend(h, lbls, 'Location', 'SouthEast')
        xlabel(sprintf('Correlation with %s', nonvisNames{nonvis}))
        ylabel(sprintf('Proportion of %s', label))
        title(sprintf('Correlation with %s during %s', nonvisNames{nonvis}, stimuli{st}))
    end
end

% Plot cumulative distribtuion of absolute correlation values
% xLimits = [0 0.85];
% for st = 1:length(stimuli)
%     if all(isnan(rhos(:,st)))
%         continue
%     end
%     for s = 1:length(neuronSets)
%         figure
%         hold on
%         h = zeros(1, length(neuronSets{s}));
%         lbls = cell(1, length(neuronSets{s}));
%         test = [];
%         for type = 1:length(neuronSets{s})
%             ind = neuronSets{s}{type} & ~isnan(rhos(:,st));
%             hs = general.plotCumHist(gca, abs(rhos(ind,st)), ...
%                 squeeze(abs(nullRhos(ind,st,:))), [0 1], colors(type,:));
%             plot(mean(abs(rhos(ind,st))), 1.1, 'v', 'MarkerEdgeColor', 'none', ...
%                 'MarkerFaceColor', colors(type,:))
%             h(type) = hs(1);
%             lbls{type} = sprintf('%s (n=%d)', setNames{s}{type}, sum(ind));
%             test = [test; [abs(rhos(ind,st)), ones(sum(ind),1).*type]];
%         end
%         xlim(xLimits)
%         ylim([0 1.2])
%         legend(h, lbls, 'Location', 'SouthEast')
%         xlabel(sprintf('Absolute correlation with %s', nonvisCorr))
%         ylabel(sprintf('Proportion of %s', label))
%         
%         [p,tbl,stats] = anova1(test(:,1), test(:,2), 'off');
%         title(sprintf('Absolute correlation with %s during %s\n(p=%.1e, ANOVA across groups)', ...
%             nonvisCorr, stimuli{st}, p))
%         figure
%         stats.gnames = setNames{s}';
%         multcompare(stats);
%         title(sprintf('Absolute correlation with %s during %s\n(p=%.1e, ANOVA across groups)', ...
%             nonvisCorr, stimuli{st}, p))
%     end
% end

%% Plot cumulative distributions of tuning modulations for cell types
if strcmp(label, 'neurons')
    neuronSets = {{isGad==-1; isGad==1&isSuppr==-1; isGad==1&isSuppr==1}, ... % exc, inh&enh, inh&suppr
        {OnOffRatios>onThr & validRF; ...
        OnOffRatios>=offThr & OnOffRatios<=onThr & validRF; ...
        OnOffRatios<offThr & validRF}, ... % ON, Absolute, OFF
        {any(p_DSI<.05,2) & all(p_OSI>=0.05,2); ...
        any(p_OSI<.05,2) & all(p_DSI>=0.05,2); ...
        any(p_DSI<.05,2) & any(p_OSI<.05,2)}}; % DS, OS, DS & OS
    setNames = {{'exc', 'inh&enh', 'inh&sup'}, {'ON','ON+OFF','OFF'}, ...
        {'DS','OS','DS & OS'}};
else % boutons
    neuronSets = {{isSuppr==-1; isSuppr==1}, ... %enh, suppr
        {OnOffRatios>onThr & validRF; ...
        OnOffRatios>=offThr & OnOffRatios<=onThr & validRF; ...
        OnOffRatios<offThr & validRF}, ... % ON, Absolute, OFF
        {any(p_DSI<.05,2) & all(p_OSI>=0.05,2); ...
        any(p_OSI<.05,2) & all(p_DSI>=0.05,2); ...
        any(p_DSI<.05,2) & any(p_OSI<.05,2)}}; % DS, OS, DS & OS
    setNames = {{'enh', 'sup'}, {'ON','ON+OFF','OFF'}, {'DS','OS','DS & OS'}};
end
colors = lines(3);

DIs = {DImaxisAll, DImaxis, DImods};
nullDIs = {nullDImaxisAll, nullDImaxis, nullDImods};
modNames = {'resp @ pref orientation (incl untuned)', ...
    'resp @ pref orientation (tuned only)', 'tuning depth'};

% Plot cumulative distribtutions of tuning modulations
xLimits = [-0.3 0.4; -0.3 0.4; -0.7 0.4];
% xLimits = [0.65 0.75];
for m = 1:length(DIs)
    for s = 1:length(neuronSets)
        figure
        hold on
        plot([0 0], [0 1.2], 'k')
        h = zeros(1, length(neuronSets{s}));
        lbls = cell(1, length(neuronSets{s}));
        test = [];
        for type = 1:length(neuronSets{s})
            ind = neuronSets{s}{type} & ~isnan(DIs{m});
            hs = general.plotCumHist(gca, DIs{m}(ind), ...
                nullDIs{m}(ind,:), [-1 1], colors(type,:));
            plot(mean(DIs{m}(ind)), 1.1, 'v', 'MarkerEdgeColor', 'none', ...
                'MarkerFaceColor', colors(type,:))
            h(type) = hs(1);
            lbls{type} = sprintf('%s (n=%d)', setNames{s}{type}, sum(ind));
            test = [test; [DIs{m}(ind), ones(sum(ind),1).*type]];
        end
        xlim(xLimits(m,:))
        ylim([0 1.2])
        legend(h, lbls, 'Location', 'SouthEast')
        xlabel(sprintf('DI of %s', modNames{m}))
        ylabel(sprintf('Proportion of %s', label))
        
        [p,tbl,stats] = anova1(test(:,1), test(:,2), 'off');
        title(sprintf('Tuning modulation in %s\n(p=%.1e, ANOVA across groups)', ...
            modNames{m}, p))
        figure
        stats.gnames = setNames{s}';
        multcompare(stats);
        title(sprintf('Tuning modulation in %s\n(p=%.1e, ANOVA across groups)', ...
            modNames{m}, p))
    end
end

% Test for significant differences between cell types
xLimits = [-0.3 0.3; -0.3 0.3; -0.4 0.4];
for m = [1 3]
    for s = 1:length(neuronSets)
        types = cell(size(rhos{nonvis},1),1);
        for type = 1:length(neuronSets{s})
            [types{neuronSets{s}{type}}] = deal(setNames{s}{type});
        end
        
        ind = ~isnan(DIs{m}) & ~cellfun(@isempty,types);
        tbl = table(DIs{m}(ind), nominal(types(ind)), dataset(ind), mice(ind), ...
            'VariableNames', {'AI','type','session','mouse'});
        lme = fitlme(tbl, 'AI ~ type + (type|session) + (type|mouse)', 'DummyVarCoding', 'effects')
        lme = fitlme(tbl, 'AI ~ type + (type|session)', 'DummyVarCoding', 'effects')
        lme = fitlme(tbl, 'AI ~ type + (1|session) + (-1 + type|session)', 'DummyVarCoding', 'effects')
        lme = fitlme(tbl, 'AI ~ type + (1|session)', 'DummyVarCoding', 'effects')
        lme = fitlme(tbl, 'AI ~ type + (1|session) + (-1 + type|session) + (type|mouse)', 'DummyVarCoding', 'effects')
        lme = fitlme(tbl, 'AI ~ type + (1|session) + (-1 + type|session) + (1|mouse) + (-1 + type|mouse)', 'DummyVarCoding', 'effects')
        lme = fitlme(tbl, 'AI ~ type + (1|session) + (-1 + type|session) + (1|mouse)', 'DummyVarCoding', 'effects')
        lme = fitlme(tbl, 'AI ~ type + (-1 + type|session) + (1|mouse) + (-1 + type|mouse)', 'DummyVarCoding', 'effects')
        lme = fitlme(tbl, 'AI ~ type + (1|session) + (1|mouse) + (-1 + type|mouse)', 'DummyVarCoding', 'effects')
        lme = fitlme(tbl, 'AI ~ type + (type|session) + (1|mouse)', 'DummyVarCoding', 'effects')
        lme = fitlme(tbl, 'AI ~ type + (type|session) + (1|mouse) + (-1 + type|mouse)', 'DummyVarCoding', 'effects')
        lme = fitlme(tbl, 'AI ~ type + (type|session) + (-1 + type|mouse)', 'DummyVarCoding', 'effects')
        anova(lme)
        coefTest(lme, [0 1 -1])
        
        [b, bnames] = fixedEffects(lme);
        bnames = bnames.Variables;
        means = NaN(1, length(neuronSets{s}));
        figure
        hold on
        plot([0 0], [0 1.2], 'k')
        h = zeros(1, length(neuronSets{s}));
        lbls = cell(1, length(neuronSets{s}));
        for type = 1:length(neuronSets{s})
            ind = ~isnan(DIs{m}) & neuronSets{s}{type};
            hs = general.plotCumHist(gca, DIs{m}(ind), ...
                nullDIs{m}(ind,:), [-1 1], colors(type,:));
            v = find(strcmp(bnames, ['type_' setNames{s}{type}]));
            if ~isempty(v)
                means(type) = b(1) + b(v);
            else
                means(type) = b(1) - sum(b(2:end));
            end
            plot(means(type), 1.1, 'v', 'MarkerEdgeColor', 'none', ...
                'MarkerFaceColor', colors(type,:))
            h(type) = hs(1);
            lbls{type} = sprintf('%s (n=%d)', setNames{s}{type}, sum(ind));
        end
        xlim(xLimits(m,:))
        ylim([0 1.2])
        set(gca, 'XTick', sort([xLimits(m,:) 0]))
        legend(h, lbls, 'Location', 'SouthEast')
        xlabel(sprintf('DI of %s', modNames{m}))
        ylabel(sprintf('Proportion of %s', label))
        title(sprintf('Tuning modulation in %s', modNames{m}))
    end
end

% % Plot cumulative distribtutions of absolute tuning modulations
% DIs = {abs(DImaxis), abs(DImaxisAll), abs(DImods)};
% nullDIs = {abs(nullDImaxis), abs(nullDImaxisAll), abs(nullDImods)};
% modNames = {'resp @ pref orientation (tuned only)', 'resp @ pref orientation (incl untuned)', 'tuning depth'};
% xLimits = [0.65 0.65 0.75];
% for m = 1:length(DIs)
%     for s = 1:length(neuronSets)
%         figure
%         hold on
%         h = zeros(1, length(neuronSets{s}));
%         lbls = cell(1, length(neuronSets{s}));
%         test = [];
%         for type = 1:length(neuronSets{s})
%             ind = neuronSets{s}{type} & ~isnan(DIs{m});
%             hs = general.plotCumHist(gca, DIs{m}(ind), ...
%                 nullDIs{m}(ind,:), [0 1], colors(type,:));
%             plot(mean(DIs{m}(ind)), 1.1, 'v', 'MarkerEdgeColor', 'none', ...
%                 'MarkerFaceColor', colors(type,:))
%             h(type) = hs(1);
%             lbls{type} = sprintf('%s (n=%d)', setNames{s}{type}, sum(ind));
%             test = [test; [DIs{m}(ind), ones(sum(ind),1).*type]];
%         end
%         xlim([0 xLimits(m)])
%         ylim([0 1.2])
%         legend(h, lbls, 'Location', 'SouthEast')
%         xlabel(sprintf('Absolute DI of %s', modNames{m}))
%         ylabel(sprintf('Proportion of %s', label))
%         
%         [p,tbl,stats] = anova1(test(:,1), test(:,2), 'off');
%         title(sprintf('Absolute tuning modulation in %s\n(p=%.1e, ANOVA across groups)', ...
%             modNames{m}, p))
%         figure
%         stats.gnames = setNames{s}';
%         multcompare(stats);
%         title(sprintf('Absolute tuning modulation in %s\n(p=%.1e, ANOVA across groups)', ...
%             modNames{m}, p))
%     end
% end

%% Plot correlations and tuning modulation against preferred orientation
binSize = 90;
prefDirs2 = prefDirs;
ind = prefDirs > 360-binSize/2;
prefDirs2(ind) = prefDirs2(ind)-360;

edges = -binSize/2 : binSize : 360+binSize/2;
[N,~,bin] = histcounts(prefDirs2, edges);
bins = 0 : binSize : 270;

ind = ~isnan(DImaxis);
tbl = table(DImaxis(ind), nominal(bin(ind)), dataset(ind), mice(ind), ...
    'VariableNames', {'DI','prefDir','session','mouse'});
lme = fitlme(tbl, 'DI ~ prefDir + (prefDir|session) + (prefDir|mouse)', ...
    'DummyVarCoding', 'effects')
lme = fitlme(tbl, 'DI ~ prefDir + (prefDir|session)', 'DummyVarCoding', 'effects')
lme = fitlme(tbl, 'DI ~ prefDir + (1|session)', 'DummyVarCoding', 'effects')

means = NaN(1, length(bins));
ses = NaN(1, length(bins));
for b = 1:length(bins)
    means(b) = nanmean(DImaxis(bin==b));
    ses(b) = nanstd(DImaxis(bin==b)) ./ sqrt(sum(~isnan(DImaxis(bin==b))));
end

yLim = [-0.45 0.45];
figure
hold on
di = DImaxis;
di(di<yLim(1)) = yLim(1);
di(di>yLim(2)) = yLim(2);
h = scatter(prefDirs2, di, 15, 'k', 'filled');
errorbar(bins, means, ses, 'r.', 'CapSize', 0, 'LineWidth', 1, 'MarkerSize', 20)
xlim([-binSize/2 360-binSize/2])
ylim(yLim)
set(gca, 'XTick', 0:90:270, 'YTick', [yLim(1) 0 yLim(2)])
xlabel('Preferred direction')
ylabel('DI @ pref. resp.')
title(sprintf('%s (n=%d)', label, sum(all(~isnan([prefDirs,DImaxis]),2))))

%% Plot percentage of large pupil trials for each presented direction
[p,~,stats] = anova1(behCondStim, directionStim);

bins = unique(directionStim);
means = NaN(1, length(bins));
ses = NaN(1, length(bins));
for b = 1:length(bins)
    means(b) = nanmean(behCondStim(directionStim==bins(b)));
    ses(b) = nanstd(behCondStim(directionStim==bins(b))) ./ sqrt(sum(directionStim==bins(b)));
end
means(end+1) = means(1);
ses(end+1) = ses(1);
bins(end+1) = 360;
means = means - 1;
figure
hold on
plot(bins, means, 'k', 'LineWidth', 2)
plot(bins, means+ses, 'k--')
plot(bins, means-ses, 'k--')
ylim([0.3 .7])
xlim([0 360])
xlabel('Direction of stimulus')
ylabel('Proportion of large pupil trials')
title(sprintf('For %s: distribution of large pupil trials', label))
set(gca, 'XTick', 0:90:360)

%% Plot correlations with pupil/running against tuning modulation for cell types
rhos = {rhosRun, rhosPupil};
rho_pValues = {pValues_rhoRun, pValues_rhoPupil};
nonvisNames = {'running','pupil'};
if strcmp(label, 'neurons')
    neuronSets = {{isGad==-1; isGad==1&isSuppr==-1; isGad==1&isSuppr==1}, ... % exc, inh&enh, inh&suppr
        {OnOffRatios>0.7; OnOffRatios>=0.3&OnOffRatios<=0.7; OnOffRatios<0.3}, ... % ON, Absolute, OFF
        {OnOffRatios>0.7&isSuppr==-1; OnOffRatios>=0.3&OnOffRatios<=0.7&isSuppr==-1; ...
        OnOffRatios<0.3&isSuppr==-1}, ... % ON&enh, Absolute&enh, OFF&enh
        {OnOffRatios>0.7&isSuppr==1; OnOffRatios>=0.3&OnOffRatios<=0.7&isSuppr==1; ...
        OnOffRatios<0.3&isSuppr==1}}; % ON&supp, Absolute&supp, OFF&supp
    setNames = {{'exc', 'inh&enh', 'inh&sup'}, {'ON','ON/OFF','OFF'}, ...
        {'ON&enh','Abs&enh','OFF&enh'}, {'ON&sup','Abs&sup','OFF&sup'}};
else % boutons
    neuronSets = {{isSuppr==-1; isSuppr==1}, ... %enh, suppr
        {OnOffRatios>onThr & validRF; ...
        OnOffRatios>=offThr & OnOffRatios<=onThr & validRF; ...
        OnOffRatios<offThr & validRF}, ... % ON, ON+OFF, OFF
        {any(p_DSI<.05,2) & all(p_OSI>=0.05,2) & ~isnan(DImaxis); ...
        any(p_OSI<.05,2) & all(p_DSI>=0.05,2) & ~isnan(DImaxis); ...
        any(p_DSI<.05,2) & any(p_OSI<.05,2) & ~isnan(DImaxis)}}; % DS, OS, DS & OS
    setNames = {{'enh', 'sup'}, {'ON','ON+OFF','OFF'}, {'DS','OS','DS & OS'}};
end

stims = [3 1]; % use responses during darkness for correlations with running
               % use responses during gratings for correlations with pupil
markers = {'o','v','<','d'};
sizes = [5 20 20 30];
for nonvis = 1:2
    st = stims(nonvis);
    if all(isnan(rhos{nonvis}(:,st)))
        continue
    end
    signif1 = [rho_pValues{nonvis}(:,st)>=0.05 & pValues_DImaxisAll>=0.05, ...
        rho_pValues{nonvis}(:,st)<0.05 & pValues_DImaxisAll>=0.05, ...
        rho_pValues{nonvis}(:,st)>=0.05 & pValues_DImaxisAll<0.05, ...
        rho_pValues{nonvis}(:,st)<0.05 & pValues_DImaxisAll<0.05];
    signif2 = [rho_pValues{nonvis}(:,st)>=0.05 & pValues_DImods>=0.05, ...
        rho_pValues{nonvis}(:,st)<0.05 & pValues_DImods>=0.05, ...
        rho_pValues{nonvis}(:,st)>=0.05 & pValues_DImods<0.05, ...
        rho_pValues{nonvis}(:,st)<0.05 & pValues_DImods<0.05];
    for s = 1:length(neuronSets)
        xmn = Inf;
        xmx = -Inf;
        y1mn = Inf;
        y1mx = -Inf;
        y2mn = Inf;
        y2mx = -Inf;
        for type = 1:length(neuronSets{s})
            xmn = min([xmn; rhos{nonvis}(neuronSets{s}{type}, st)]);
            xmx = max([xmx; rhos{nonvis}(neuronSets{s}{type}, st)]);
            y1mn = min([y1mn; DImaxisAll(neuronSets{s}{type})]);
            y1mx = max([y1mx; DImaxisAll(neuronSets{s}{type})]);
            y2mn = min([y2mn; DImods(neuronSets{s}{type})]);
            y2mx = max([y2mx; DImods(neuronSets{s}{type})]);
        end
        for type = 1:length(neuronSets{s})
            ind = neuronSets{s}{type};
            
            figure
            hold on
            for j = 1:4
                h = scatter(rhos{nonvis}(ind & signif1(:,j), st), ...
                    DImaxisAll(ind & signif1(:,j)), sizes(j), 'k', 'filled', ...
                    markers{j}, 'MarkerFaceAlpha', 0.3);
                row = dataTipTextRow('Subject', ...
                    {corrsPupil(dataset(ind & signif1(:,j))).subject}');
                h.DataTipTemplate.DataTipRows(end+1) = row;
                row = dataTipTextRow('Date', ...
                    {corrsPupil(dataset(ind & signif1(:,j))).date}');
                h.DataTipTemplate.DataTipRows(end+1) = row;
                row = dataTipTextRow('Plane',plane(ind & signif1(:,j)));
                h.DataTipTemplate.DataTipRows(end+1) = row;
                row = dataTipTextRow('Neuron',neuron(ind & signif1(:,j)));
                h.DataTipTemplate.DataTipRows(end+1) = row;
                h.DataTipTemplate.Interpreter = 'none';
            end
            plot([0 0], [y1mn y1mx], 'k:')
            plot([xmn xmx], [0 0], 'k:')
            b = regress(DImaxisAll(ind), [ones(sum(ind),1), rhos{nonvis}(ind, st)]);
            plot([xmn xmx], b(1) + b(2).*[xmn xmx], 'r')
            axis([xmn xmx y1mn y1mx])
            axis square
            ind1 = ind & ~isnan(rhos{nonvis}(:,st)) & ~isnan(DImaxisAll);
            [rho, p] = corr(rhos{nonvis}(ind1, st), DImaxisAll(ind1));
            xlabel(sprintf('Corr. with %s during %s', nonvisNames{nonvis}, ...
                stimuli{st}))
            ylabel('DI of resp @ pref direction')
            title(sprintf('%s (n = %d, rho = %.3f, p = %.3f)', setNames{s}{type}, ...
                sum(ind1), rho, p))
            
            figure
            hold on
            for j = 1:4
                h = scatter(rhos{nonvis}(ind & signif2(:,j), st), ...
                    DImods(ind & signif2(:,j)), sizes(j), 'k', 'filled', ...
                    markers{j}, 'MarkerFaceAlpha', 0.3);
                row = dataTipTextRow('Subject', ...
                    {corrsPupil(dataset(ind & signif2(:,j))).subject}');
                h.DataTipTemplate.DataTipRows(end+1) = row;
                row = dataTipTextRow('Date', ...
                    {corrsPupil(dataset(ind & signif2(:,j))).date}');
                h.DataTipTemplate.DataTipRows(end+1) = row;
                row = dataTipTextRow('Plane',plane(ind & signif2(:,j)));
                h.DataTipTemplate.DataTipRows(end+1) = row;
                row = dataTipTextRow('Neuron',neuron(ind & signif2(:,j)));
                h.DataTipTemplate.DataTipRows(end+1) = row;
                h.DataTipTemplate.Interpreter = 'none';
            end
            plot([0 0], [y2mn y2mx], 'k:')
            plot([xmn xmx], [0 0], 'k:')
            b = regress(DImods(ind), [ones(sum(ind),1), rhos{nonvis}(ind, st)]);
            plot([xmn xmx], b(1) + b(2).*[xmn xmx], 'r')
            axis([xmn xmx y2mn y2mx])
            axis square
            ind1 = ind & ~isnan(rhos{nonvis}(:,st)) & ~isnan(DImods);
            [rho, p] = corr(rhos{nonvis}(ind1, st), DImods(ind1));
            xlabel(sprintf('Corr. with %s during %s', nonvisNames{nonvis}, ...
                stimuli{st}))
            ylabel('DI of tuning depth')
            title(sprintf('%s (n = %d, rho = %.3f, p = %.3f)', setNames{s}{type}, ...
                sum(ind1), rho, p))
        end
    end
end

%% Plot correlation with pupil against correlation with running for cell types
if strcmp(label, 'neurons')
    neuronSets = {{isGad==-1; isGad==1&isSuppr==-1; isGad==1&isSuppr==1}, ... % exc, inh&enh, inh&suppr
        {OnOffRatios>0.7; OnOffRatios>=0.3&OnOffRatios<=0.7; OnOffRatios<0.3}, ... % ON, Absolute, OFF
        {OnOffRatios>0.7&isSuppr==-1; OnOffRatios>=0.3&OnOffRatios<=0.7&isSuppr==-1; ...
        OnOffRatios<0.3&isSuppr==-1}, ... % ON&enh, Absolute&enh, OFF&enh
        {OnOffRatios>0.7&isSuppr==1; OnOffRatios>=0.3&OnOffRatios<=0.7&isSuppr==1; ...
        OnOffRatios<0.3&isSuppr==1}}; % ON&supp, Absolute&supp, OFF&supp
    setNames = {{'exc', 'inh&enh', 'inh&sup'}, {'ON','ON/OFF','OFF'}, ...
        {'ON&enh','Abs&enh','OFF&enh'}, {'ON&sup','Abs&sup','OFF&sup'}};
else % boutons
    neuronSets = {{isSuppr==-1; isSuppr==1}, ... %enh, suppr
        {OnOffRatios>onThr & validRF; ...
        OnOffRatios>=offThr & OnOffRatios<=onThr & validRF; ...
        OnOffRatios<offThr & validRF}, ... % ON, ON+OFF, OFF
        {any(p_DSI<.05,2) & all(p_OSI>=0.05,2) & ~isnan(DImaxis); ...
        any(p_OSI<.05,2) & all(p_DSI>=0.05,2) & ~isnan(DImaxis); ...
        any(p_DSI<.05,2) & any(p_OSI<.05,2) & ~isnan(DImaxis)}}; % DS, OS, DS & OS
    setNames = {{'enh', 'sup'}, {'ON','ON+OFF','OFF'}, {'DS','OS','DS & OS'}};
end

stims = [3 1]; % use responses during darkness for correlations with running
               % use responses during gratings for correlations with pupil
markers = {'o','v','<','d'};
sizes = [5 20 20 30];
signif = [pValues_rhoRun(:,stims(1))>=0.05 & pValues_rhoPupil(:,stims(2))>=0.05, ...
    pValues_rhoRun(:,stims(1))<0.05 & pValues_rhoPupil(:,stims(2))>=0.05, ...
    pValues_rhoRun(:,stims(1))>=0.05 & pValues_rhoPupil(:,stims(2))<0.05, ...
    pValues_rhoRun(:,stims(1))<0.05 & pValues_rhoPupil(:,stims(2))<0.05];
for s = 1:length(neuronSets)
    mn = Inf;
    mx = -Inf;
    for type = 1:length(neuronSets{s})
        mn = min([mn; rhosRun(neuronSets{s}{type}, stims(1)); ...
            rhosPupil(neuronSets{s}{type}, stims(2))]);
        mx = max([mx; rhosRun(neuronSets{s}{type}, stims(1)); ...
            rhosPupil(neuronSets{s}{type}, stims(2))]);
    end
    for type = 1:length(neuronSets{s})
        ind = neuronSets{s}{type} & ~any(isnan([rhosRun(:,stims(1)), ...
            rhosPupil(:,stims(2))]),2);
        
        figure
        hold on
        for j = 1:4
            h = scatter(rhosRun(ind & signif(:,j), stims(1)), ...
                rhosPupil(ind & signif(:,j), stims(2)), sizes(j), 'k', 'filled', ...
                markers{j}, 'MarkerFaceAlpha', 0.3);
            row = dataTipTextRow('Subject', ...
                {corrsPupil(dataset(ind & signif(:,j))).subject}');
            h.DataTipTemplate.DataTipRows(end+1) = row;
            row = dataTipTextRow('Date', ...
                {corrsPupil(dataset(ind & signif(:,j))).date}');
            h.DataTipTemplate.DataTipRows(end+1) = row;
            row = dataTipTextRow('Plane',plane(ind & signif(:,j)));
            h.DataTipTemplate.DataTipRows(end+1) = row;
            row = dataTipTextRow('Neuron',neuron(ind & signif(:,j)));
            h.DataTipTemplate.DataTipRows(end+1) = row;
            h.DataTipTemplate.Interpreter = 'none';
        end
        plot([0 0], [mn mx], 'k:')
        plot([mn mx], [0 0], 'k:')
        plot([mn mx], [mn mx], 'k')
        b = regress(rhosPupil(ind, stims(2)), [ones(sum(ind),1), rhosRun(ind, stims(1))]);
        plot([mn mx], b(1) + b(2).*[mn mx], 'r')
        axis([mn mx mn mx])
        axis square
        [rho, p] = corr(rhosRun(ind, stims(1)), rhosPupil(ind, stims(2)));
        xlabel(sprintf('Corr. with running during %s', stimuli{stims(1)}))
        ylabel(sprintf('Corr. with pupil during %s', stimuli{stims(2)}))
        title(sprintf('%s (n = %d, rho = %.3f, p = %.3f)', setNames{s}{type}, ...
            sum(ind), rho, p))
    end
end

%% Plot DI of minimum response versus DI of maximum response for cell types
if strcmp(label, 'neurons')
    neuronSets = {{isGad==-1; isGad==1&isSuppr==-1; isGad==1&isSuppr==1}, ... % exc, inh&enh, inh&suppr
        {OnOffRatios>0.7; OnOffRatios>=0.3&OnOffRatios<=0.7; OnOffRatios<0.3}, ... % ON, Absolute, OFF
        {OnOffRatios>0.7&isSuppr==-1; OnOffRatios>=0.3&OnOffRatios<=0.7&isSuppr==-1; ...
        OnOffRatios<0.3&isSuppr==-1}, ... % ON&enh, Absolute&enh, OFF&enh
        {OnOffRatios>0.7&isSuppr==1; OnOffRatios>=0.3&OnOffRatios<=0.7&isSuppr==1; ...
        OnOffRatios<0.3&isSuppr==1}}; % ON&supp, Absolute&supp, OFF&supp
    setNames = {{'exc', 'inh&enh', 'inh&sup'}, {'ON','ON/OFF','OFF'}, ...
        {'ON&enh','Abs&enh','OFF&enh'}, {'ON&sup','Abs&sup','OFF&sup'}};
else % boutons
    neuronSets = {{isSuppr==-1; isSuppr==1}, ... %enh, suppr
        {OnOffRatios>onThr & validRF; ...
        OnOffRatios>=offThr & OnOffRatios<=onThr & validRF; ...
        OnOffRatios<offThr & validRF}, ... % ON, Absolute, OFF
        {any(p_DSI<.05,2) & all(p_OSI>=0.05,2) & ~isnan(DImaxis); ...
        any(p_OSI<.05,2) & all(p_DSI>=0.05,2) & ~isnan(DImaxis); ...
        any(p_DSI<.05,2) & any(p_OSI<.05,2) & ~isnan(DImaxis)}}; % DS, OS, DS & OS
    setNames = {{'enh', 'sup'}, {'ON','ON+OFF','OFF'}, {'DS','OS','DS & OS'}};
end

for s = 1:length(neuronSets)
    mn = Inf;
    mx = -Inf;
    for type = 1:length(neuronSets{s})
        mn = min([mn; DIminis(neuronSets{s}{type}); ...
            DImaxis(neuronSets{s}{type})]);
        mx = max([mx; DIminis(neuronSets{s}{type}); ...
            DImaxis(neuronSets{s}{type})]);
    end
    for type = 1:length(neuronSets{s})
        ind = neuronSets{s}{type} & ~isnan(DIminis) & ~isnan(DImaxis);
        
        figure
        hold on
        h = scatter(DIminis(ind), DImaxis(ind), [], 'k', 'filled', ...
            'MarkerFaceAlpha', 0.3);
        row = dataTipTextRow('Subject', ...
            {corrsPupil(dataset(ind)).subject}');
        h.DataTipTemplate.DataTipRows(end+1) = row;
        row = dataTipTextRow('Date', ...
            {corrsPupil(dataset(ind)).date}');
        h.DataTipTemplate.DataTipRows(end+1) = row;
        row = dataTipTextRow('Plane',plane(ind));
        h.DataTipTemplate.DataTipRows(end+1) = row;
        row = dataTipTextRow('Neuron',neuron(ind));
        h.DataTipTemplate.DataTipRows(end+1) = row;
        h.DataTipTemplate.Interpreter = 'none';
        
        plot([0 0], [mn mx], 'k:')
        plot([mn mx], [0 0], 'k:')
        b = regress(DImaxis(ind), [ones(sum(ind),1), DIminis(ind)]);
        plot([mn mx], b(1) + b(2).*[mn mx], 'r')
        axis([mn mx mn mx])
        axis square
        [rho, p] = corr(DIminis(ind), DImaxis(ind));
        xlabel('DI of resp @ null direction')
        ylabel('DI of resp @ pref direction')
        title(sprintf('%s (n = %d, rho = %.3f, p = %.3f)', setNames{s}{type}, ...
            sum(ind), rho, p))
    end
end

%% Select example boutons of different RF types
inds = find(validRF);
[validRatios,validOrder]=sort(OnOffRatios(inds));
validInfo = [validRatios, ...
    dataset(inds(validOrder)), plane(inds(validOrder)), neuron(inds(validOrder))];


onInds = find(OnOffRatios > onThr & validRF);
[onRhos,onOrder] = sort(rhosPupil(onInds,strcmp(stimuli,'gratings')));
onInfo = [OnOffRatios(onInds(onOrder)), ...
    onRhos, pValues_rhoPupil(onInds(onOrder),strcmp(stimuli,'gratings')), ...
    rhosPupil(onInds(onOrder),strcmp(stimuli,'grayScreen')), ...
    pValues_rhoPupil(onInds(onOrder),strcmp(stimuli,'grayScreen')), ...
    DImaxisAll(onInds(onOrder)), pValues_DImaxisAll(onInds(onOrder)), ...
    dataset(onInds(onOrder)), plane(onInds(onOrder)), neuron(onInds(onOrder))];
offInds = find(OnOffRatios < offThr & validRF);
[offRhos,offOrder] = sort(rhosPupil(offInds,strcmp(stimuli,'gratings')));
offInfo = [OnOffRatios(offInds(offOrder)), ...
    offRhos, pValues_rhoPupil(offInds(offOrder),strcmp(stimuli,'gratings')), ...
    rhosPupil(offInds(offOrder),strcmp(stimuli,'grayScreen')), ...
    pValues_rhoPupil(offInds(offOrder),strcmp(stimuli,'grayScreen')), ...
    DImaxisAll(offInds(offOrder)), pValues_DImaxisAll(offInds(offOrder)), ...
    dataset(offInds(offOrder)), plane(offInds(offOrder)), neuron(offInds(offOrder))];
onOffInds = find(OnOffRatios >= offThr & OnOffRatios <= onThr & validRF);
[onOffRhos,onOffOrder] = sort(rhosPupil(onOffInds,strcmp(stimuli,'gratings')));
onOffInfo = [OnOffRatios(onOffInds(onOffOrder)), ...
    onOffRhos, pValues_rhoPupil(onOffInds(onOffOrder),strcmp(stimuli,'gratings')), ...
    rhosPupil(onOffInds(onOffOrder),strcmp(stimuli,'grayScreen')), ...
    pValues_rhoPupil(onOffInds(onOffOrder),strcmp(stimuli,'grayScreen')), ...
    DImaxisAll(onOffInds(onOffOrder)), pValues_DImaxisAll(onOffInds(onOffOrder)), ...
    dataset(onOffInds(onOffOrder)), plane(onOffInds(onOffOrder)), neuron(onOffInds(onOffOrder))];

%% Compare number of tuned vs untuned units for suppressed and driven units
isTuned = ~any(isnan(maxis),2);
N = [isSuppr==-1 & isTuned, isSuppr==-1 & ~isTuned, isSuppr==1 & isTuned, ...
    isSuppr==1 & ~isTuned];
figure
bar(sum(N))
set(gca, 'XTickLabel', {'driv+ tun+','driv+ tun-','driv- tun+','driv- tun-'})
fprintf('%.2f%% driven neurons are tuned\n', sum(N(:,1))/sum(sum(N(:,1:2)))*100)
fprintf('%.2f%% suppressed neurons are not tuned\n', sum(N(:,4))/sum(sum(N(:,3:4)))*100)

% determine mean values while accounting for sessions as random effects
ind = ~isnan(isSuppr);
tbl = table(nominal(isSuppr(ind)), double(isTuned(ind)), dataset(ind), ...
    'VariableNames', {'suppr','tuned','session'});
lme = fitlme(tbl, 'tuned ~ suppr + (suppr|session)', 'DummyVarCoding', 'effects');

%% OLD: Collect relevant variables from visual noise data
file = 'receptiveFields.mat';
data = load(fullfile(folderRFResults, file));
RFs = data.RFs;
file = 'receptiveFields_stim+run_4RFmodels.mat';
data = load(fullfile(folderRFResults, file));
RFs_old = data.RFs;

evStim = NaN(0,4);
confIntStimOnly = NaN(0,4,2);
lambdasStim = NaN(0,4);
sizeRun = NaN(0,1);
signRF = NaN(0,4);
onOffPeaks = NaN(0,2);
for k = 1:length(corrsRun)
%     disp(size(evStim,1))
    krf = find(strcmp(corrsRun(k).subject, {RFs.subject}) & ...
        strcmp(corrsRun(k).date, {RFs.date}));
    if isempty(krf) || isempty(RFs(krf).RFTimes) || sum(dataset == k)==0
        numCells = sum(dataset == k);
        evStim(end+(1:numCells),:) = NaN;
        confIntStimOnly(end+(1:numCells),:,:) = NaN;
        lambdasStim(end+(1:numCells),:) = NaN;
        sizeRun(end+(1:numCells),1) = NaN;
        signRF(end+(1:numCells),1) = NaN;
        onOffPeaks(end+(1:numCells),:) = NaN;
        continue
    end
    runBin = diff(RFs(krf).runWindow(1:2));
    for iPlane = 1:length(RFs(krf).planes)
        evStim = [evStim; RFs(krf).plane(iPlane).explainedVariances_stimOnly];
        confIntStimOnly = [confIntStimOnly; ...
            RFs(krf).plane(iPlane).explainedVariances_stimOnly_confInt];
        lambdasStim = [lambdasStim; RFs(krf).plane(iPlane).lambdas];
        
        runKs = RFs(krf).plane(iPlane).runningKernels; %[t x neuron]
        sizeRun = [sizeRun; mean(runKs,1)'./runBin];
        
        rfs = RFs(krf).plane(iPlane).receptiveFields; %[row x column x t x neuron x model]
        rfs = reshape(rfs, [], size(rfs,4), size(rfs,5)); %[pixels x neuron x model]
        [~,maxPx] = max(abs(rfs), [], 1);
        linindsLin = sub2ind(size(rfs), maxPx(:), ...
            repmat((1:size(rfs,2))',size(rfs,3),1), ...
            reshape(repmat(1:size(rfs,3),size(rfs,2),1),[],1));
        signRF = [signRF; reshape(sign(rfs(linindsLin)), size(rfs,2), size(rfs,3))];
    end
end
% switch sign of "black" stimulus model because negative value means the
% neuron is driven by black
signRF(:,end) = -signRF(:,end);

% s explained variance to 0 if not significant, or smaller 0.015, or if
% lambda >= 1
ind = evStim < confIntStimOnly(:,:,2) | evStim < 0.015 | ...
    lambdasStim >= 1;
trueEVstim = evStim;
trueEVstim(ind) = 0;

% get explained variance and sign of RF for best stimulus model (for each
% neuron)
[evStim_bestModel, bestModel] = max(trueEVstim, [], 2);
ind = sub2ind(size(trueEVstim), (1:size(trueEVstim,1))', bestModel);
signStim_bestModel = signRF(ind);
ind = evStim_bestModel == 0;
bestModel(ind) = NaN;
signStim_bestModel(ind) = NaN;

% get ON-OFF-ratio (if 1, completely ON, if 0, completely OFF)
OnOffRatios = NaN(length(dataset),1);
OnOffRatios(bestModel == 3) = 1;
OnOffRatios(bestModel == 4) = 0;
ind = bestModel == 2;
OnOffRatios(ind) = evStim(ind,3) ./ (evStim(ind,3) + evStim(ind,4));