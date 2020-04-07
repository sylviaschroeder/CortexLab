%% Define data

label = 'neurons';
% label = 'boutons';

if strcmp(label, 'neurons')
    % % Neurons in sSC
    examples = [];
    % OLD:
%     exampleSet = 7; % SS047, 2015-11-23, 1
%     examples = [1 168; 1 207]; %; 2 71];
    % exampleSet = 10; % SS048, 2015-12-02, 1
    % examples = [1 105; 2 6; 1 37; 1 44]; % [plane, cellID]
else
    % Boutons in sSC
    % (1) SS076, 2017-09-28 (OS), (2) SS077, 2017-10-05 (DS), (3) SS069,
    % 2016-10-13 (suppressed), (4) SS077, 2017-10-03 (driven, pos corr
    % during gratings, neg DI for resp @ pref dir)
    examples = [12 1 137; 15 1 24; 3 3 4; 14 1 56];
    % OLD:
%     exampleSet = 18; % SS078, 2017-10-05, 2
%     examples = [1 125; 1 87; 1 15];
    % OLD:
%     exampleSet = 4; % SS069, 2016-10-21, 1
%     examples = [1 5; 1 120; 1 164];
end

stimuli = {'gratings', 'grayScreen', 'dark'};

%% Folders
% file location depending on PCs
folderPC = 'C:\STORAGE\OneDrive - University College London';
% folderPC = 'C:\Users\skgtchr\OneDrive - University College London'; % ZEN notebook
% folderPC = 'C:\Storage\OneDrive - University College London'; % Lenovo laptop
% data folders
folderROIData = fullfile(folderPC, 'Lab\DATA\InfoStructs');
if strcmp(label, 'neurons')
    folderResults = fullfile(folderPC, 'Lab\RESULTS\nonvisualEffects\modelGratingResp');
    
    corrections = [];
else
    folderResults = fullfile(folderPC, 'Lab\RESULTS\boutons\nonvisualEffects');
    
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

exampleIndices = zeros(size(examples,1),1);
for ex = 1:size(examples,1)
    s = find(strcmp(tuning(examples(ex,1)).subject, {corrsRun.subject}) & ...
        strcmp(tuning(examples(ex,1)).date, {corrsRun.date}));
    exampleIndices(ex) = find(dataset==s & plane==examples(ex,2) & ...
        neuron==examples(ex,3));
end

%% Plot correlations with pupil/running against tuning modulation for cell types
rhos = {rhosRun, rhosPupil};
rho_pValues = {pValues_rhoRun, pValues_rhoPupil};
nonvisNames = {'running','pupil'};

stims = [3 1]; % use responses during darkness for correlations with running
               % use responses during gratings for correlations with pupil
colors = lines(4);
xLimits = [-.65 .85];
yLimits = [-.45 .45];
for nonvis = 1:2
    st = stims(nonvis);
    if all(isnan(rhos{nonvis}(:,st)))
        continue
    end
    
    figure
    hold on
    r = rhos{nonvis}(:, st);
    r(r > xLimits(2)) = xLimits(2);
    r(r < xLimits(1)) = xLimits(1);
    d = DImaxisAll;
    d(d > yLimits(2)) = yLimits(2);
    d(d < yLimits(1)) = yLimits(1);
    scatter(r, d, [], 'k', 'filled', 'MarkerFaceAlpha', 0.3);
    for ex = 1:length(exampleIndices)
        plot(rhos{nonvis}(exampleIndices(ex), st), ...
            DImaxisAll(exampleIndices(ex)), 'o', 'MarkerSize', 10, ...
            'MarkerFaceColor', colors(ex,:), 'MarkerEdgeColor', 'none')
    end
    plot([0 0], yLimits, 'k:')
    plot(xLimits, [0 0], 'k:')
    axis([xLimits yLimits])
    axis square
    xlabel(sprintf('Corr. with %s during %s', nonvisNames{nonvis}, ...
        stimuli{st}))
    ylabel('DI of resp @ pref direction')
    title(sprintf('n = %d', sum(~isnan(rhos{nonvis}(:, st)) & ...
        ~isnan(DImaxisAll))))
end
