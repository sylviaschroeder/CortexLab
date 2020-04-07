%% Parameters
stimuli = {'gratings', 'grayScreen', 'dark'};
% label = 'neurons';
label = 'boutons';
nonvisCorr = 'pupil';
% nonvisCorr = 'running';
nonvisTun = 'pupil';

%% Folders
% file location depending on PCs
folderPC = 'C:\STORAGE\OneDrive - University College London';
% folderPC = 'C:\Users\skgtchr\OneDrive - University College London'; % ZEN notebook
% folderPC = 'C:\Storage\OneDrive - University College London'; % Lenovo laptop
% data folders
if strcmp(label, 'neurons')
    folderResults = fullfile(folderPC, 'Lab\RESULTS\nonvisualEffects\modelGratingResp');
    folderRFResults = fullfile(folderPC, 'Lab\RESULTS\receptiveFields\SC neurons');
else
    folderResults = fullfile(folderPC, 'Lab\RESULTS\boutons\nonvisualEffects');
    folderRFResults = fullfile(folderPC, 'Lab\RESULTS\receptiveFields\boutons');
end

%% Collect correlation values from analyses of continuous traces
data = load(fullfile(folderResults, nonvisCorr, ...
    'corrsDuringGratingsAndGrayScreen_sigma1.00.mat'));
corrs = data.corrs;

dataset = [];
plane = [];
neuron = [];
rhos = [];
nullRhos = [];
resp = [];
isGad = [];
n = 0;
for k = 1:length(corrs)
    for iPlane = 1:length(corrs(k).plane)
        fields = fieldnames(corrs(k).plane(iPlane));
        st = find(~structfun(@isempty, corrs(k).plane(iPlane)),1);
        numCells = length(corrs(k).plane(iPlane).(fields{st}).rhos);
        dataset = [dataset; ones(numCells,1).*k];
        plane = [plane; ones(numCells,1).*iPlane];
        neuron = [neuron; (1:numCells)'];
        rhos = [rhos; NaN(numCells, length(stimuli))];
        nullRhos = [nullRhos; NaN(numCells, length(stimuli), ...
            size(corrs(1).plane(1).gratings.nullRhos,2))];
        for st = 1:length(stimuli)
            if ~isfield(corrs(k).plane(iPlane), stimuli{st})
                continue
            end
            rhos((1:numCells) + n, st) = ...
                corrs(k).plane(iPlane).(stimuli{st}).rhos;
            nullRhos((1:numCells) + n, st, :) = ...
                corrs(k).plane(iPlane).(stimuli{st}).nullRhos;
        end
        if isfield(corrs(k).plane(iPlane), 'responsive')
            resp = [resp; corrs(k).plane(iPlane).responsive];
        else resp = [resp; NaN(numCells,1)];
        end
        if isfield(corrs(k).plane, 'isGad')
            isGad = [isGad; corrs(k).plane(iPlane).isGad];
        else isGad = [isGad; NaN(numCells,1)];
        end
        
        n = n + numCells;
    end 
end

%% Collect tuning data
data = load(fullfile(folderResults, nonvisTun, ...
    'tuning_prefDirSigmaDIFixed_lineFit.mat'));
tuning = data.tuning;
degrees = data.x;
data = load(fullfile(folderResults, nonvisTun, ...
    'nullTuning_prefDirSigmaDIFixed.mat'));
null = data.null;

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
fprintf('Dataset (of %d): ', length(corrs))
for k = 1:length(corrs)
    fprintf('%d ', k)
    ktun = find(strcmp(corrs(k).subject, {tuning.subject}) & ...
        strcmp(corrs(k).date, {tuning.date}));
    if isempty(ktun) || isempty(tuning(ktun).plane)
        numCells = sum(dataset==k);
        minis = [minis; NaN(numCells,2)];
        maxis = [maxis; NaN(numCells,2)];
        stimMeans = [stimMeans; NaN(numCells,2)];
        nullMinis = [nullMinis; NaN(numCells,2,draws)];
        nullMaxis = [nullMaxis; NaN(numCells,2,draws)];
        nullStimMeans = [nullStimMeans; NaN(numCells,2,draws)];
        isSuppr = [isSuppr; NaN(numCells,1)];
        continue
    end
    for iPlane = 1:length(tuning(ktun).plane)
        numCells = sum(dataset==k & plane==iPlane);
        data = tuning(ktun).plane(iPlane);
        mn = NaN(numCells,2);
        mx = NaN(numCells,2);
        sm = NaN(numCells,2);
        nmn = NaN(numCells,2,draws);
        nmx = NaN(numCells,2,draws);
        nsm = NaN(numCells,2,draws);
        spp = NaN(numCells,1);
        if isempty(data.cellIDs)
            minis = [minis; mn];
            maxis = [maxis; mx];
            stimMeans = [stimMeans; sm];
            nullMinis = [nullMinis; nmn];
            nullMaxis = [nullMaxis; nmx];
            nullStimMeans = [nullStimMeans; nsm];
            isSuppr = [isSuppr; spp];
            continue
        end
        
        spp(data.cellIDs) = data.isSuppressed;
        neurons = find(~cellfun(@isempty, {data.cond(1).cell.parameters}));
        for j = 1:length(neurons)
            iCell = neurons(j);
            ID = data.cellIDs(iCell);
            for c = 1:2
                pars = data.cond(c).cell(iCell).parameters;
                curve = data.cond(c).cell(iCell).curve;
                if length(pars) == 1 % not tuned
                    sm(ID,c) = pars;
                else
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
        end
        minis = [minis; mn];
        maxis = [maxis; mx];
        stimMeans = [stimMeans; sm];
        nullMinis = [nullMinis; nmn];
        nullMaxis = [nullMaxis; nmx];
        nullStimMeans = [nullStimMeans; nsm];
        isSuppr = [isSuppr; spp];
    end
%     disp(size(minis,1))
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
DImaxisAll = modFun(mx(:,1), mx(:,2));
nullDImaxisAll = modFun(squeeze(nmx(:,1,:)), squeeze(nmx(:,2,:)));
DImods = modFun(mods(:,1), mods(:,2));
nullDImods = modFun(squeeze(nullMods(:,1,:)), squeeze(nullMods(:,2,:)));

% for suppressed units, change sign of DIs of preferred responses so that
% a more negative response with arousal results in a positive DI
DImaxis(isSuppr==1) = -DImaxis(isSuppr==1);
nullDImaxis(isSuppr==1,:) = -nullDImaxis(isSuppr==1,:);
DImaxisAll(isSuppr==1) = -DImaxisAll(isSuppr==1);
nullDImaxisAll(isSuppr==1,:) = -nullDImaxisAll(isSuppr==1,:);

%% Collect relevant variables from visual noise data
file = 'receptiveFields.mat';
data = load(fullfile(folderRFResults, file));
RFs = data.RFs;
file = 'receptiveFields_stim+run_4RFmodels.mat';
data = load(fullfile(folderRFResults, file));
RFs_old = data.RFs;

evStimOnly = NaN(0,4);
confIntStimOnly = NaN(0,4,2);
lambdasStim = NaN(0,1);
onOffPeaks = NaN(0,2);
for k = 1:length(corrs)
%     disp(size(evStimOnly,1))
    krf = find(strcmp(corrs(k).subject, {RFs.subject}) & ...
        strcmp(corrs(k).date, {RFs.date}));
    if isempty(krf) || isempty(RFs(krf).RFTimes) || sum(dataset == k)==0
        numCells = sum(dataset == k);
        evStimOnly(end+(1:numCells),:) = NaN;
        confIntStimOnly(end+(1:numCells),:,:) = NaN;
        lambdasStim(end+(1:numCells),:) = NaN;
        onOffPeaks(end+(1:numCells),:) = NaN;
        continue
    end
    for iPlane = 1:length(RFs(krf).planes)
        evStimOnly = [evStimOnly; RFs_old(krf).plane(iPlane).explainedVariances_stimOnly];
        confIntStimOnly = [confIntStimOnly; ...
            RFs_old(krf).plane(iPlane).explainedVariances_stimOnly_confInt];
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

goodRF = any(evStimOnly > confIntStimOnly(:,:,2) & ...
    evStimOnly > 0.015, 2) & lambdasStim < 1;

% get ON-OFF-ratio
[~,type] = max(abs(onOffPeaks),[],2); % determine whether On or Off is stronger
inds = sub2ind(size(onOffPeaks), (1:size(onOffPeaks,1))', type);
signs = sign(onOffPeaks(inds)); % determine whether neuron is driven or suppressed
OnOffratio = onOffPeaks;
OnOffratio(signs<0,:) = -OnOffratio(signs<0,:); % if suppressed, reverse values
OnOffratio(OnOffratio<0) = 0; % set values that have opposite drive to dominant RF type to zero
OnOffratio = OnOffratio(:,1) ./ sum(OnOffratio,2);

%% Plot correlations against tuning modulations
% Plot correlations against modulation of preferred responses (do not
% distinguish between cell types other than enhanced and suppressed)
% groupNames = {'All', 'Enhanced', 'Suppressed'};
% groupInds = [true(size(rhos,1),1), isSuppr==-1, isSuppr==1];
% groupColors = [[0 0 0];lines(2)];

groupNames = {'ON', 'ON/OFF', 'OFF'};
groupInds = [OnOffratio>0.8, OnOffratio>=0.2&OnOffratio<=0.8, OnOffratio<0.2];
groupColors = [1 0 0;.5 0 .5;0 0 1];


measures = {DImaxisAll, DImaxis, DImods; nullDImaxisAll, nullDImaxis, nullDImods};
measureNames = {'DI resp @ pref stim (all)', 'DI resp @ pref stim (tuned)', ...
    'DI tuning depth'};

figSizes = [[370 680 1540 420]; [370 290 1540 800]; [620 42 1287 1074]];
dotSizes = [10 20];
xLimits = [floor(min(rhos(:))*20)/20, ceil(max(rhos(:))*20)/20];

% do NOT distinguish between cell types
for stim = 1:length(stimuli)
    if all(isnan(rhos(:,stim)))
        continue
    end
    rhoCnf = squeeze(prctile(nullRhos(:,stim,:), [2.5 97.5], 3));
    for m = 1:length(measures)
        mCnf = prctile(measures{2,m}, [2.5 97.5], 2);
        isSgnf = (rhos(:,stim)<rhoCnf(:,1) | rhos(:,stim)>rhoCnf(:,2)) & ...
            (measures{1,m}<mCnf(:,1) | measures{1,m}>mCnf(:,2));
        isSgnf = [~isSgnf, isSgnf];
        yLimits = [floor(min(measures{1,m})*20)/20, ceil(max(measures{1,m})*20)/20];
        figure('Position', figSizes(1,:))
        annotation('textbox', [0.02,0.85,0.1,0.03], 'String', ...
            sprintf('%s\n%s\n%s', label, stimuli{stim}, measureNames{m}), ...
            'LineStyle', 'none', 'FontSize', 12, 'FontWeight', 'bold')
        for group = 1:length(groupNames)
            subplot(1,3,group)
            hold on
            plot([-1 1], [0 0], 'k')
            plot([0 0], [-1 1], 'k')
            for k = 1:2
                h = scatter(rhos(groupInds(:,group) & isSgnf(:,k),stim), ...
                    measures{1,m}(groupInds(:,group) & isSgnf(:,k)), ...
                    dotSizes(k), groupColors(group,:), 'filled', 'MarkerFaceAlpha', 0.3);
                row = dataTipTextRow('Subject',{corrs(dataset(groupInds(:,group) & isSgnf(:,k))).subject}');
                h.DataTipTemplate.DataTipRows(end+1) = row;
                row = dataTipTextRow('Date',{corrs(dataset(groupInds(:,group) & isSgnf(:,k))).date}');
                h.DataTipTemplate.DataTipRows(end+1) = row;
                row = dataTipTextRow('Plane',plane(groupInds(:,group) & isSgnf(:,k)));
                h.DataTipTemplate.DataTipRows(end+1) = row;
                row = dataTipTextRow('Neuron',neuron(groupInds(:,group) & isSgnf(:,k)));
                h.DataTipTemplate.DataTipRows(end+1) = row;
                h.DataTipTemplate.Interpreter = 'none';
            end
            xlim(xLimits)
            ylim(yLimits)
            title(sprintf('%s (n=%d)', groupNames{group}, ...
                sum(~any(isnan([rhos(groupInds(:,group),stim), measures{1,m}(groupInds(:,group))]),2))))
            if group == 1
                xlabel(sprintf('Correlation with %s during %s', nonvisCorr, stimuli{stim}))
                ylabel(measureNames{m})
            end
            axis square
        end
    end
end


cellTypes = {[isGad==-1, isGad==1], [OnOffratio<0.3, ...
    OnOffratio>=0.3&OnOffratio<=0.7, OnOffratio>0.7]};
cellTypeNames = {{'Excitatory', 'Inhibitory'},{'OFF','Absolute','ON'}};
cellTypeColors = lines(5);
cellTypeColors = cellTypeColors(3:end,:);

% DO distinguish between cell types
for stim = 1:length(stimuli)
    if all(isnan(rhos(:,stim)))
        continue
    end
    for m = 1:length(measures)
        yLimits = [floor(min(measures{m})*20)/20, ceil(max(measures{m})*20)/20];
        for class = 1:length(cellTypes)
            if ~any(cellTypes{class}(:))
                continue
            end
            rows = size(cellTypes{class},2);
            figure('Position', figSizes(rows,:))
            annotation('textbox', [0.02,0.97,0.2,0.03], 'String', ...
                sprintf('%s\n%s\n%s', label, stimuli{stim}, measureNames{m}), ...
                'LineStyle', 'none', 'FontSize', 11, 'FontWeight', 'bold')
            for type = 1:rows
                for group = 1:length(groupNames)
                    subplot(rows,3,(type-1)*length(groupNames)+group)
                    hold on
                    plot([-1 1], [0 0], 'k')
                    plot([0 0], [-1 1], 'k')
                    inds = groupInds(:,group) & cellTypes{class}(:,type);
                    h = scatter(rhos(inds,stim), measures{m}(inds), ...
                        15, cellTypeColors(type,:), 'filled', 'MarkerFaceAlpha', 0.3);
                    xlim(xLimits)
                    ylim(yLimits)
                    title(sprintf('%s (n=%d)', groupNames{group}, ...
                        sum(~any(isnan([rhos(inds,stim), measures{m}(inds)]),2))))
                    if group == 1
                        ylabel(cellTypeNames{class}{type}, 'FontWeight', 'bold')
                    end
                    if group == length(groupNames) && type == rows
                        xlabel(sprintf('Correlation with %s during %s', nonvisCorr, stimuli{stim}))
                        ylabel(measureNames{m})
                    end
                    row = dataTipTextRow('Subject',{corrs(dataset(inds)).subject}');
                    h.DataTipTemplate.DataTipRows(end+1) = row;
                    row = dataTipTextRow('Date',{corrs(dataset(inds)).date}');
                    h.DataTipTemplate.DataTipRows(end+1) = row;
                    row = dataTipTextRow('Plane',plane(inds));
                    h.DataTipTemplate.DataTipRows(end+1) = row;
                    row = dataTipTextRow('Neuron',neuron(inds));
                    h.DataTipTemplate.DataTipRows(end+1) = row;
                    h.DataTipTemplate.Interpreter = 'none';
                    axis square
                end
            end
        end
    end
end

%% Compare suppression-/enhancement-by-contrast across gratings and visual noise
figure
bar([sum(isSuppr==-1&signStim_bestModel==1), ...
    sum(isSuppr==1&signStim_bestModel==-1), ...
    sum(isSuppr==-1&signStim_bestModel==-1), ...
    sum(isSuppr==1&signStim_bestModel==1)])
ylabel(['#' label])
title(sprintf('%s: enhanced/suppressed by gratings/visual noise', label))
set(gca, 'XTickLabel', {'grat+ RF+','grat- RF-','grat+ RF-','grat- RF+'})
% make sure: bestModel set to NaN when no model is good