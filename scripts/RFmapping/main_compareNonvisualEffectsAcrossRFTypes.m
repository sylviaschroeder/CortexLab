%% Folders
% file location depending on PCs
% folderPC = 'C:\STORAGE\OneDrive - University College London';
folderPC = 'C:\Users\skgtchr\OneDrive - University College London'; % ZEN notebook
% folderPC = 'C:\Storage\OneDrive - University College London'; % Lenovo laptop
% data folders
folderRFResults = fullfile(folderPC, 'Lab\RESULTS\receptiveFields\SC neurons');
% folderRFResults = fullfile(folderPC, 'Lab\RESULTS\receptiveFields\boutons');
folderCorrResults = fullfile(folderPC, 'Lab\RESULTS\nonvisCorrualEffects\modelGratingResp');
% folderCorrResults = fullfile(folderPC, 'Lab\RESULTS\boutons\nonvisCorrualEffects');

%% Parameters
stimuli = {'gratings', 'grayScreen', 'dark'};
modelNames = {'no RF','linear','absolute','On','Off'};
% label = 'neurons';
label = 'boutons';
nonvisCorr = 'pupil';
% nonvisCorr = 'running';
nonVisTun = 'pupil';

%% Collect relevant variables from visual noise data
% for running + stimulus, select best stimulus model (linear, absolute,
% white, or black)
file = 'receptiveFields.mat';
% file = 'receptiveFields1-5.mat';
data = load(fullfile(folderRFResults, file));
RFs = data.RFs;

dataset = [];
plane = [];
neuron = [];
evRun = [];
evTotal = [];
evStimOnly = [];
confIntStimOnly = [];
lambdasStim = [];
absSizeRun = [];
sizeRun = [];
signRun = [];
signStim = [];

for k = 1:length(RFs)
    if isempty(RFs(k).RFTimes)
        continue
    end
    runBin = diff(RFs(k).runWindow(1:2));
    for iPlane = 1:length(RFs(k).planes)
        evRun = [evRun; RFs(k).plane(iPlane).explainedVariances_runOnly];
        evTotal = [evTotal; RFs(k).plane(iPlane).explainedVariances];
        evStimOnly = [evStimOnly; RFs(k).plane(iPlane).explainedVariances_stimOnly];
        confIntStimOnly = [confIntStimOnly; ...
            RFs(k).plane(iPlane).explainedVariances_stimOnly_confInt];
        lambdasStim = [lambdasStim; RFs(k).plane(iPlane).lambdas];
        
        runKs = RFs(k).plane(iPlane).runningKernels; %[t x neuron]
        absSizeRun = [absSizeRun; mean(abs(runKs),1)'./runBin];
        sizeRun = [sizeRun; mean(runKs,1)'./runBin];
        signRun = [signRun; sign(sum(runKs,1))'];
        
        rfs = RFs(k).plane(iPlane).receptiveFields; %[row x column x t x neuron x model]
        rfs = reshape(rfs, [], size(rfs,4), size(rfs,5)); %[pixels x neuron x model]
        [~,maxPx] = max(abs(rfs), [], 1);
        linindsLin = sub2ind(size(rfs), maxPx(:), ...
            repmat((1:size(rfs,2))',size(rfs,3),1), ...
            reshape(repmat(1:size(rfs,3),size(rfs,2),1),[],1));
        signStim = [signStim; reshape(sign(rfs(linindsLin)), size(rfs,2), size(rfs,3))];
        
        dataset = [dataset; ones(size(runKs,2),1).*k];
        plane = [plane; ones(size(runKs,2),1).*iPlane];
        neuron = [neuron; (1:size(runKs,2))'];
    end
end
% switch sign of "black" stimulus model because negative value means the
% neuron is driven by black
signStim(:,end) = -signStim(:,end);

% set explained variance to 0 if not significant
ind = evStimOnly < confIntStimOnly(:,:,2);
trueEVstim = evStimOnly;
trueEVstim(ind) = 0;
trueEVTotal = evTotal;
r = repmat(evRun,1,size(evTotal,2));
trueEVTotal(ind) = r(ind);

% get exaplained variance and sign of RF for best stimulus model (for each
% neuron)
[maxEVStim, model] = max(trueEVstim, [], 2);
ind = sub2ind(size(trueEVstim), (1:size(trueEVTotal,1))', model);
maxTotal = trueEVTotal(ind);
maxSignStim = signStim(ind);
bestLams = lambdasStim(ind);

% for "good RF", lambda for RF should be < 1 and explained variance should
% be > 1.5%
model(bestLams >= 1) = 0;
model(maxEVStim < 0.015) = 0;

% % find units where EV for best model is with confidence interval of
% % shuffled data but best lambda is < 1 and EV is larger 1.5%
% % in other words: shuffled data says RF is bad, best lambda and EV say it's
% % good
% n = [];
% [mx, mdl] = max(evStimOnly, [], 2);
% for m = 1:4
%     ind = mdl==m & mx<confIntStimOnly(:,m,2) & lambdasStim(:,m)<1 & mx>0.015;
%     n = [n; find(ind)];
% end
% n = sort(n);
% % find units where EV for best model is above confidence interval of
% % shuffled data and EV is larger 1.5%, but best lambda is >= 1 
% % in other words: shuffled data says RF is good, best lambda says it's bad
% n = [];
% [mx, mdl] = max(evStimOnly, [], 2);
% for m = 1:4
%     ind = mdl==m & mx>confIntStimOnly(:,m,2) & lambdasStim(:,m)>=1 & mx>0.015;
%     n = [n; find(ind)];
% end
% n = sort(n);

%% Collect correlation values from analyses of continuous traces
data = load(fullfile(folderCorrResults, nonvisCorr, ...
    'corrsDuringGratingsAndGrayScreen_sigma1.00.mat'));
corrs = data.corrs;

rhos = [];
nullRhos = [];
resp = [];
isGad = [];
n = 0;
for k = 1:length(RFs)
    kcorr = find(strcmp(RFs(k).subject, {corrs.subject}) & ...
        strcmp(RFs(k).date, {corrs.date}));
    if isempty(kcorr)
        numCells = sum(cellfun(@length, {RFs(k).plane.explainedVariances_runOnly}));
        rhos = [rhos; NaN(numCells, length(stimuli))];
        nullRhos = [nullRhos; NaN(numCells, length(stimuli), ...
            size(corrs(1).plane(1).gratings.nullRhos,2))];
        resp = [resp; NaN(numCells, 1)];
        if isfield(corrs(1).plane, 'isGad')
            isGad = [isGad; NaN(numCells, 1)];
        end
        n = n + numCells;
        continue
    end
    for iPlane = 1:length(corrs(kcorr).plane)
        fields = fieldnames(corrs(kcorr).plane(iPlane));
        st = find(~structfun(@isempty, corrs(kcorr).plane(iPlane)),1);
        numCells = length(corrs(kcorr).plane(iPlane).(fields{st}).rhos);
        rhos = [rhos; NaN(numCells, length(stimuli))];
        nullRhos = [nullRhos; NaN(numCells, length(stimuli), ...
            size(corrs(1).plane(1).gratings.nullRhos,2))];
        for st = 1:length(stimuli)
            if ~isfield(corrs(kcorr).plane(iPlane), stimuli{st})
                continue
            end
            rhos((1:numCells) + n, st) = ...
                corrs(kcorr).plane(iPlane).(stimuli{st}).rhos;
            nullRhos((1:numCells) + n, st, :) = ...
                corrs(kcorr).plane(iPlane).(stimuli{st}).nullRhos;
        end
        if isfield(corrs(kcorr).plane(iPlane), 'responsive')
            resp = [resp; corrs(kcorr).plane(iPlane).responsive];
        else resp = [resp; NaN(numCells,1)];
        end
        if isfield(corrs(kcorr).plane, 'isGad')
            g = corrs(kcorr).plane(iPlane).isGad;
            isGad = [isGad; g];
        end
        
        n = n + numCells;
    end 
end

%% Collect tuning data
data = load(fullfile(folderCorrResults, nonvisTun, ...
    'tuning_prefDirSigmaDIFixed_lineFit.mat'));
tuning = data.tuning;
degrees = data.x;
data = load(fullfile(folderCorrResults, nonvisTun, ...
    'nullTuning_prefDirSigmaDIFixed.mat'));
null = data.null;

minis = []; % resp to non-pref stimulus
maxis = []; % resp to pref stimulus
stimMeans = [];
nullMinis = [];
nullMaxis = [];
nullStimMeans = [];
draws = cellfun(@size,squeeze(struct2cell(null(1).plane(1).cond(1).cell)),'UniformOutput',false);
draws = max(cell2mat(draws),[],1);
draws = draws(1);
fprintf('Dataset (of %d): ', length(RFs))
for k = 1:length(RFs)
    fprintf('%d ', k)
    if isempty(RFs(k).RFTimes)
        continue
    end
    ktun = find(strcmp(RFs(k).subject, {tuning.subject}) & ...
        strcmp(RFs(k).date, {tuning.date}));
    for iPlane = 1:length(tuning(ktun).plane)
        numCells = sum(dataset==k & plane==iPlane);
        data = tuning(ktun).plane(iPlane);
        mn = NaN(numCells,2);
        mx = NaN(numCells,2);
        sm = NaN(numCells,2);
        nmn = NaN(numCells,2,draws);
        nmx = NaN(numCells,2,draws);
        nsm = NaN(numCells,2,draws);
        if isempty(data.cellIDs)
            minis = [minis; mn];
            maxis = [maxis; mx];
            stimMeans = [stimMeans; sm];
            nullMinis = [nullMinis; nmn];
            nullMaxis = [nullMaxis; nmx];
            nullStimMeans = [nullStimMeans; nsm];
            continue
        end
        
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
    end
end
fprintf('\n')

mods = abs(maxis - minis);
nullMods = abs(nullMaxis - nullMinis);

%% Plot size of running kernel versus explained variance
figure
scatter(absSizeRun, evRun, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.1)
xlabel('Size of running kernel (mean absolute value per s)')
ylabel('Explained variance by running')

figure
scatter(sizeRun, evRun, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.1)
xlabel('Size of running kernel (mean value per s)')
ylabel('Explained variance by running')

%% Plot explained variance: running only versus running+stimulus or +stimulus

f1 = figure;
f2 = figure;
cols = [0 0 0;lines(size(evTotal,2))];
limx1 = [min(evRun .* signRun), max(evRun .* signRun)];
limx2 = [min(0, min(evRun)), max(evRun)];
limy1 = [min(maxEVStim .* maxSignStim), max(maxEVStim .* maxSignStim)];
limy2 = [min(maxEVStim), max(maxEVStim)];
plotPos = [3 1 2 4 5];
for m = 0:max(model)
    ind = model == m;
    
    figure(f1)
    subplot(2,3,plotPos(m+1))
    hold on
    plot(limx1, [0 0], 'k')
    plot([0 0], limy1, 'k')
    h = scatter(evRun(ind) .* signRun(ind), maxEVStim(ind) .* maxSignStim(ind), ...
        15, cols(m+1,:), 'filled', 'MarkerFaceAlpha', .3);
    xlim(limx1)
    ylim(limy1)
    title(modelNames{m+1})
    if m == 3
        xlabel('EV running only')
        ylabel('EV added by stim')
    end
    row = dataTipTextRow('Dataset',dataset(ind));
    h.DataTipTemplate.DataTipRows(end+1) = row;
    row = dataTipTextRow('Plane',plane(ind));
    h.DataTipTemplate.DataTipRows(end+1) = row;
    row = dataTipTextRow('Neuron',neuron(ind));
    h.DataTipTemplate.DataTipRows(end+1) = row;
    
    figure(f2)
    subplot(2,3,plotPos(m+1))
    hold on
    plot(limx1, [0 0], 'k')
    plot([0 0], limy1, 'k')
    h = scatter(evRun(ind), maxEVStim(ind), 15, cols(m+1,:), 'filled', ...
        'MarkerFaceAlpha', 0.3);
    xlim(limx2)
    ylim(limy2)
    title(modelNames{m+1})
    if m == 3
        xlabel('EV running only')
        ylabel('EV added by stim')
    end
    row = dataTipTextRow('Dataset',dataset(ind));
    h.DataTipTemplate.DataTipRows(end+1) = row;
    row = dataTipTextRow('Plane',plane(ind));
    h.DataTipTemplate.DataTipRows(end+1) = row;
    row = dataTipTextRow('Neuron',neuron(ind));
    h.DataTipTemplate.DataTipRows(end+1) = row;
end
savefig(f1, fullfile(folderRFResults, 'EVrunningPerStimModel_datasets.fig'))

%% Plot running (expl. var. or kernel size) versus On-Off-continuum
% On-Off-continuum: EV(on) / (EV(on) + EV(off)) -> if 0, then completely
% off; if 1, then completely on
% mark whether RF is more linear or complex (absolute) depending on which
% EV for two models is larger
onOffRatio = trueEVstim(:,3) ./ sum(trueEVstim(:,3:4),2);
hasRF = any(trueEVstim > 0, 2);
isLin = trueEVstim(:,1) > trueEVstim(:,2);

cols = lines(2);
indsLin = {hasRF & ~isLin, isLin};
indsStim = {maxSignStim > 0, maxSignStim < 0};
titles = {'Enhanced neurons', 'Suppressed neurons'};
for s = 1:2
    figure
    hold on
    for k = 1:2
        h = scatter(evRun(indsLin{k} & indsStim{s}) .* ...
            signRun(indsLin{k} & indsStim{s}), ...
            onOffRatio(indsLin{k} & indsStim{s}), ...
            15, cols(k,:), 'filled', 'MarkerFaceAlpha', .3);
        row = dataTipTextRow('Dataset',dataset(indsLin{k} & indsStim{s}));
        h.DataTipTemplate.DataTipRows(end+1) = row;
        row = dataTipTextRow('Plane',plane(indsLin{k} & indsStim{s}));
        h.DataTipTemplate.DataTipRows(end+1) = row;
        row = dataTipTextRow('Neuron',neuron(indsLin{k} & indsStim{s}));
        h.DataTipTemplate.DataTipRows(end+1) = row;
    end
    xlabel('Explained variance running (* sign of kernel)')
    ylabel('On/(On+Off)-ratio of RF type')
    legend('absolute RF','linear RF')
    title(titles{s})
end

% histograms for RF types On, Off, absolute, linear
cols = [0 0 0; lines(4)];
figure
hold on
plot([0 0],[0 1],'k')
h = zeros(1, size(evStimOnly,2)+1);
for m = 1:size(evStimOnly,2)+1
    ind = model == (m-1) & ~isnan(evRun);
    x = sort(evRun(ind) .* signRun(ind));
    y = (1:length(x))' ./ length(x);
    x = [min(evRun .* signRun); x; max(evRun .* signRun)];
    y = [0; y; 1];
    h(m) = plot(x, y, 'Color', cols(m,:), 'LineWidth', 2);
end
xlim([min(evRun .* signRun) max(evRun .* signRun)])
ylim([0 1])
xlabel('Explained variance running only')
ylabel('Proportion of units')
legend(h, modelNames, 'Location', 'SouthEast')

% histograms for RF types On enhanced, On suppressed, Off enhanced, Off
% suppressed, absolute enhanced, absolute suppressed
cols = lines(4);
linType = {'-',':'};
names = {'Abs enh', 'Abs supp', 'On enh', 'On supp', 'Off enh', 'Off supp', 'no RF'};
h = zeros(1, 7);
figure
hold on
plot([0 0],[0 1],'k')
x = sort(evRun(model==0) .* signRun(model==0));
y = (1:length(x))' ./ length(x);
x = [min(evRun .* signRun); x; max(evRun .* signRun)];
y = [0; y; 1];
h(end) = plot(x, y, 'k', 'LineWidth', 2);
k = 1;
for m = [2 3 4]
    for r = 1:2
        ind = model == m & indsStim{r} & ~isnan(evRun);
        x = sort(evRun(ind) .* signRun(ind));
        y = (1:length(x))' ./ length(x);
        x = [min(evRun .* signRun); x; max(evRun .* signRun)];
        y = [0; y; 1];
        h(k) = plot(x, y, linType{r}, 'Color', cols(m,:), 'LineWidth', 2);
        k = k + 1;
    end
end
xlim([min(evRun .* signRun) max(evRun .* signRun)])
ylim([0 1])
xlabel('Explained variance running only')
ylabel('Proportion of units')
legend(h, names, 'Location', 'SouthEast')

%% Plot RF types for excitatory and inhibitory neurons
onOffRatio = trueEVstim(:,3) ./ sum(trueEVstim(:,3:4),2);
indsStim = {maxSignStim > 0, maxSignStim < 0};
indCellType = {~(isGad==1), isGad==1};
bins = 0 : 0.05 : 1;
edges = [bins, 1.1] - 0.025;
cols = lines(2);
linType = {'-',':'};
figure
hold on
h = [0 0 0 0];
for s = 1:2
    for t = 1:2
        n = histcounts(onOffRatio(indCellType{t} & indsStim{s}), edges);
        n = n ./ sum(n);
        h((s-1)*2+t) = plot(bins, n, linType{s}, 'Color', cols(t,:), ...
            'LineWidth', 1);
    end
end
legend(h, {'exc+enh','inh+enh','exc+supp','inh+supp'})
xlabel('RF type: (On/(On+Off))')
ylabel('Proportion of neurons')

%% Plot running/pupil correlation for different RF types
indsEnh = {maxSignStim > 0, maxSignStim < 0};

% Cum. histogram of correlation coefficients (data and null distribution)
% for each RF type (On enh., On supp., Off enh., Off supp., ...)
RFTypes = {'Linear', 'Absolute', 'On', 'Off'};
respTypes = {'enhanced', 'suppressed'};
xLimits = [floor(min(rhos(:))*10)/10, ceil(max(rhos(:))*10)/10];
for stim = 1:length(stimuli)
    if all(isnan(rhos(:,stim)))
        continue
    end
    xAll = sort(rhos(~isnan(rhos(:,stim)),stim));
    yAll = (1:length(xAll))' ./ length(xAll);
    xAll = [-1; xAll; 1];
    yAll = [0; yAll; 1];
    
    ind = model==0 & ~isnan(rhos(:,stim));
    general.plotCumHist(rhos(ind,stim), squeeze(nullRhos(ind,stim,:)), ...
        xAll, yAll)
    xlim(xLimits)
    xlabel(['Correlation with ' nonvisCorr])
    ylabel(['Proportion of ' label])
    title(sprintf('No RF, correlation during %s (n = %d)', stimuli{stim}, sum(ind)))
    for m = 1:4
        for r = 1:2
            ind = model==m & indsEnh{r} & ~isnan(rhos(:,stim));
            general.plotCumHist(rhos(ind,stim), squeeze(nullRhos(ind,stim,:)), ...
                xAll, yAll)
            xlim(xLimits)
            xlabel(['Correlation with ' nonvisCorr])
            ylabel(['Proportion of ' label])
            title(sprintf('%s %s, correlation during %s (n = %d)', RFTypes{m}, ...
                respTypes{r}, stimuli{stim}, sum(ind)))
        end
    end
end


% On-Off-continuum: EV(on) / (EV(on) + EV(off)) -> if 0, then completely
% off; if 1, then completely on
% mark whether RF is more linear or complex (absolute) depending on which
% EV for two models is larger
onOffRatio = trueEVstim(:,3) ./ sum(trueEVstim(:,3:4),2);

cols = lines(4);
cols(1:2,:) = [];
titles = {'Enhanced neurons', 'Suppressed neurons'};
for stim = 1:length(stimuli)
    figure
    hold on
    for s = 1:2
        h = scatter(rhos(indsEnh{s},stim), ...
            onOffRatio(indsEnh{s}), ...
            15, cols(s,:), 'filled', 'MarkerFaceAlpha', .3);
        row = dataTipTextRow('Dataset',dataset(indsEnh{s}));
        h.DataTipTemplate.DataTipRows(end+1) = row;
        row = dataTipTextRow('Plane',plane(indsEnh{s}));
        h.DataTipTemplate.DataTipRows(end+1) = row;
        row = dataTipTextRow('Neuron',neuron(indsEnh{s}));
        h.DataTipTemplate.DataTipRows(end+1) = row;
    end
    xlabel('Correlation with pupil')
    ylabel('On/(On+Off)-ratio of RF type')
    legend('enhanced','suppressed')
    title(stimuli{stim})
end

% histograms for RF types On, Off, absolute, linear
cols = [0 0 0; lines(4)];
binEdg = floor(min(rhos(:))*10)/10 : 0.02 : ceil(max(rhos(:))*10)/10;
binCen = binEdg(1:end-1) + 0.01;
for stim = 1:length(stimuli)
    figure
    hold on
    plot([0 0],[0 1],'k')
    h = zeros(1, size(evStimOnly,2)+1);
    for m = 1:size(evStimOnly,2)+1
        ind = model == (m-1);
        n = histcounts(rhos(ind,stim), binEdg);
        n = cumsum(n ./ sum(n));
        h(m) = plot(binCen, n, 'Color', cols(m,:), 'LineWidth', 2);
    end
    xlim(binEdg([1 end]))
    ylim([0 1])
    xlabel('Correlation with pupil')
    ylabel('Proportion of units')
    legend(h, modelNames, 'Location', 'SouthEast')
    title(stimuli{stim})
end

% histograms for RF types On enhanced, On suppressed, Off enhanced, Off
% suppressed, absolute enhanced, absolute suppressed
cols = lines(4);
binEdg = floor(min(rhos(:))*10)/10 : 0.02 : ceil(max(rhos(:))*10)/10;
binCen = binEdg(1:end-1) + 0.01;
linType = {'-',':'};
names = {'Abs enh', 'Abs supp', 'On enh', 'On supp', 'Off enh', 'Off supp', 'no RF'};
for stim = 1:length(stimuli)
    h = zeros(1, 7);
    figure
    hold on
    plot([0 0],[0 1],'k')
    n = histcounts(rhos(model==0,stim), binEdg);
    n = cumsum(n ./ sum(n));
    h(end) = plot(binCen, n, 'k', 'LineWidth', 2);
    k = 1;
    for m = [2 3 4]
        for r = 1:2
            n = histcounts(rhos(model==m & indsEnh{r},stim), binEdg);
            n = cumsum(n ./ sum(n));
            h(k) = plot(binCen, n, linType{r}, 'Color', cols(m,:), 'LineWidth', 2);
            k = k + 1;
        end
    end
    xlim(binEdg([1 end]))
    ylim([0 1])
    xlabel('Correlation with pupil')
    ylabel('Proportion of units')
    legend(h, names, 'Location', 'SouthEast')
    title(stimuli{stim})
end

%% Plot tuning modulation for different RF types
indsEnh = {maxSignStim > 0, maxSignStim < 0};
RFTypes = {'Linear', 'Absolute', 'On', 'Off'};
respTypes = {'enhanced', 'suppressed'};
limits = [0.6 0.7];

% Cum. histogram of difference index (DI) for prefered orientation and tuning depth
% for each RF type (On enh., On supp., Off enh., Off supp., ...)
ind = all(isnan(maxis),2); % untuned neurons
mx = maxis;
mx(ind,:) = stimMeans(ind,:);
nmx = nullMaxis;
nmx(ind,:,:) = nullStimMeans(ind,:,:);
measures = {mx, mods};
nullMeasures = {nmx, nullMods};
modFun = @(a,b) (b-a)./(abs(a)+abs(b));
funLabel = 'DI (large-small)/(large+small)';
measureLabels = {'Resp @ pref dir','Tuning depth'};

for m = 1:length(measures)
    DIs = modFun(measures{m}(:,1), measures{m}(:,2));
    nullDIs = modFun(squeeze(nullMeasures{m}(:,1,:)), ...
        squeeze(nullMeasures{m}(:,2,:)));
    
    xAll = sort(DIs(~isnan(DIs)));
    yAll = (1:length(xAll))' ./ length(xAll);
    xAll = [-1; xAll; 1];
    yAll = [0; yAll; 1];
    
    ind = model==0 & ~isnan(DIs);
    general.plotCumHist(DIs(ind), nullDIs(ind,:), xAll, yAll)
    xlim([-limits(m) limits(m)])
    xlabel(funLabel)
    ylabel(['Proportion of ' label])
    title(sprintf('No RF, DI for %s (n = %d)', measureLabels{m}, sum(ind)))
    for rf = 1:4
        for r = 1:2
            ind = model==rf & indsEnh{r} & ~isnan(DIs);
            general.plotCumHist(DIs(ind), nullDIs(ind,:), xAll, yAll)
            xlim([-limits(m) limits(m)])
            xlabel(funLabel)
            ylabel(['Proportion of ' label])
            title(sprintf('%s %s, DI for %s (n = %d)', RFTypes{rf}, ...
                respTypes{r}, measureLabels{m}, sum(ind)))
        end
    end
end

%% Pie chart of RF types
labels = {'Linear','Abs+','Abs-','On+','On-','Off+','Off-','no RF'};
indsEnh = {maxSignStim > 0, maxSignStim < 0};
nums = zeros(length(labels)+1,1);
for m = 1:4
    for r = 1:2
        nums((m-1)*2+r) = sum(model == m & indsEnh{r});
    end
end
nums(end-1) = sum(model == 0);
nums(2) = sum(nums(1:2));
nums(1) = [];
figure
pie(nums, labels)

labels = {'RF+ gratings+','RF+ gratings-','RF- gratings+','RF- gratings-'};
nums = [sum(model>0&resp==1), sum(model>0&resp==0), sum(model==0&resp==1), ...
    sum(model==0&resp==0)];
figure
pie(nums, labels)
