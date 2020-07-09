%% Define data

% label = 'neurons';
label = 'boutons';

if strcmp(label, 'neurons')
    % (1) SS038, 2015-02-17 (OS, neg. DSI and corr, inhib.)
    % (2) SS041, 2015-04-23 (DS, pos. DSI and corr, exc.)
    % (3) SS041, 2015-04-23 (OS+DS, pos. DSI and corr, exc.)
    % (4) SS038, 2015-02-17 (not tuned, suppressed, neg. DSI and corr, inhib.)
    examples = [1 1 227; 2 2 38; 2 2 151; 1 1 203];
    % % Neurons in sSC
%     exampleSet = 7; % SS047, 2015-11-23, 1
%     examples = [1 168; 1 207]; %; 2 71];
    % OLD:
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

data = load(fullfile(folderResults, 'pupil', ...
    'tuning_prefDirSigmaDIFixed_lineFit.mat'));
tuning = data.tuning;
degrees = data.x;

data = load(fullfile(folderResults, 'kernelFit', 'results.mat'));
results = data.results;

expField = 'exp';
if ~isfield(results,expField)
    expField = 'expGratings';
end

%% Plot traces of examples cell aligned to stimulus and tuning curves
buffer = 2; % in sec (before and after stim period)
cols = [0 0 0; 1 0 0];
lin = {'-','-'};

for ex = 1:size(examples,1)
    exSet = examples(ex,1);
    pl = examples(ex,2);
    cellID = examples(ex,3);
    cellInd = find(results(exSet).plane(pl).cellIDs == cellID);
    % load meta
    data=load(fullfile(folderROIData, results(exSet).subject, ...
        results(exSet).date, num2str(results(exSet).(expField)), ...
        sprintf('%s_%d_%s_2P_plane%03d_ROI.mat', results(exSet).date, ...
        results(exSet).(expField), results(exSet).subject, ...
        results(exSet).planes(pl))));
    meta = data.meta;
    meta.folderTL = strrep(meta.folderTL, 'ioo.ucl.ac.uk', ...
        'cortexlab.net');
    
    [stimTimes, stimSeq, stimMatrix, frameTimes, samplingRate] = ...
        ssLocal.getStimulusResponseInfo(meta);
    repetitions = sum(stimSeq.seq == 1);
    stimDurInFrames = round(sum(stimMatrix(1,:)) / repetitions);
    frameDur = median(diff(frameTimes));
    stimDur = stimDurInFrames * frameDur;
    offset = ceil(buffer / frameDur);
    respTime = (-offset:stimDurInFrames+offset-1) .* frameDur;
    
    trace = meta.F_final(:,cellID);
    prediction = results(exSet).plane(pl).kernelFit(cellInd).prediction;
    a = 1;
    b = 0;
    if ~isempty(corrections)
        a = corrections(exSet).plane(pl).a{tuning(exSet).exp}(cellID);
        b = corrections(exSet).plane(pl).b{tuning(exSet).exp}(cellID);
        trace = doCorrect(a,b,trace);
        prediction = doCorrect(a,b,prediction);
    end
    
    resp = squeeze(ssLocal.getTracesPerStimulus(trace, stimMatrix, [1 1] .* offset)); % [stimulus x trial x time]
    predResp = squeeze(ssLocal.getTracesPerStimulus(prediction, stimMatrix, [1 1] .* offset));
    
    respMean = NaN(size(resp,1)-1, 2, size(resp,3));
    predMean = NaN(size(resp,1)-1, 2, size(resp,3));
    for c = 1:2
        ind = ~isnan(tuning(exSet).plane(pl).cond(c) ...
            .cell(cellInd).responses);
        ind = repmat(ind, 1, 1, size(resp,3));
        temp = resp(1:end-1,:,:);
        temp(~ind) = NaN;
        respMean(:,c,:) = nanmean(temp, 2);
        temp = predResp(1:end-1,:,:);
        temp(~ind) = NaN;
        predMean(:,c,:) = nanmean(temp, 2);
    end
    
    mini = min([respMean(:); predMean(:)]);
    maxi = max([respMean(:); predMean(:)]);
    rng = maxi - mini;
    mini = mini - 0.05*rng;
    maxi = maxi + 0.05*rng;
    xDist = .5;
    traceDur = stimDur + 2*buffer;
    
    figure('Position',[3 640 1915 420])
    hold on
    h = [0 0 0 0];
    for c = 1:2
        x0 = 0;
        for st = 1:size(respMean,1)
            if c == 1
                fill([0 stimDur stimDur 0] + x0, ...
                    [mini mini maxi maxi], 'k', 'FaceColor', cols(c,:), ...
                    'FaceAlpha', 0.1, 'EdgeColor', 'none')
                plot(respTime([1 end]) + x0, [0 0], 'k:')
            end
            h1 = plot(respTime + x0, squeeze(respMean(st,c,:)), 'Color', cols(c,:));
            h2 = plot(respTime + x0, squeeze(predMean(st,c,:)), 'Color', cols(c,:), 'LineWidth', 2);
            x0 = x0 + traceDur + xDist;
            if st==1 && c==1
                h(1) = h2;
                h(2) = h1;
            elseif st==2
                if c==1
                    h(3) = h2;
                else
                    h(4) = h2;
                end
            end
        end
    end
    axis tight
    set(gca, 'XTick', [0 stimDur])
    legend(h, {'Prediction from kernel fit','Data','Small pupil trials', 'Large pupil trials'})
    xlabel('Stimuli')
    ylabel('\DeltaF/F')
    title(sprintf('%s, %s, plane %d, cell %d', tuning(exSet).subject, ...
        tuning(exSet).date, tuning(exSet).planes(pl), cellID))
    
    figure
    hold on
    h = [0 0];
    for c = 1:2
        curve = tuning(exSet).plane(pl).cond(c).cell(cellInd).curve;
        resp = tuning(exSet).plane(pl).cond(c).cell(cellInd).responses;
        if ~isempty(corrections)
            curve = doCorrect(a, b, curve);
            resp = doCorrect(a,b,resp);
        end
        h(c) = plot(degrees, curve, lin{c}, 'Color', cols(c,:), 'LineWidth',2);
        m = nanmean(resp,2);
        s = nanstd(resp,0,2) ./ sqrt(sum(~isnan(resp),2));
        errorbar([tuning(exSet).plane(pl).cond(c).cell(cellInd) ...
            .directions; 360], m([1:end 1]), s([1:end 1]), 'o', ...
            'Color', cols(c,:), 'CapSize', 2, 'MarkerFaceColor', cols(c,:))
    end
    plot([0 360], [0 0], 'k:', 'LineWidth', 2)
    set(gca,'XTick',0:90:360)
    legend(h, {tuning(exSet).plane(pl).cond.name})
    title(sprintf('%s, %s, plane %d, cell %d', tuning(exSet).subject, ...
        tuning(exSet).date, tuning(exSet).planes(pl), cellID))
    xlim([0 360])
    ylim([mini maxi])
    xlabel('Direction (in degrees)')
    ylabel('\DeltaF/F')
end

%% Collect variables of tuning curves
leg = {'data','null distribution'};

data = load(fullfile(folderResults, 'pupil', ...
    'nullTuning_prefDirSigmaDIFixed.mat'));
null = data.null;

prefDirs = [];
responses = cell(1, length(tuning));

dataset = [];
plane = [];
neuron = [];
minis = []; % resp to non-pref stimulus
maxis = []; % resp to pref stimulus
stimMeans =[];
nullMinis = [];
nullMaxis = [];
nullStimMeans = [];
isSuppr = [];
isGad = [];
draws = cellfun(@size,squeeze(struct2cell(null(1).plane(1).cond(1).cell)),'UniformOutput',false);
draws = max(cell2mat(draws),[],1);
draws = draws(1);
for iExp = 1:length(tuning)
    fprintf('Dataset %d: %s %s exp.: %d\n', iExp, tuning(iExp).subject, ...
        tuning(iExp).date, tuning(iExp).exp);
    respD = [];
    for iPlane = 1:length(tuning(iExp).plane)
        data = tuning(iExp).plane(iPlane);
        if isempty(data.cellIDs)
            continue
        end
        neurons = find(~cellfun(@isempty, {data.cond(1).cell.parameters}));
        
        for j = 1:length(neurons)
            iCell = neurons(j);
            prefDir = NaN;
            resp = NaN(1,12,size(data.cond(1).cell(iCell).responses,2),2);
            means = NaN(1,2);
            prefs = NaN(1,2);
            nonprefs = NaN(1,2);
            nullMeans = NaN(1,2,draws);
            nullPrefs = NaN(1,2,draws);
            nullNonprefs = NaN(1,2,draws);
            for c = 1:2
                resp(1,:,:,c) = data.cond(c).cell(iCell).responses;
                
                pars = data.cond(c).cell(iCell).parameters;
                curve = data.cond(c).cell(iCell).curve;
                if length(pars) == 1 % not tuned
                    means(c) = pars;
                else
                    prefDir = pars(1);
                    means(c) = mean(curve);
                    oris = mod(pars(1) + [0 90 180], 360);
                    sr = gratings.orituneWrappedConditions(pars, oris);
                    prefs(c) = sr(1);
                    if sr(1)-sr(2)>0
                        [nonprefs(c), ind] = min(sr(2:3));
                    else % suppressed neurons
                        [nonprefs(c), ind] = max(sr(2:3));
                    end
                    ind = ind+1;
                end
                
                pars = null(iExp).plane(iPlane).cond(c).cell(iCell).parameters;
                if size(pars,2) == 1 % not tuned
                    nullMeans(1,c,:) = pars;
                else
                    sr = NaN(size(pars,1), 3);
                    curves = NaN(size(pars,1), length(degrees));
                    for p = 1:size(pars,1)
                        oris = mod(pars(p,1) + [0 90 180], 360);
                        sr(p,:) = gratings.orituneWrappedConditions(pars(p,:), oris);
                        curves(p,:) = gratings.orituneWrappedConditions(pars(p,:), degrees);
                    end
                    nullMeans(1,c,:) = mean(curves,2);
                    nullPrefs(1,c,:) = sr(:,1);
                    nullNonprefs(1,c,:) = sr(:,ind);
                end
            end
            
            if ~isempty(corrections)
                a = corrections(iExp).plane(iPlane).a{tuning(iExp).exp} ...
                    (data.cellIDs(iCell));
                b = corrections(iExp).plane(iPlane).b{tuning(iExp).exp} ...
                    (data.cellIDs(iCell));
                resp = doCorrect(a,b,resp);
                nonprefs = doCorrect(a,b,nonprefs);
                prefs = doCorrect(a,b,prefs);
                means = doCorrect(a,b,means);
                nullNonprefs = doCorrect(a,b,nullNonprefs);
                nullPrefs = doCorrect(a,b,nullPrefs);
                nullMeans = doCorrect(a,b,nullMeans);
            end
            
            prefDirs(end+1,1) = prefDir;
            respD(end+1,:,:,:) = resp;
            
            minis(end+1,:) = nonprefs;
            maxis(end+1,:) = prefs;
            stimMeans(end+1,:) = means;
            nullMinis(end+1,:,:) = nullNonprefs;
            nullMaxis(end+1,:,:) = nullPrefs;
            nullStimMeans(end+1,:,:) = nullMeans;
            
            isSuppr(end+1,:) = data.isSuppressed(iCell);
            if isfield(data, 'isGad')
                isGad(end+1,:) = data.isGad(iCell);
            end
            
            dataset(end+1,1) = iExp;
            plane(end+1,1) = iPlane;
            neuron(end+1,1) = data.cellIDs(iCell);
        end
        responses{iExp} = respD;
    end
end

mods = abs(maxis - minis);
nullMods = abs(nullMaxis - nullMinis);

exampleIndices = zeros(size(examples,1),1);
for ex = 1:size(examples,1)
    exampleIndices(ex) = find(dataset==examples(ex,1) & plane==examples(ex,2) & ...
        neuron==examples(ex,3));
end

%% Plot histogram of maximum to differentiate suppressed-by-contrast from
% enhanced-by-contrast neurons
cols = lines(4);
labels = {'small pupil','large pupil'};
isTuned = ~any(isnan(minis),2);
for c = 1 % small pupil (2: large pupil)
    if isempty(isGad)
        types = 1;
        cellTypes = {'all'};
%         cols = [0 0 0];
        bs = 10;
        xLim = [9 10];
        yLim = [3 5];
        zLim = NaN;
    else
        types = [-1 0 1];
        cellTypes = {'exc.','NK','inh.'};
%         cols = [lines(2); 0 0 0];
%         cols = cols([2 3 1],:);
        bs = 1;
        xLim = [4 4];
        yLim = [3 6];
        zLim = 10.46;
    end
    for ty = 1:length(types)
        if isempty(isGad)
            ind = true(size(stimMeans,1),1);
        else
            ind = isGad == types(ty);
        end
    
        mi = minis(:,c);
        ma = maxis(:,c);
        j = isnan(mi);
        mi(j) = stimMeans(j,c);
        ma(j) = stimMeans(j,c);
        j = ma - mi < 0;
        m = mi(j);
        mi(j) = ma(j);
        ma(j) = m;
        b_mi = [-flip(2.^(-3:xLim(1))), 0, 2.^(-3:xLim(2))];
        b_ma = [-flip(2.^(-3:yLim(1))), 0, 2.^(-3:yLim(2))];
        
        figure
        N1 = histcounts(ma(ind&isTuned), b_ma);
        N2 = histcounts(ma(ind&~isTuned), b_ma);
        b = bar(1.5:length(b_ma), [N1',N2'], 'stacked');
        b(1).FaceColor = 'k';
        b(2).FaceColor = 'w';
        hold on
        for ex = 1:size(examples,1)
            j = exampleIndices(ex);
            if ~ind(j)
                continue
            end
            exMax = find(b_ma<ma(j),1,'last');
            if ma(j)>0
                exMax = exMax + rem(log2(ma(j)),1);
            else
                exMax = exMax - rem(log2(-ma(j)),1);
            end
            plot(exMax,10,'v','Color',cols(ex,:))
        end
        xlabel('Maximum')
        ylabel(['#' label])
        title(sprintf('%s - %s (n = %d)', cellTypes{ty}, labels{c},sum(ind)))
        set(gca,'XTick',1:length(b_ma),'XTickLabel',b_ma,'box','off')
        legend(b,{'tuned','not tuned'},'Location','NorthWest')
    end
end

%% Plot tuning modulation measures
% invert sign of responses of suppressed cells
minis(isSuppr==1,:) = -minis(isSuppr==1,:);
maxis(isSuppr==1,:) = -maxis(isSuppr==1,:);
stimMeans(isSuppr==1,:) = -stimMeans(isSuppr==1,:);
nullMinis(isSuppr==1,:,:) = -nullMinis(isSuppr==1,:,:);
nullMaxis(isSuppr==1,:,:) = -nullMaxis(isSuppr==1,:,:);
nullStimMeans(isSuppr==1,:,:) = -nullStimMeans(isSuppr==1,:,:);

ind = all(isnan(maxis),2);
maxis2 = maxis;
maxis2(ind,:) = stimMeans(ind,:);
nullMaxis2 = nullMaxis;
nullMaxis2(ind,:,:) = nullStimMeans(ind,:,:);
measures = {minis, maxis, mods, stimMeans, maxis2};
nullMeasures = {nullMinis, nullMaxis, nullMods, nullStimMeans, nullMaxis2};
% exMeasures = {exMinis, exMaxis, abs(exMaxis - exMinis), exStimMeans};
labels = {'Minimum', 'Maximum', 'Modulation depth', 'Mean', 'Maximum (w/ untuned)'};
modFuns = @(a,b) (b-a)./(abs(a)+abs(b));
funLabels = 'Diff-index (large-small)/(large+small)';
binSizes = 0.1;

cols = lines(4);

% for boutons: get data for neurons for comparison ------------------------
% if strcmp(label, 'boutons')
%     leg = {'data','null distribution','neurons'};
%     % leg = {'data','null distribution','neurons (exc)','neurons (inh&enh)','neurons (inh&sup)'};
%     fR = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\nonVisualEffects\modelGratingResp\';
%     data = load(fullfile(fR, 'pupil', 'tuning_prefDirSigmaDIFixed_lineFit.mat'));
%     tuning_neurons = data.tuning;
%     minis_n = [];
%     maxis_n = [];
%     stimMeans_n =[];
%     isSuppr_n = [];
%     isGad_n = [];
%     for iExp = 1:length(tuning_neurons)
%         fprintf('Dataset %d: %s %s exp.: %d\n', iExp, tuning_neurons(iExp).subject, ...
%             tuning_neurons(iExp).date, tuning_neurons(iExp).exp);
%         for iPlane = 1:length(tuning_neurons(iExp).plane)
%             data = tuning_neurons(iExp).plane(iPlane);
%             if isempty(data.cellIDs)
%                 continue
%             end
%             neurons = find(~cellfun(@isempty, {data.cond(1).cell.parameters}));
%             
%             for j = 1:length(neurons)
%                 iCell = neurons(j);
%                 means = NaN(1,2);
%                 prefs = NaN(1,2);
%                 nonprefs = NaN(1,2);
%                 for c = 1:2
%                     pars = data.cond(c).cell(iCell).parameters;
%                     curve = data.cond(c).cell(iCell).curve;
%                     if length(pars) == 1 % not tuned
%                         means(c) = pars;
%                     else
%                         means(c) = mean(curve);
%                         oris = mod(pars(1) + [0 90 180], 360);
%                         sr = gratings.orituneWrappedConditions(pars, oris);
%                         prefs(c) = sr(1);
%                         if sr(1)-sr(2)>0
%                             [nonprefs(c), ind] = min(sr(2:3));
%                         else % suppressed neurons
%                             [nonprefs(c), ind] = max(sr(2:3));
%                         end
%                         ind = ind+1;
%                     end
%                 end
%                 
%                 minis_n(end+1,:) = nonprefs;
%                 maxis_n(end+1,:) = prefs;
%                 stimMeans_n(end+1,:) = means;
%                 isSuppr_n(end+1,:) = data.isSuppressed(iCell);
%                 isGad_n(end+1,:) = data.isGad(iCell);
%             end
%         end
%     end
%     mods_n = abs(maxis_n - minis_n);
%     ind = all(isnan(maxis_n),2);
%     maxis2_n = maxis_n;
%     maxis2_n(ind,:) = stimMeans_n(ind,:);
%     % invert sign of responses of suppressed cells
%     minis_n(isSuppr_n==1,:) = -minis_n(isSuppr_n==1,:);
%     maxis_n(isSuppr_n==1,:) = -maxis_n(isSuppr_n==1,:);
%     stimMeans_n(isSuppr_n==1,:) = -stimMeans_n(isSuppr_n==1,:);
%     
%     measures_n = {minis_n, maxis_n, mods_n, stimMeans_n, maxis2_n};
% end

% for neurons: get data for boutons for comparison ------------------------
if strcmp(label, 'neurons')
    leg = {'data','null distribution','boutons'};
    % leg = {'data','null distribution','neurons (exc)','neurons (inh&enh)','neurons (inh&sup)'};
    fR = fullfile(folderPC, 'Lab\RESULTS\boutons\nonvisualEffects');
    data = load(fullfile(fR, 'pupil', 'tuning_prefDirSigmaDIFixed_lineFit.mat'));
    tuning_boutons = data.tuning;
    minis_n = [];
    maxis_n = [];
    stimMeans_n =[];
    isSuppr_n = [];
    for iExp = 1:length(tuning_boutons)
        fprintf('Dataset %d: %s %s exp.: %d\n', iExp, tuning_boutons(iExp).subject, ...
            tuning_boutons(iExp).date, tuning_boutons(iExp).exp);
        for iPlane = 1:length(tuning_boutons(iExp).plane)
            data = tuning_boutons(iExp).plane(iPlane);
            if isempty(data.cellIDs)
                continue
            end
            neurons = find(~cellfun(@isempty, {data.cond(1).cell.parameters}));
            
            for j = 1:length(neurons)
                iCell = neurons(j);
                means = NaN(1,2);
                prefs = NaN(1,2);
                nonprefs = NaN(1,2);
                for c = 1:2
                    pars = data.cond(c).cell(iCell).parameters;
                    curve = data.cond(c).cell(iCell).curve;
                    if length(pars) == 1 % not tuned
                        means(c) = pars;
                    else
                        means(c) = mean(curve);
                        oris = mod(pars(1) + [0 90 180], 360);
                        sr = gratings.orituneWrappedConditions(pars, oris);
                        prefs(c) = sr(1);
                        if sr(1)-sr(2)>0
                            [nonprefs(c), ind] = min(sr(2:3));
                        else % suppressed neurons
                            [nonprefs(c), ind] = max(sr(2:3));
                        end
                        ind = ind+1;
                    end
                end
                
                minis_n(end+1,:) = nonprefs;
                maxis_n(end+1,:) = prefs;
                stimMeans_n(end+1,:) = means;
                isSuppr_n(end+1,:) = data.isSuppressed(iCell);
            end
        end
    end
    mods_n = abs(maxis_n - minis_n);
    ind = all(isnan(maxis_n),2);
    maxis2_n = maxis_n;
    maxis2_n(ind,:) = stimMeans_n(ind,:);
    % invert sign of responses of suppressed cells
    minis_n(isSuppr_n==1,:) = -minis_n(isSuppr_n==1,:);
    maxis_n(isSuppr_n==1,:) = -maxis_n(isSuppr_n==1,:);
    stimMeans_n(isSuppr_n==1,:) = -stimMeans_n(isSuppr_n==1,:);
    
    measures_n = {minis_n, maxis_n, mods_n, stimMeans_n, maxis2_n};
end
% -------------------------------------------------------------------------

% Plot histograms and cumulative distributions of difference indices (+
% examples)
limits = [1 .4 .7 .5 .4];
significant = NaN(size(stimMeans,1),length(measures));
groupPValues = NaN(1,length(measures));
diffIndices = NaN(size(stimMeans,1),length(measures));
nullIndices = NaN(size(stimMeans,1),size(nullMeasures{1},3),length(measures));
for m = 1:length(measures)
    diffIndices(:,m) = modFuns(measures{m}(:,1), measures{m}(:,2));
    nullIndices(:,:,m) = modFuns(squeeze(nullMeasures{m}(:,1,:)), ...
        squeeze(nullMeasures{m}(:,2,:)));
    confInt = prctile(nullIndices(:,:,m), [2.5 97.5], 2);
    significant(:,m) = diffIndices(:,m) < confInt(:,1) | ...
        diffIndices(:,m) > confInt(:,2);
    pVals = sum(nullIndices(:,:,m) < diffIndices(:,m), 2) ./ size(nullIndices,2);
    ind = pVals > 0.5;
    pVals(ind) = 1 - pVals(ind);
    pVals = 2 .* pVals;
    pVals(pVals==0) = 1/size(nullIndices,2);
    % Fisher's method to combine p-values
    ind = ~isnan(diffIndices(:,m));
    chi_vals = -2.*log(pVals(ind));
    groupPValues(m) = 1 - chi2cdf(sum(chi_vals),2*length(chi_vals));
end

for m = 1:length(measures)
    % histogram
    figure
    mini = -limits(m); %round(min(real) / binSizes) * binSizes;
    maxi = limits(m); %round(max(real) / binSizes) * binSizes;
    bins = mini:binSizes:maxi;
    edges = [bins-binSizes/2, maxi+binSizes/2];
    n1 = histcounts(diffIndices(significant(:,m)==1,m), edges);
    n2 = histcounts(diffIndices(significant(:,m)==0,m), edges);
    b = bar(bins, [n1',n2'], 'stacked');
    b(1).FaceColor = 'k';
    b(2).FaceColor = 'w';
    hold on
    for ex = 1:length(exampleIndices)
        plot(diffIndices(exampleIndices(ex),m), 10, 'v', 'Color', cols(ex,:))
    end
    xlim(edges([1 end]))
    title(sprintf('%s (n=%d, %d%% signif.)', labels{m}, sum(~isnan(diffIndices(:,m))), ...
        round(sum(significant(:,m))/sum(~isnan(diffIndices(:,m)))*100)))
    xlabel(funLabels)
    ylabel(['#' label])
    ax = gca;
    ax.Box = 'off';
    ax.XTick = [mini 0 maxi];
    
    % cumulative distribution
    figure('Position', [1250 680 560 420])
    hold on
    h = [0 0];
    plot([0 0], [0 1], 'k:')
%     % for boutons: plot date of neurons for comparison --------------------
%     if strcmp(label, 'boutons')
%         real_n = modFuns(measures_n{m}(:,1), measures_n{m}(:,2));
%         ind = ~isnan(real_n);
%         x_n = sort(real_n(ind), 'ascend');
%         y_n = (1:sum(ind)) ./ sum(ind);
%         x_n = [-1; x_n; 1];
%         y_n = [0 y_n 1];
%         h(3) = plot(x_n, y_n, 'Color', [1 1 1].*0.7, 'LineWidth', 2);
%         %     x_n = sort(real_n(ind&isGad_n==-1), 'ascend');
%         %     y_n = (1:sum(ind&isGad_n==-1)) ./ sum(ind&isGad_n==-1);
%         %     x_n = [-1; x_n; 1];
%         %     y_n = [0 y_n 1];
%         %     h(3) = plot(x_n, y_n, 'Color', 'm', 'LineWidth', 2);
%         %     x_n = sort(real_n(ind&isGad_n==1&isSuppr_n==-1), 'ascend');
%         %     y_n = (1:sum(ind&isGad_n==1&isSuppr_n==-1)) ./ sum(ind&isGad_n==1&isSuppr_n==-1);
%         %     x_n = [-1; x_n; 1];
%         %     y_n = [0 y_n 1];
%         %     h(4) = plot(x_n, y_n, 'Color', [0 0 .8], 'LineWidth', 2);
%         %     x_n = sort(real_n(ind&isGad_n==1&isSuppr_n==1), 'ascend');
%         %     y_n = (1:sum(ind&isGad_n==1&isSuppr_n==1)) ./ sum(ind&isGad_n==1&isSuppr_n==1);
%         %     x_n = [-1; x_n; 1];
%         %     y_n = [0 y_n 1];
%         %     h(5) = plot(x_n, y_n, 'Color', [0 .6 0], 'LineWidth', 2);
%     end
%     % ---------------------------------------------------------------------
    % for neurons: plot date of boutons for comparison --------------------
    if strcmp(label, 'neurons')
        real_n = modFuns(measures_n{m}(:,1), measures_n{m}(:,2));
        ind = ~isnan(real_n);
        x_n = sort(real_n(ind), 'ascend');
        y_n = (1:sum(ind)) ./ sum(ind);
        x_n = [-1; x_n; 1];
        y_n = [0 y_n 1];
        h(3) = plot(x_n, y_n, 'Color', [1 1 1].*0.7, 'LineWidth', 2);
    end
    % ---------------------------------------------------------------------

    ind = ~isnan(diffIndices(:,m));
    x = sort(diffIndices(ind,m), 'ascend');
    y = (1:sum(ind)) ./ sum(ind);
    x = [-1; x; 1];
    y = [0 y 1];
    pseudo = nullIndices(:,:,m);
    pseudo(all(isnan(pseudo),2),:) = [];
    
    xNull = sort(pseudo, 1, 'ascend');
    xNull = sort(xNull, 2, 'ascend');
    limNull = prctile(xNull, [2.5 97.5], 2);
    limNull = [[-1 -1]; limNull; [1 1]];
    xNull = nanmedian(xNull, 2);
%     ind = ~isnan(pseudo);
%     xNull = sort(pseudo(ind), 'ascend');
    
    xNull = [min(x); xNull; max(x)];
    yNull = (1:length(xNull)) ./ length(xNull);
    
    fill([limNull(:,1);flip(limNull(:,2))], [yNull, flip(yNull)], ...
        'k', 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.2)
    h(1) = plot(x, y, 'k', 'LineWidth', 2);
    h(2) = plot(xNull, yNull, ':k', 'LineWidth', 2);
    exReal = diffIndices(exampleIndices,m);
    exReal(exReal>edges(end)) = edges(end);
    exReal(exReal<edges(1)) = edges(1);
    exSign = significant(exampleIndices,m);
    [xUni, ind] = unique(x);
    yUni = y(ind);
    heights = interp1(xUni, yUni, exReal);
    for j = 1:length(exampleIndices)
        if exSign(j)==1
            plot(exReal(j), heights(j), 'o', 'MarkerEdgeColor', cols(j,:), ...
                'MarkerFaceColor', cols(j,:), 'LineWidth', 2)
        else
            plot(exReal(j), heights(j), 'o', 'MarkerEdgeColor', cols(j,:), ...
                'LineWidth', 2)
        end
    end
    xlim(edges([1 end]))
    xlabel(funLabels)
    ylabel(['Proportion of ' label])
    title(sprintf('%s (n=%d)', labels{m}, sum(~isnan(diffIndices(:,m)))))
    legend(h, leg, 'Location', 'NorthWest')
    legend('boxoff')
    set(gca, 'XTick', [mini 0 maxi]);
end
figure
hold on
for j = 1:length(exIndices)
    plot(j, 1, 'o', 'MarkerEdgeColor', cols(j,:), ...
        'MarkerFaceColor', cols(j,:), 'MarkerSize', 30)
end

% plot excitatory, inhibitory enhanced and inhibitory suppressed (1. only
% confirmed inhibitory, 2. all suppressed (assuming inhibitory))
for m = 1:length(measures)
    % excitatory
    figure
    mini = -limits(m); %round(min(real) / binSizes) * binSizes;
    maxi = limits(m); %round(max(real) / binSizes) * binSizes;
    bins = mini:binSizes:maxi;
    n1 = hist(diffIndices(significant(:,m)==1 & isGad==-1,m), bins);
    n2 = hist(diffIndices(significant(:,m)==0 & isGad==-1,m), bins);
    b = bar(bins, [n1',n2'], 'stacked');
    b(1).FaceColor = 'r';
    b(2).FaceColor = 'w';
    xlim([mini-.5*binSizes maxi+.5*binSizes])
    n = ~isnan(diffIndices(:,m)) & isGad==-1;
    set(gca, 'XTick', [mini 0 maxi]);
    title(sprintf('%s - excitatory (n=%d, %d%% signif.)', labels{m}, ...
        sum(n), round(sum(significant(:,m) & n) / sum(n)*100)))
    xlabel(funLabels)
    ylabel('#Neurons')
    ax = gca;
    ax.Box = 'off';
    % inhibitory & enhanced
    figure
    n1 = hist(diffIndices(significant(:,m)==1 & isGad==1 & isSuppr==-1,m), bins);
    n2 = hist(diffIndices(significant(:,m)==0 & isGad==1 & isSuppr==-1,m), bins);
    b = bar(bins, [n1',n2'], 'stacked');
    b(1).FaceColor = 'b';
    b(2).FaceColor = 'w';
    xlim([mini-.5*binSizes maxi+.5*binSizes])
    n = ~isnan(diffIndices(:,m)) & isGad==1 & isSuppr==-1;
    set(gca, 'XTick', [mini 0 maxi]);
    title(sprintf('%s - inhibitory & enhanced (n=%d, %d%% signif.)', labels{m}, ...
        sum(n), round(sum(significant(:,m) & n) / sum(n)*100)))
    xlabel(funLabels)
    ylabel('#Neurons')
    ax = gca;
    ax.Box = 'off';
    % inhibitory & suppressed
    figure
    n1 = hist(diffIndices(significant(:,m)==1 & isGad==1 & isSuppr==1,m), bins);
    n2 = hist(diffIndices(significant(:,m)==0 & isGad==1 & isSuppr==1,m), bins);
    b = bar(bins, [n1',n2'], 'stacked');
    b(1).FaceColor = cols(3,:);
    b(2).FaceColor = 'w';
    xlim([mini-.5*binSizes maxi+.5*binSizes])
    n = ~isnan(diffIndices(:,m)) & isGad==1 & isSuppr==1;
    set(gca, 'XTick', [mini 0 maxi]);
    title(sprintf('%s - inhibitory & suppressed (n=%d, %d%% signif.)', labels{m}, ...
        sum(n), round(sum(significant(:,m) & n) / sum(n)*100)))
    xlabel(funLabels)
    ylabel('#Neurons')
    ax = gca;
    ax.Box = 'off';
    % suppressed
    figure
    n1 = hist(diffIndices(significant(:,m)==1 & isSuppr==1,m), bins);
    n2 = hist(diffIndices(significant(:,m)==0 & isSuppr==1,m), bins);
    b = bar(bins, [n1',n2'], 'stacked');
    b(1).FaceColor = cols(4,:);
    b(2).FaceColor = 'w';
    xlim([mini-.5*binSizes maxi+.5*binSizes])
    n = ~isnan(diffIndices(:,m)) & isSuppr==1;
    set(gca, 'XTick', [mini 0 maxi]);
    title(sprintf('%s - suppressed (n=%d, %d%% signif.)', labels{m}, ...
        sum(n), round(sum(significant(:,m) & n) / sum(n)*100)))
    xlabel(funLabels)
    ylabel('#Neurons')
    ax = gca;
    ax.Box = 'off';
    
    figure('Position', [1250 680 560 420])
    hold on
    h = [0 0];
    % excitatory
    ind = ~isnan(diffIndices(:,m)) & isGad==-1;
    xEx = sort(diffIndices(ind,m), 'ascend');
    yEx = (1:sum(ind)) ./ sum(ind);
    xEx = [-1; xEx; 1];
    yEx = [0 yEx 1];
    % inhibitory & enhanced
    ind = ~isnan(diffIndices(:,m)) & isGad==1 & isSuppr==-1;
    xInhEnh = sort(diffIndices(ind,m), 'ascend');
    yInhEnh = (1:sum(ind(:))) ./ sum(ind(:));
    xInhEnh = [-1; xInhEnh; 1];
    yInhEnh = [0 yInhEnh 1];
    % inhibitory & suppressed
    ind = ~isnan(diffIndices(:,m)) & isGad==1 & isSuppr==1;
    xInhSup = sort(diffIndices(ind,m), 'ascend');
    yInhSup = (1:sum(ind(:))) ./ sum(ind(:));
    xInhSup = [-1; xInhSup; 1];
    yInhSup = [0 yInhSup 1];
    % suppressed
    ind = ~isnan(diffIndices(:,m)) & isSuppr==1;
    xSup = sort(diffIndices(ind,m), 'ascend');
    ySup = (1:sum(ind(:))) ./ sum(ind(:));
    xSup = [-1; xSup; 1];
    ySup = [0 ySup 1];
    h(1) = plot(xEx, yEx, 'r', 'LineWidth', 2);
    h(2) = plot(xInhEnh, yInhEnh, 'b', 'LineWidth', 2);
    h(3) = plot(xInhSup, yInhSup, 'Color', cols(3,:), 'LineWidth', 2);
    h(4) = plot(xSup, ySup, 'Color', cols(4,:), 'LineWidth', 2);
    xlim([mini maxi])
    set(gca, 'XTick', [mini 0 maxi]);
    xlabel(funLabels)
    ylabel(['Proportion of ' label])
    title(sprintf('%s', labels{m}))
    legend(h, {'exc','inh+enh','inh+sup','sup'}, 'Location', 'NorthWest')
    legend('boxoff')
    % difference: excitatory versus inhibitory&enhanced
    [~,p] = kstest2(diffIndices(isGad==-1,m), diffIndices(isGad==1 & isSuppr==-1,m));
    fprintf('KS-test, %s, exc. vs inh. & enh.: p = %.2e\n', labels{m}, p)
    % difference: excitatory versus inhibitory&suppressed
    [~,p] = kstest2(diffIndices(isGad==-1,m), diffIndices(isGad==1 & isSuppr==1,m));
    fprintf('KS-test, %s, exc. vs inh. & suppr.: p = %.2e\n', labels{m}, p)
    % difference: inhibitory&enhanced versus inhibitory&suppressed
    [~,p] = kstest2(diffIndices(isGad==1 & isSuppr==-1,m), diffIndices(isGad==1 & isSuppr==1,m));
    fprintf('KS-test, %s, inh. & enh. vs inh. & suppr.: p = %.2e\n', labels{m}, p)
end

% plot enhanced and suppressed (boutons)
for m = 1:length(measures)
    % enhanced
    figure
    mini = -limits(m); %round(min(real) / binSizes) * binSizes;
    maxi = limits(m); %round(max(real) / binSizes) * binSizes;
    bins = mini:binSizes:maxi;
    n1 = hist(diffIndices(significant(:,m)==1 & isSuppr==-1,m), bins);
    n2 = hist(diffIndices(significant(:,m)==0 & isSuppr==-1,m), bins);
    b = bar(bins, [n1',n2'], 'stacked');
    b(1).FaceColor = cols(1,:);
    b(2).FaceColor = 'w';
    xlim([mini-.5*binSizes maxi+.5*binSizes])
    n = ~isnan(diffIndices(:,m)) & isSuppr==-1;
    set(gca, 'XTick', [mini 0 maxi]);
    title(sprintf('%s - enhanced (n=%d, %d%% signif.)', labels{m}, ...
        sum(n), round(sum(significant(:,m) & n) / sum(n)*100)))
    xlabel(funLabels)
    ylabel('#Neurons')
    ax = gca;
    ax.Box = 'off';
    % suppressed
    figure
    n1 = hist(diffIndices(significant(:,m)==1 & isSuppr==1,m), bins);
    n2 = hist(diffIndices(significant(:,m)==0 & isSuppr==1,m), bins);
    b = bar(bins, [n1',n2'], 'stacked');
    b(1).FaceColor = cols(2,:);
    b(2).FaceColor = 'w';
    xlim([mini-.5*binSizes maxi+.5*binSizes])
    n = ~isnan(diffIndices(:,m)) & isSuppr==1;
    set(gca, 'XTick', [mini 0 maxi]);
    title(sprintf('%s - suppressed (n=%d, %d%% signif.)', labels{m}, ...
        sum(n), round(sum(significant(:,m) & n) / sum(n)*100)))
    xlabel(funLabels)
    ylabel('#Neurons')
    ax = gca;
    ax.Box = 'off';
    
    figure('Position', [1250 680 560 420])
    hold on
    h = [0 0];
    % enhanced
    ind = ~isnan(diffIndices(:,m)) & isSuppr==-1;
    xEnh = sort(diffIndices(ind,m), 'ascend');
    yEnh = (1:sum(ind(:))) ./ sum(ind(:));
    xEnh = [-1; xEnh; 1];
    yEnh = [0 yEnh 1];
    % suppressed
    ind = ~isnan(diffIndices(:,m)) & isSuppr==1;
    xSup = sort(diffIndices(ind,m), 'ascend');
    ySup = (1:sum(ind(:))) ./ sum(ind(:));
    xSup = [-1; xSup; 1];
    ySup = [0 ySup 1];
    h(1) = plot(xEnh, yEnh, 'Color', cols(1,:), 'LineWidth', 2);
    h(2) = plot(xSup, ySup, 'Color', cols(2,:), 'LineWidth', 2);
    % neurons
    real_n = modFuns(measures_n{m}(:,1), measures_n{m}(:,2));
    ind = ~isnan(real_n);
    x_n = sort(real_n(ind&isGad_n==-1), 'ascend');
    y_n = (1:sum(ind&isGad_n==-1)) ./ sum(ind&isGad_n==-1);
    x_n = [-1; x_n; 1];
    y_n = [0 y_n 1];
    h(3) = plot(x_n, y_n, 'Color', 'm', 'LineWidth', 2);
    x_n = sort(real_n(ind&isGad_n==1&isSuppr_n==-1), 'ascend');
    y_n = (1:sum(ind&isGad_n==1&isSuppr_n==-1)) ./ sum(ind&isGad_n==1&isSuppr_n==-1);
    x_n = [-1; x_n; 1];
    y_n = [0 y_n 1];
    h(4) = plot(x_n, y_n, 'Color', [0 0 .8], 'LineWidth', 2);
    x_n = sort(real_n(ind&isGad_n==1&isSuppr_n==1), 'ascend');
    y_n = (1:sum(ind&isGad_n==1&isSuppr_n==1)) ./ sum(ind&isGad_n==1&isSuppr_n==1);
    x_n = [-1; x_n; 1];
    y_n = [0 y_n 1];
    h(5) = plot(x_n, y_n, 'Color', [0 .6 0], 'LineWidth', 2);
    
    xlim([mini maxi])
    set(gca, 'XTick', [mini 0 maxi]);
    xlabel(funLabels)
    ylabel(['Proportion of ' label])
    legend(h, {'enh','sup','neurons (exc)','neurons (inh&enh)','neurons (inh&sup)'}, 'Location', 'NorthWest')
    legend('boxoff')
    % difference: enhanced versus suppressed
    [~,p] = kstest2(diffIndices(isSuppr==-1,m), diffIndices(isSuppr==1,m));
    title(sprintf('%s (p = %.2e)', labels{m}, p))
end

% plot tuned versus untuned neurons
isTuned = ~any(isnan(minis),2);
for m = 4:5
    % tuned
    figure
    mini = -limits(m); %round(min(real) / binSizes) * binSizes;
    maxi = limits(m); %round(max(real) / binSizes) * binSizes;
    bins = mini:binSizes:maxi;
    n1 = hist(diffIndices(significant(:,m)==1 & isTuned,m), bins);
    n2 = hist(diffIndices(significant(:,m)==0 & isTuned,m), bins);
    b = bar(bins, [n1',n2'], 'stacked');
    b(1).FaceColor = 'r';
    b(2).FaceColor = 'w';
    xlim([mini-.5*binSizes maxi+.5*binSizes])
    n = isTuned;
    set(gca, 'XTick', [mini 0 maxi]);
    title(sprintf('%s - tuned (n=%d, %d%% signif.)', labels{m}, ...
        sum(n), round(sum(significant(:,m) & n) / sum(n)*100)))
    xlabel(funLabels)
    ylabel('#Neurons')
    ax = gca;
    ax.Box = 'off';
    % untuned
    figure
    n1 = hist(diffIndices(significant(:,m)==1 & ~isTuned,m), bins);
    n2 = hist(diffIndices(significant(:,m)==0 & ~isTuned,m), bins);
    b = bar(bins, [n1',n2'], 'stacked');
    b(1).FaceColor = 'b';
    b(2).FaceColor = 'w';
    xlim([mini-.5*binSizes maxi+.5*binSizes])
    n = ~isnan(diffIndices(:,m)) & ~isTuned;
    set(gca, 'XTick', [mini 0 maxi]);
    title(sprintf('%s - untuned (n=%d, %d%% signif.)', labels{m}, ...
        sum(n), round(sum(significant(:,m) & n) / sum(n)*100)))
    xlabel(funLabels)
    ylabel('#Neurons')
    ax = gca;
    ax.Box = 'off';
    
    figure('Position', [1250 680 560 420])
    hold on
    h = [0 0];
    % tuned
    ind = ~isnan(diffIndices(:,m)) & isTuned;
    xT = sort(diffIndices(ind,m), 'ascend');
    yT = (1:sum(ind)) ./ sum(ind);
    xT = [-1; xT; 1];
    yT = [0 yT 1];
    % untuned
    ind = ~isnan(diffIndices(:,m)) & ~isTuned;
    xU = sort(diffIndices(ind,m), 'ascend');
    yU = (1:sum(ind(:))) ./ sum(ind(:));
    xU = [-1; xU; 1];
    yU = [0 yU 1];
    h(1) = plot(xT, yT, 'r', 'LineWidth', 2);
    h(2) = plot(xU, yU, 'b', 'LineWidth', 2);
    xlim([mini maxi])
    set(gca, 'XTick', [mini 0 maxi]);
    xlabel(funLabels)
    ylabel(['Proportion of ' label])
    title(sprintf('%s', labels{m}))
    legend(h, {'tuned','untuned'}, 'Location', 'NorthWest')
    legend('boxoff')
    % difference: excitatory versus inhibitory&enhanced
    [~,p] = kstest2(diffIndices(isTuned,m), diffIndices(~isTuned,m));
    fprintf('KS-test, %s, tuned vs untuned: p = %.2e\n', labels{m}, p)
end

% Plot scatter plots
% for s = 1:2
m1 = 5;
m2 = 3;
limits(5) = 0.7;
h = [0 0];
%     m1 = (s-1)*2 + 1;
%     m2 = (s-1)*2 + 2;
    real1 = modFuns(measures{m1}(:,1), measures{m1}(:,2));
    real1(real1 > limits(m1)) = limits(m1);
    real1(real1 < -limits(m1)) = -limits(m1);
    real2 = modFuns(measures{m2}(:,1), measures{m2}(:,2));
    real2(real2 > limits(m2)) = limits(m2);
    real2(real2 < -limits(m2)) = -limits(m2);
    figure
    hold on
    
    h(1) = scatter(real1(isSuppr==-1), real2(isSuppr==-1), 'k', 'filled');
    h(2) = scatter(real1(isSuppr==1), real2(isSuppr==1), 'b', 'filled');
    legend(h,{'enh-by-contrast','supp-by-contrast'},'Location','NorthWest')
%     scatter(real1, real2, 20, 'k', 'filled')
    
    alpha(.3)
    hold on
    for j = 1:length(exIndices)
        plot(real1(exIndices(j)), real2(exIndices(j)), 'o', ...
            'MarkerEdgeColor', cols(j,:), 'MarkerFaceColor', cols(j,:), 'LineWidth', 2)
    end
    xlim([-limits(m1) limits(m1)])
    ylim([-limits(m2) limits(m2)])
    real1 = modFuns(measures{m1}(:,1), measures{m1}(:,2));
    real2 = modFuns(measures{m2}(:,1), measures{m2}(:,2));
    ind = ~isnan(real1) & ~isnan(real2);
    [rho,p] = corr(real1(ind),real2(ind));
    title(sprintf('rho = %.3f (p = %.4f, n = %d)', rho, p, sum(~isnan(real1)&~isnan(real2))))
    xlabel(labels{m1})
    ylabel(labels{m2})
    axis square
% end

%% Determine mean responses and DSIs and OSIs
cond = 1; % small or large pupil
shuffles = 1000;
directions = 0:30:330;
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
meanResp = meanResp ./ sum(meanResp,2);
shuffleResp = shuffleResp ./ sum(shuffleResp,2);

% Determine DSIs
vects = sum(dirVectors .* meanResp, 2);
shuffleVects = squeeze(sum(dirVectors .* shuffleResp, 2));
DSIs = abs(vects);
nullDSIs = abs(shuffleVects);
p_DSI = sum(nullDSIs > DSIs,2) ./ shuffles;

% Determine OSIs
vects = sum(oriVectors .* meanResp, 2);
shuffleVects = squeeze(sum(oriVectors .* shuffleResp, 2));
OSIs = abs(vects);
nullOSIs = abs(shuffleVects);
p_OSI = sum(nullOSIs > OSIs,2) ./ shuffles;

% Determine pref directions
vects = sum(dirVectors .* meanResp, 2);
prefDirs_vectAvg = mod(angle(vects) ./ pi .* 180, 360);

%% Plot histogram of preferred directions
% fit first harmonic to distribution of preferred directions
% formula from: "Robust Temporal Coding of Contrast by V1 Neurons for
% Transient But Not for Steady-State Stimuli" by Mechler, ..., Shapley (J
% Neurosci, 1998)
% harmonic = sum(exp((-1i*2*pi/360).*(1:4).*prefDirs(~isnan(prefDirs) & ...
%     (p_DSI<0.05 | p_OSI<0.05))));

cols = lines(4);
binSize = 10;
edges = 0:binSize:360;
bins = edges(1:end-1) + binSize/2;
N1 = histcounts(prefDirs(p_DSI < 0.05), edges);
N2 = histcounts(prefDirs(p_DSI >= 0.05 & p_OSI < 0.05), edges);
figure
% b = bar(bins, [N1; N2], 'stacked');
% b(1).FaceColor = 'k';
% b(2).FaceColor = 'w';
b = bar(bins, [N2; N1], 'stacked');
b(1).FaceColor = 'w';
b(2).FaceColor = 'k';
hold on
Y = fft(N1+N2);
bins10 = bins + (0:9)'./10.*binSize;
bins10 = bins10(:)';
Y4 = zeros(1,length(bins10));
Y4(1) = Y(1)*10;
Y4([1+4 end-4+1]) = Y([1+4 end-4+1]).*10;
plot(bins10, ifft(Y4), 'r', 'LineWidth', 1)
for ex = 1:length(exampleIndices)
    plot(prefDirs(exampleIndices(ex)), 10, 'v', 'Color', cols(ex,:))
end
xlim([0 360])
set(gca, 'XTick', 0:90:360, 'box', 'off')
legend('direction selective','only orientation selective')
xlabel('Preferred direction (deg)')
ylabel(sprintf('#%s', label))
title(sprintf('Preferred directions (fitted with Gaussians) (n=%d)', ...
    sum(N1+N2)))

% Plot one histogram for each dataset
binSize = 22.5;
minNum = 10;
edges = 0:binSize:360;
bins = edges(1:end-1) + binSize/2;
dSets = unique(dataset);
hists = NaN(max(dSets), length(bins));
for d = 1:length(dSets)
    ind = dataset==dSets(d) & (p_DSI<0.05 | p_OSI<0.05);
    if sum(ind) < minNum
        continue
    end
    hists(dSets(d),:) = histcounts(prefDirs(ind), edges);
end
hists = hists ./ sum(hists, 2);
m = nanmean(hists,1);
se = nanstd(hists,0,1)./sqrt(sum(~isnan(hists),1));
figure
hold on
set(gca, 'ColorOrder', jet(max(dSets)))
plot(bins, hists)
fill([bins flip(bins)], [m+se, flip(m-se)], 'k', 'EdgeColor', 'none', ...
    'FaceColor', 'k', 'FaceAlpha', 0.2)
plot(bins, m, 'k', 'LineWidth', 2)
xlim([0 360])
set(gca, 'XTick', 0:90:360, 'box', 'off')
xlabel('Preferred direction')
ylabel(['Proportion of ' label])
title(sprintf('%s (%d datasets)', label, min(sum(~isnan(hists),1))))

%% Plot histograms and scatter of OSIs and DSIs
indPref = ~isnan(prefDirs);
cols = lines(4);
binSize = 0.02;
mx = ceil(max([DSIs; OSIs]) / binSize) * binSize;

% Plot histograms of DSIs and OSIs
edges = 0:binSize:mx;
bins = edges(1:end-1) + binSize/2;
N1 = histcounts(DSIs(p_DSI < 0.05), edges);
N2 = histcounts(DSIs(p_DSI >= 0.05), edges);
figure
b = bar(bins, [N1; N2], 'stacked');
b(1).FaceColor = 'k';
b(2).FaceColor = 'w';
hold on
for ex = 1:length(exampleIndices)
    if p_DSI(exampleIndices(ex)) < 0.05
        plot(DSIs(exampleIndices(ex)), 10, 'v', 'Color', cols(ex,:), ...
            'MarkerFaceColor', cols(ex,:))
    else
        plot(DSIs(exampleIndices(ex)), 10, 'v', 'Color', cols(ex,:))
    end
end
xlim([0 mx])
set(gca, 'box','off','XTick',0:.1:mx)
xlabel('DSI')
ylabel(sprintf('#%s',label))
title(sprintf('DSI (n=%d)', sum(~isnan(DSIs))))

N1 = histcounts(OSIs(p_OSI < 0.05), edges);
N2 = histcounts(OSIs(p_OSI >= 0.05), edges);
figure
b = bar(bins, [N1; N2], 'stacked');
b(1).FaceColor = 'k';
b(2).FaceColor = 'w';
hold on
for ex = 1:length(exampleIndices)
    if p_OSI(exampleIndices(ex)) < 0.05
        plot(OSIs(exampleIndices(ex)), 10, 'v', 'Color', cols(ex,:), ...
            'MarkerFaceColor', cols(ex,:))
    else
        plot(OSIs(exampleIndices(ex)), 10, 'v', 'Color', cols(ex,:))
    end
end
xlim([0 mx])
set(gca, 'box','off','XTick',0:.1:mx)
xlabel('OSI')
ylabel(sprintf('#%s',label))
title(sprintf('OSI (n=%d)', sum(~isnan(OSIs))))

% Plot scatter of DSI versus OSI
colors = [1 0 0; 0 0 1; .5 0 .5];
inds = [p_DSI<.05&p_OSI>=.05, p_DSI>=.05&p_OSI<.05, p_DSI<.05&p_OSI<.05];
h = [0 0 0];
figure
hold on
plot([0 mx], [0 mx], 'k')
for k = 1:3
    h(k) = scatter(DSIs(indPref & inds(:,k)), OSIs(indPref & inds(:,k)), [], colors(k,:), 'filled');
end
for ex = 1:length(exampleIndices)
    scatter(DSIs(exampleIndices(ex)), OSIs(exampleIndices(ex)), [], ...
        cols(ex,:), 'LineWidth', 2)
end
legend(h,'dir-sel','ori-sel','dir&ori')
axis square
axis([0 mx 0 mx])
set(gca, 'box','off','XTick',0:.1:mx,'YTick',0:.1:mx)
xlabel('DSI')
ylabel('OSI')
title(sprintf('%s (n=%d)', label, sum(indPref&(p_DSI<0.05|p_OSI<0.05))))

% Test significance for separation between high OSIs and high DSIs
goodDSIs = DSIs(~isnan(DSIs) & ~isnan(OSIs) & indPref & (p_DSI<.05|p_OSI<.05));
goodOSIs = OSIs(~isnan(DSIs) & ~isnan(OSIs) & indPref & (p_DSI<.05|p_OSI<.05));
permutations = NaN(length(goodDSIs), 10000);
for p = 1:size(permutations,2)
    permutations(:,p) = randperm(size(permutations,1));
end
permOSIs = goodOSIs(permutations);
meanAbsDiff = mean(abs(goodDSIs - goodOSIs));
permAbsDiffs = mean(abs(goodDSIs - permOSIs), 1);
p = sum(permAbsDiffs > meanAbsDiff) / 1000;

% Plot preferred directions
figure
h=histogram(prefDirs_vectAvg(indPref&(p_DSI<0.05|p_OSI<0.05)), 0:10:360, 'FaceColor', 'k');
hold on
for ex = 1:length(exampleIndices)
    plot(prefDirs_vectAvg(exampleIndices(ex)), 10, 'v', 'Color', cols(ex,:))
end
set(gca, 'box', 'off', 'XTick', 0:90:360)
xlim([0 360])
xlabel('Preferred direction')
ylabel(sprintf('#%s',label))
title(sprintf('Preferred directions (vector average) (n=%d)', sum(h.Values)))

%% Plot population direction tuning curve
groups = [true(size(meanResp,1),1), ... % all
    p_DSI<.05|p_OSI<.05, ...            % exclude non-selective
    p_DSI<.05&p_OSI>=.05, ...           % only dir- but NOT ori-selective
    p_DSI>=.05&p_OSI<.05, ...           % only ori- but NOT dir-selective
    p_DSI<.05&p_OSI<.05];               % boht dir- AND ori-selective
groupNames = {'all','DS or OS','DS only','OS only','DS and OS'};
colors = [0 0 0; 0 0 0; 1 0 0; 0 0 1; .5 0 .5];

respScaled = meanResp ./ max(meanResp,[],2);

for g = 1:size(groups,2)
    figure
    hold on
    m = nanmean(respScaled(groups(:,g),[1:end,1]),1);
    sem = nanstd(respScaled(groups(:,g),[1:end,1]),0,1) ./ sqrt(sum(groups(:,g)));
    fill([directions 360 360 flip(directions)], [m+sem flip(m-sem)], ...
        'k', 'EdgeColor', 'none', 'FaceColor', colors(g,:), 'FaceAlpha', .2)
    plot([directions 360], m, 'Color', colors(g,:), 'LineWidth', 2)
    xlim([0 360])
%     ylim([0.4 0.8])
    set(gca, 'box', 'off', 'XTick', 0:90:360)
    xlabel('Direction')
    ylabel('Mean population response')
    title(sprintf('%s, %s (n = %d)', label, groupNames{g}, sum(groups(:,g))))
end

% Plot one histogram for each dataset
dSets = unique(dataset);
minNum = 10;
curves = cell(1, size(groups,2));
for g = 1:size(groups,2)
    curves{g} = NaN(max(dSets), length(directions)+1);
    for d = 1:length(dSets)
        ind = dataset==dSets(d) & groups(:,g);
        if sum(ind) < minNum
            continue
        end
        curves{g}(dSets(d),:) = nanmean(respScaled(ind,[1:end,1]),1);
    end
    
    figure
    hold on
    set(gca, 'ColorOrder', jet(max(dSets)))
    plot([directions 360], curves{g}, 'LineWidth', 2)
    m = nanmean(curves{g},1);
    sem = nanstd(curves{g},0,1) ./ sqrt(sum(~isnan(curves{g}),1));
    fill([directions 360 360 flip(directions)], [m+sem flip(m-sem)], ...
        'k', 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', .2)
    plot([directions 360], m, 'Color', 'k', 'LineWidth', 2)
    xlim([0 360])
%     ylim([0.4 0.8])
    set(gca, 'box', 'off', 'XTick', 0:90:360)
    xlabel('Direction')
    ylabel('Mean population response')
    title(sprintf('%s, %s (%d datasets)', label, groupNames{g}, min(sum(~isnan(curves{g}),1))))
end

% Plot population direction tuning curve on top of population orientation
% tuning curve for each dataset
for d = 1:length(dSets)
    figure
    hold on
    for g = [4 3 5]
        plot([directions 360], curves{g}(dSets(d),:), 'LineWidth', 2)
    end
    legend(groupNames([4 3 5]))
    xlim([0 360])
    %     ylim([0.4 0.8])
    set(gca, 'box', 'off', 'XTick', 0:90:360)
    xlabel('Direction')
    ylabel('Mean population response')
    title(sprintf('%s, datasets %d', label, dSets(d)))
end

%% OLD: Plot traces of examples cell aligned to stimulus
buffer = 2; % in sec (before and after stim period)
cols = [0 0 0; 1 0 0];

for ex = 1:size(examples,1)
    pl = examples(ex,1);
    if ex==1 || pl~=examples(ex-1,1)
        % load meta
        data=load(fullfile(folderROIData, results(exampleSet).subject, ...
            results(exampleSet).date, num2str(results(exampleSet).(expField)), ...
            sprintf('%s_%d_%s_2P_plane%03d_ROI.mat', results(exampleSet).date, ...
            results(exampleSet).(expField), results(exampleSet).subject, ...
            results(exampleSet).planes(pl))));
        meta = data.meta;
        meta.folderTL = strrep(meta.folderTL, 'ioo.ucl.ac.uk', ...
            'cortexlab.net');
        
        [stimTimes, stimSeq, stimMatrix, frameTimes, samplingRate] = ...
            ssLocal.getStimulusResponseInfo(meta);
        nStim = max(stimSeq.seq);
        nRep = sum(stimSeq.seq==1);
        stimOnsetFrames = [0 diff(sum(stimMatrix,1))];
        stimOnsetFrames = find(stimOnsetFrames == 1);
        stimDuration = mean(stimTimes.offset - stimTimes.onset);
        seq = reshape(stimSeq.seq, nStim, []);
        seq = reshape(seq + (0:size(seq,2)-1).*size(seq,1), [], 1);
        sr = 1/median(diff(frameTimes));
        respTime = [-flip(1/sr:1/sr:buffer) 0:1/sr:stimDuration+buffer];
        framesRel = round(respTime.*sr);
        
        SigLen = size(meta.F_final,1);
        PulseLen = length(results(exampleSet).plane(pl).kernelTime);
        nPulses = length(stimTimes.onset);
        xs = repmat(1:nPulses,PulseLen,1);
        ys = (0:PulseLen-1)' + stimOnsetFrames;
        sUseMe = find(ys<=SigLen); % this makes sure we don't write outside the matrix
    end
    
    trace = meta.F_final(:,results(exampleSet).plane(pl) ...
        .cellIDs(examples(ex,2)));
    
    ss = repmat(results(exampleSet).plane(pl).kernelFit(examples(ex,2)) ...
        .kernel, nPulses, 1);
    if ~isempty(ss)
        ShiftedPulses = sparse(ys(sUseMe), xs(sUseMe), ss(sUseMe), ...
            SigLen, nPulses);
        amps = reshape(results(exampleSet).plane(pl).kernelFit(examples(ex,2)) ...
            .alphaEachTrial',[],1);
        amps = amps(seq);
    end
    
    response = NaN(nStim,nRep,length(respTime)); % [stimulus x trial x time]
    for tr = 1:length(stimSeq.seq)
        if isempty(ss)
            r = trace;
        else
            a = amps;
            a(tr) = 0;
            pred = ShiftedPulses * a;
            r = trace - pred;
        end
        fr = stimOnsetFrames(tr)+framesRel;
        fr(fr<1 | fr>length(trace)) = NaN;
        response(stimSeq.seq(tr),ceil(tr/nStim),~isnan(fr)) = ...
            r(fr(~isnan(fr)));
    end
    
    conditions = NaN(size(response,1)-1, size(response,2));
    for c = 1:2
        conditions(~isnan(tuning(exampleSet).plane(pl).cond(c) ...
            .cell(examples(ex,2)).responses)) = c;
    end
    conditions(end+1,:) = NaN;
    
    min0 = min(response(:));
    max0 = max(response(:));
    mini = Inf;
    maxi = -Inf;
    xDist = .5;
%     yDist = .1 * (max0 - min0);
    traceDur = stimDuration + 2*buffer;
    
    figure('Position',[3 640 1915 420])
    hold on
    y0 = 0;
    marks = 0;
%     marks = linspace(0, maxi, 5);
%     marks(end) = [];
    for c = 1:2
        x0 = 0;
        for st = 1:size(response,1)
            if c == 1
                fill([0 stimDuration stimDuration 0] + x0, ...
                    [min0 min0 max0 max0] + y0, 'k', 'FaceColor', cols(c,:), ...
                    'FaceAlpha', 0.1, 'EdgeColor', 'none')
            end
            ind = find(conditions(st,:) == c);
%             for tr = 1:length(ind)
%                 plot(responseTime + x0, ...
%                     squeeze(response(st,ind(tr),:)) + y0, ...
%                     'Color', cols(c,:).*.3 + [1 1 1].*.7)
%             end
            tr = squeeze(nanmean(response(st,ind,:),2));
            mini = min([mini; tr]);
            maxi = max([maxi; tr]);
            plot(respTime + x0, tr + y0, 'Color', cols(c,:), ...
                'LineWidth', 2)
            plot(repmat(respTime([1 end])',1,4) + x0, ...
                repmat(marks,2,1) + y0, 'k:')
            x0 = x0 + traceDur + xDist;
        end
%         y0 = y0 + mini - maxi - yDist;
    end
    rng = maxi - mini;
    maxi = maxi + .01*rng;
    mini = mini - 0.01*rng;
    axis tight
    ylim([mini maxi])
    set(gca, 'XTick', [0 stimDuration])
    xlabel('Stimuli')
    ylabel('\DeltaF/F')
    title(sprintf('Plane %d, cell %d', tuning(exampleSet).planes(examples(ex,1)), ...
        examples(ex,2)))
end

%% OLD: Plot tuning curves of selected examples
lin = {'-','-'};
cols = {'k', 'r'};

% examples
for n = 1:size(examples,1)
    pl = examples(ex,1);
    cellInd = find(results(exampleSet).plane(pl).cellIDs == examples(ex,2));
    a = 1;
    b = 0;
    if ~isempty(corrections)
        a = corrections(exampleSet).plane(pl).a{tuning(exampleSet).exp}(examples(ex,2));
        b = corrections(exampleSet).plane(pl).b{tuning(exampleSet).exp}(examples(ex,2));
    end
    
    figure
    hold on
    h = [0 0];
    m = cell(1,2);
    s = cell(1,2);
    for c = 1:2
        h(c) = plot(degrees, tuning(exampleSet).plane(examples(n,1)).cond(c) ...
            .cell(examples(n,2)).curve ./ normPar, lin{c}, ...
            'Color', cols{c}, 'LineWidth',2);
        m{c} = nanmean(tuning(exampleSet).plane(examples(n,1)).cond(c) ...
            .cell(examples(n,2)).responses,2) ./ normPar;
        s{c} = nanstd(tuning(exampleSet).plane(examples(n,1)).cond(c) ...
            .cell(examples(n,2)).responses ./ normPar,0,2) ./ ...
            sqrt(sum(~isnan(tuning(exampleSet).plane(examples(n,1)).cond(c) ...
            .cell(examples(n,2)).responses),2));
        errorbar([tuning(exampleSet).plane(examples(n,1)).cond(c) ...
            .cell(examples(n,2)).directions; 360], m{c}([1:end 1]), ...
            s{c}([1:end 1]), 'o', ...
            'Color', cols{c}, 'CapSize', 2, 'MarkerFaceColor', cols{c})
    end
    mini = min([0; cat(1, m{:}) - cat(1, s{:})]);
    maxi = max([0; cat(1, m{:}) + cat(1, s{:})]);
    rng = maxi-mini;
    mini = mini-.05*rng;
    maxi = maxi+.05*rng;
    plot([0 360], [0 0], 'k:', 'LineWidth', 2)
    set(gca,'XTick',0:90:360)
%     legend(h, {tuning(k).plane(examples(n,1)).cond.name})
    title(sprintf('Plane %d, cell %d', tuning(exampleSet).planes(examples(n,1)), ...
        examples(n,2)))
    xlim([0 360])
    ylim([mini maxi])
    xlabel('Direction (in degrees)')
    ylabel('\DeltaF/F')
    
    figure('Position', [1250 680 560 420])
    hold on
%     if tuning(exampleSet).plane(examples(n,1)).isSuppressed(examples(n,2)) == 1
%         for c = 1:2
%             m{c} = -m{c};
%         end
%         tmp = maxi;
%         maxi = -mini;
%         mini = -tmp;
%     end
    plot(m{1}, m{2}, 'ko', 'MarkerFaceColor', 'k')
    intercept = tuning(exampleSet).plane(examples(n,1)).lineFit(examples(n,2)) ...
        .intercept / normPar;
    slope = tuning(exampleSet).plane(examples(n,1)).lineFit(examples(n,2)) ...
        .slope;
    if tuning(exampleSet).plane(examples(n,1)).isSuppressed(examples(n,2)) == 1
        intercept = -intercept;
    end
    plot([mini maxi],[mini maxi].*slope+intercept, ...
        'Color','k','LineWidth',2)
    plot([mini maxi],[mini maxi],'k:')
    axis([mini maxi mini maxi])
    axis square
    xlabel(tuning(exampleSet).plane(examples(n,1)).cond(1).name)
    ylabel(tuning(exampleSet).plane(examples(n,1)).cond(2).name)
    title(sprintf('Plane %d, cell %d', tuning(exampleSet).planes(examples(n,1)), ...
        examples(n,2)))
    set(gca,'box','off')
end

%% OLD: Plot maximum (vs minimum)
% Plot maximum for each dataset separately
isTuned = ~any(isnan(minis),2);
if isempty(isGad)
    types = 1;
    cellTypes = {'all'};
    cols = [0 0 0];
    bs = 10;
    xLim = [9 10];
    yLim = [9 11];
else
    types = [-1 0 1];
    cellTypes = {'exc.','NK','inh.'};
    cols = [lines(2); 0 0 0];
    cols = cols([2 3 1],:);
    bs = 1;
    xLim = [4 4];
    yLim = [3 6];
    zLim = 10.46;
end
for ty = 1:length(types)
    if isempty(isGad)
        ind = true(size(stimMeans,1),1);
    else
        ind = isGad == types(ty);
    end
    
    mi = minis(:,c);
    ma = maxis(:,c);
    j = isnan(mi);
    mi(j) = stimMeans(j,c);
    ma(j) = stimMeans(j,c);
%     j = ma - mi < 0;
%     m = mi(j);
%     mi(j) = ma(j);
%     ma(j) = m;
    b_mi = [-flip(2.^(-3:xLim(1))), 0, 2.^(-3:xLim(2))];
    b_ma = [-flip(2.^(-3:yLim(1))), 0, 2.^(-3:yLim(2))];
    
    for k = unique(dataset)'
        figure
        N1 = histcounts(ma(ind&isTuned&dataset==k), b_ma);
        N2 = histcounts(ma(ind&~isTuned&dataset==k), b_ma);
%         [~,~,b_ex] = histcounts(ma(exIndices), b_ma);
        b = bar(1.5:length(b_ma), [N1',N2'], 'stacked');
        b(1).FaceColor = 'k';
        b(2).FaceColor = 'w';
        hold on
%         for ex = 1:length(exIndices)
%             plot(b_ex(ex)+0.5,10,'v','Color',cols(ex,:))
%         end
        xlabel('Maximum')
        ylabel(['#' label])
        title(sprintf('%s %s: %s (n = %d)', tuning(k).subject, ...
            tuning(k).date, cellTypes{ty}, sum(ind&dataset==k)), 'interpreter', 'none')
        set(gca,'XTick',1:length(b_ma),'XTickLabel',b_ma,'box','off')
        legend(b,{'tuned','not tuned'})
    end
end

% Plot maximum during small versus large pupil for excitatory and
% inhibitory neurons
ma = maxis;
j = isnan(maxis(:,1));
ma(j,:) = stimMeans(j,:);
ma_log = log2(abs(ma));
ma_log(ma_log<-4) = -4;
ma_axis = NaN(size(ma));
ind = ma<0;
ma_axis(ind) = -ma_log(ind)+4;
ma_axis(ind & ma_axis<1) = 1;
ind = ma>0;
ma_axis(ind) = ma_log(ind)+12;
ma_axis(ind & ma_axis>18) = 18;
for ty = 1:length(types)
    ind = isGad == types(ty);
    
    figure
    hold on
    j = ma(:,1)<0;
    plot(ma_axis(ind&j,1), ma_axis(ind&j,2), 'k.')
    plot(ma_axis(ind&~j,1), ma_axis(ind&~j,2), '.', 'Color', [1 1 1].*.7)
    plot([1 18],[1 18],'k:')
    for ex = 1:length(exIndices)
        plot(ma_axis(exIndices(ex),1), ma_axis(exIndices(ex),2), 'o', ...
            'Color', cols(ex,:))
    end
    axis square
    xlim([1 length(b_ma)])
    ylim([1 length(b_ma)])
    set(gca,'XTick',1:length(b_ma),'XTickLabel',b_ma, ...
        'YTick',1:length(b_ma),'YTickLabel',b_ma)
    xlabel('Maximum (small pupil)')
    ylabel('Maximum (large pupil)')
    title(cellTypes{ty})
end

%% OLD: Plot probablity density for preferred directions
binSize = 100;

pd = prefDirs;
pd(isnan(pd)) = [];
n = length(pd);
pd = sort(pd);
pd = [pd(end-binSize:end)-360; pd; pd(1:binSize)+360];
rngs = [pd(1:end-binSize+1), pd(binSize:end)];
y = binSize/n ./ diff(rngs,1,2);
x = mean(rngs, 2);

figure
plot(x, y, 'k', 'LineWidth', 1)
xlim([0 360])
set(gca, 'XTick', 0:90:360, 'box', 'off')
xlabel('Preferred direction (deg)')
ylabel('Probability density')
title(sprintf('%s (n=%d)', label, sum(~isnan(prefDirs))))

