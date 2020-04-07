%% Define data
exSet = 1;
exUnits = [4 7 95];
exCols = lines(length(exUnits));

%% Load data
folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\Electrophys_SC\nonVisualEffects\';

data = load(fullfile(folderResults, 'pupil', ...
    'tuning_prefDirSigmaDIFixed_lineFit.mat'));
tuning = data.tuning;
degrees = data.x;
data = load(fullfile(folderResults, 'pupil', ...
    'nullTuning_prefDirSigmaDIFixed_behaviour.mat'));
nullBeh = data.null;
data = load(fullfile(folderResults, 'pupil', ...
    'nullTuning_prefDirSigmaDIFixed_laser.mat'));
nullLaser = data.null;

%% Extract relevant data
minis = cell(1,4); % response to non-preferred stimulus (only tuned neurons)
maxis = cell(1,4); % response to preferred stimulus (only tuned neurons)
stimMeans = cell(1,4);
blank = cell(1,4);
nullBehMinis = cell(1,4);
nullBehMaxis = cell(1,4);
nullBehStimMeans = cell(1,4);
nullLaserMinis = cell(1,4);
nullLaserMaxis = cell(1,4);
nullLaserStimMeans = cell(1,4);
isSuppr = [];
exIDs = NaN(1, length(exUnits));

for iExp = 1:length(tuning)
    fprintf('Dataset %d: %s %s exp.: %d\n', iExp, tuning(iExp).subject, ...
        tuning(iExp).date, tuning(iExp).exp);
    if isempty(tuning(iExp).cond)
        continue
    end
    for iCell = 1:length(tuning(iExp).cellIDs)
        if isnan(tuning(iExp).isTuned(iCell))
            continue
        end
        for c = 1:4
            pars = tuning(iExp).cond(c).cell(iCell).parameters;
            curve = tuning(iExp).cond(c).cell(iCell).curve;
            if length(pars) == 1
                m = pars;
                mi = NaN;
                ma = NaN;
            else
                m = mean(curve);
                oris = mod(pars(1) + [0 90 180], 360);
                sr = gratings.orituneWrappedConditions(pars, oris);
                ma = sr(1);
                if sr(1)-sr(2)>0
                    [mi, ind] = min(sr(2:3));
                else % suppressed
                    [mi, ind] = max(sr(2:3));
                end
                ind = ind+1;
            end
            minis{c}(end+1,1) = mi;
            maxis{c}(end+1,1) = ma;
            stimMeans{c}(end+1,1) = m;
            blank{c}(end+1,1) = nanmean(tuning(iExp).cond(c).cell(iCell).blankResponses);
            
            pars = nullBeh(iExp).cond(c).cell(iCell).parameters;
            if size(pars,2) == 1
                m = pars;
                mi = NaN(size(m));
                ma = NaN(size(m));
            else
                sr = NaN(size(pars,1), 3);
                curves = NaN(size(pars,1), length(degrees));
                for p = 1:size(pars,1)
                    oris = mod(pars(p,1) + [0 90 180], 360);
                    sr(p,:) = gratings.orituneWrappedConditions(pars(p,:), oris);
                    curves(p,:) = gratings.orituneWrappedConditions(pars(p,:), degrees);
                end
                m = mean(curves,2);
                ma = sr(:,1);
                mi = sr(:,ind);
            end
            nullBehMinis{c}(end+1,:) = mi;
            nullBehMaxis{c}(end+1,:) = ma;
            nullBehStimMeans{c}(end+1,:) = m;
            
            pars = nullLaser(iExp).cond(c).cell(iCell).parameters;
            if size(pars,2) == 1
                m = pars;
                mi = NaN(size(m));
                ma = NaN(size(m));
            else
                sr = NaN(size(pars,1), 3);
                curves = NaN(size(pars,1), length(degrees));
                for p = 1:size(pars,1)
                    oris = mod(pars(p,1) + [0 90 180], 360);
                    sr(p,:) = gratings.orituneWrappedConditions(pars(p,:), oris);
                    curves(p,:) = gratings.orituneWrappedConditions(pars(p,:), degrees);
                end
                m = mean(curves,2);
                ma = sr(:,1);
                mi = sr(:,ind);
            end
            nullLaserMinis{c}(end+1,:) = mi;
            nullLaserMaxis{c}(end+1,:) = ma;
            nullLaserStimMeans{c}(end+1,:) = m;
        end
        if iExp == exSet && ismember(iCell,exUnits)
            exIDs(ismember(exUnits,iCell)) = length(minis{1});
        end
    end
    isSuppr = [isSuppr; tuning(iExp).isSuppressed(~isnan(tuning(iExp).isTuned))];
end

mods = cell(1,4);
nullBehMods = cell(1,4);
nullLaserMods = cell(1,4);
maxis2 = maxis; % response to preferred stimulus, includes untuned neurons
maxis3 = maxis; % includes untuned neurons and is the maximum response (rather than response at pref. direction)
nullBehMaxis2 = nullBehMaxis;
nullLaserMaxis2 = nullLaserMaxis;
nullLaserMaxis3 = nullLaserMaxis;
isTuned = ~isnan(maxis{1});
j = isSuppr==1;
for c = 1:4
    mods{c} = abs(maxis{c} - minis{c});
    nullBehMods{c} = abs(nullBehMaxis{c} - nullBehMinis{c});
    nullLaserMods{c} = abs(nullLaserMaxis{c} - nullLaserMinis{c});
    maxis2{c}(~isTuned) = stimMeans{c}(~isTuned);
    maxis3{c}(~isTuned) = stimMeans{c}(~isTuned);
    maxis3{c}(j & isTuned) = minis{c}(j & isTuned);
    nullBehMaxis2{c}(~isTuned,:) = nullBehStimMeans{c}(~isTuned,:);
    nullLaserMaxis2{c}(~isTuned,:) = nullLaserStimMeans{c}(~isTuned,:);
    nullLaserMaxis3{c}(~isTuned,:) = nullLaserStimMeans{c}(~isTuned,:);
    nullLaserMaxis3{c}(j & isTuned,:) = nullLaserMinis{c}(j & isTuned,:);
    
%     m = minis{c}(j);
%     minis{c}(j) = maxis{c}(j);
%     maxis{c}(j) = m;    
end

%% Make plots (13.03.2020)

measures = {minis, maxis, mods, stimMeans, maxis2};
nullBehMeasures = {nullBehMinis, nullBehMaxis, nullBehMods, ...
    nullBehStimMeans, nullBehMaxis2};
nullLaserMeasures = {nullLaserMinis, nullLaserMaxis, nullLaserMods, ...
    nullLaserStimMeans, nullLaserMaxis2};
labels = {'Minimum', 'Maximum', 'Modulation', 'Mean', 'Maximum (w/ untuned)'};
modFuns = @(a,b) (b-a)./(abs(a)+abs(b));
funLabels = 'Diff-index (large-small)/(large+small)';
binSizes = 0.1;
mini = -1;%-0.95;
maxi = 1;%0.95;
bins = mini : binSizes : maxi;

validUnits = all(~isnan(cat(2,maxis3{:})),2) & ...
    ~isnan(modFuns(measures{5}{1}, measures{5}{2})) & ...
    ~isnan(modFuns(measures{5}{3}, measures{5}{4}));

% Plot max. resp. and blank resp. when laser off versus laser on
measure = maxis3;
nullMeas = nullLaserMaxis3;
conds = {'Small','Large'};
for c = 1:2
    mini = min([measure{c}(validUnits); measure{c+2}(validUnits)]);
    maxi = max([measure{c}(validUnits); measure{c+2}(validUnits)]);
    diffs = modFuns(measure{c}, measure{c+2});
    pseudoDiffsLaser = modFuns(nullMeas{c}, nullMeas{c+2});
    confInt = prctile(pseudoDiffsLaser, [2.5 97.5], 2);
    sgnf = diffs<confInt(:,1) | diffs>confInt(:,2);
    h = [0 0];
    figure
    hold on
    h(1) = plot(measure{c}(validUnits & ~sgnf), measure{c+2}(validUnits & ~sgnf), ...
        'o','Color',[1 1 1].*.8,'MarkerFaceColor',[1 1 1].*.8);
    h(2) = plot(measure{c}(validUnits & sgnf), measure{c+2}(validUnits & sgnf), ...
        'ok','MarkerFaceColor','k');
    plot([mini maxi],[mini maxi],'r')
    axis image
    set(gca,'box','off')
    xlabel('Maximum (laser off)')
    ylabel('Maximum (laser on)')
    p = signrank(measure{c}(validUnits), measure{c+2}(validUnits));
    title(sprintf('%s pupil: (off - on)/off = %.2f (p = %.3f, n = %d)', ...
        conds{c}, nanmedian((measure{c}(validUnits) - ...
        measure{c+2}(validUnits)) ./ measure{c}(validUnits)), p, ...
        sum(validUnits)))
    legend(h,{'p>=0.05','p<0.05'},'Location','NorthWest')
    p = signrank(modFuns(measure{c}(validUnits), measure{c+2}(validUnits)))
    
    mini = min([blank{c}(validUnits); blank{c+2}(validUnits)]);
    maxi = max([blank{c}(validUnits); blank{c+2}(validUnits)]);
    figure
    hold on
    plot(blank{c}(validUnits),blank{c+2}(validUnits),'ok','MarkerFaceColor','k')
    plot([mini maxi],[mini maxi],'r')
    axis image
    set(gca,'box','off')
    xlabel('Blank (laser off)')
    ylabel('Blank (laser on)')
    p = signrank(blank{c}(validUnits), blank{c+2}(validUnits));
    title(sprintf('%s pupil: (off - on)/off = %.2f (p = %.3f, n = %d)', ...
        conds{c}, nanmedian((blank{c}(validUnits) - ...
        blank{c+2}(validUnits)) ./ blank{c}(validUnits)), p, ...
        sum(validUnits)))
    
    ratiosMax = modFuns(measure{c+2}(validUnits), measure{c}(validUnits));
    ratiosBlank = modFuns(blank{c+2}(validUnits), blank{c}(validUnits));
%     ratiosMax = (measure{c}(validUnits) - measure{c+2}(validUnits)) ...
%         ./ measure{c}(validUnits);
%     ratiosBlank = (blank{c}(validUnits) - blank{c+2}(validUnits)) ...
%         ./ blank{c}(validUnits);
    figure
    hold on
    plot(ratiosBlank,ratiosMax,'ok','MarkerFaceColor','k')
    axis square
    axis([-1 1 -1 1])
    set(gca,'box','off')
    xlabel('DI Blank')
    ylabel('DI Maximum')
    ind = ~isnan(ratiosMax) & ~isnan(ratiosBlank) & ~isinf(ratiosBlank);
    [rho,p] = corr(ratiosMax(ind),ratiosBlank(ind));
    title(sprintf('%s pupil; rho = %.2f (p = %.3f, n = %d)', ...
        conds{c}, rho, p, sum(ind)))
end

% Plot scatters and histograms: diff. index - laser off vs. on
for m = 5 % 1:length(measures)
    validUnits = ...
        ~isnan(modFuns(measures{m}{1}, measures{m}{2})) & ...
        ~isnan(modFuns(measures{m}{3}, measures{m}{4}));

    diffOff = modFuns(measures{m}{1}, measures{m}{2});
    diffOn = modFuns(measures{m}{3}, measures{m}{4});
    pseudoBehOff = modFuns(nullBehMeasures{m}{1}, nullBehMeasures{m}{2});
    pseudoBehOn = modFuns(nullBehMeasures{m}{3}, nullBehMeasures{m}{4});
    pseudoLaserOff = modFuns(nullLaserMeasures{m}{1}, nullLaserMeasures{m}{2});
    pseudoLaserOn = modFuns(nullLaserMeasures{m}{3}, nullLaserMeasures{m}{4});
    
    % for suppressed units, change sign of DIs of preferred responses so that
    % a more negative response with arousal results in a positive DI
    diffOff(isSuppr==1) = -diffOff(isSuppr==1);
    diffOn(isSuppr==1) = -diffOn(isSuppr==1);
    pseudoBehOff(isSuppr==1) = -pseudoBehOff(isSuppr==1);
    pseudoBehOn(isSuppr==1) = -pseudoBehOn(isSuppr==1);
    pseudoLaserOff(isSuppr==1) = -pseudoLaserOff(isSuppr==1);
    pseudoLaserOn(isSuppr==1) = -pseudoLaserOn(isSuppr==1);
    
    confIntBehOff = prctile(pseudoBehOff, [2.5 97.5], 2);
    signBeh = diffOff < confIntBehOff(:,1) | diffOff > confIntBehOff(:,2);
    pVals = sum(pseudoBehOff < diffOff, 2) ./ size(pseudoBehOff,2);
    ind = pVals > 0.5;
    pVals(ind) = 1 - pVals(ind);
    pVals = 2 .* pVals;
    pVals(pVals==0) = 1/size(pseudoBehOff,2);
    % Fisher's method to combine p-values
    chi_vals = -2.*log(pVals(validUnits));
    group_pval_behOff = 1 - chi2cdf(sum(chi_vals),2*length(chi_vals));
    
    confIntBehOn = prctile(pseudoBehOn, [2.5 97.5], 2);
    signBehOn = diffOn < confIntBehOn(:,1) | diffOn > confIntBehOn(:,2);
    pVals = sum(pseudoBehOn < diffOn, 2) ./ size(pseudoBehOn,2);
    ind = pVals > 0.5;
    pVals(ind) = 1 - pVals(ind);
    pVals = 2 .* pVals;
    pVals(pVals==0) = 1/size(pseudoBehOn,2);
    % Fisher's method to combine p-values
    chi_vals = -2.*log(pVals(validUnits));
    group_pval_behOn = 1 - chi2cdf(sum(chi_vals),2*length(chi_vals));
    
    % changed on 13.03.2020
    diffs = diffOn - diffOff;
    pseudoDiffsLaser = pseudoLaserOn - pseudoLaserOff;
    confIntLaser = prctile(pseudoDiffsLaser, [2.5 97.5], 2);
    signLaser = diffs < confIntLaser(:,1) | diffs > confIntLaser(:,2);
    pVals = sum(pseudoDiffsLaser < diffs, 2) ./ size(pseudoDiffsLaser,2);
    ind = pVals > 0.5;
    pVals(ind) = 1 - pVals(ind);
    pVals = 2 .* pVals;
    pVals(pVals==0) = 1/size(pseudoDiffsLaser,2);
    % Fisher's method to combine p-values
    chi_vals = -2.*log(pVals(validUnits));
    group_pval_laser = 1 - chi2cdf(sum(chi_vals),2*length(chi_vals));

    % scatter for laser off versus laser on condition
    n = [0 0 0 0];
    h = [0 0 0 0];
    figure('Position', [635 685 660 420])
    hold on
    % (1) significant for laser (large dots)
    % (1a) not significant for pupil when laser off (black dots)
    ind = signLaser & ~signBeh;
    n(1) = sum(ind);
    h(1) = plot(diffOff(ind), diffOn(ind), 'o', 'MarkerSize', 7, ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w');
    % (1b) significant for pupil when laser off (red dots)
    ind = signLaser & signBeh;
    n(2) = sum(ind);
    h(2) = plot(diffOff(ind), diffOn(ind), 'o', 'MarkerSize', 7, ...
        'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'w');
    % (2) not significant for laser (small dots)
    % (2a) not significant for pupil when laser off (black dots)
    ind = ~signLaser & ~signBeh;
    n(3) = sum(ind);
    h(3) = plot(diffOff(ind), diffOn(ind), 'o', 'MarkerSize', 5, ...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k');
    % (2b) significant for pupil when laser off (red dots)
    ind = ~signLaser & signBeh;
    n(4) = sum(ind);
    h(4) = plot(diffOff(ind), diffOn(ind), 'o', 'MarkerSize', 5, ...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r');
    for ex = 1:length(exUnits)
        plot(diffOff(exIDs(ex)), diffOn(exIDs(ex)), 'o', ...
            'MarkerSize', 10, 'MarkerEdgeColor', exCols(ex,:), ...
            'LineWidth', 2, 'MarkerFaceColor', 'none')
    end
    plot([-1 1], [-1 1], 'k:')
    legend(h, ['p_{laser}<.05, p_{pupil}\geq.05 (n=' num2str(n(1)) ')'], ...
        ['p_{laser}<.05, p_{pupil}<.05 (n=' num2str(n(2)) ')'], ...
        ['p_{laser}\geq.05, p_{pupil}\geq.05 (n=' num2str(n(3)) ')'], ...
        ['p_{laser}\geq.05, p_{pupil}<.05 (n=' num2str(n(4)) ')'], ...
        'Location', 'northeastoutside')
    axis([-1 1 -1 1])
    axis square
    xlabel('laser off')
    ylabel('laser on')
    title(['\Delta' sprintf('%s: off - on = %.3f (p = %.3f, n = %d)', labels{m}, ...
        nanmedian(diffs), ...
        signrank(diffOff,diffOn), sum(~isnan(diffs)))])
    set(gca, 'box', 'off', 'XTick', -1:1, 'YTick', -1:1)
    
    % histogram for laser off condition
    figure
    n1 = hist(diffOff(signBeh & validUnits), bins);
    n2 = hist(diffOff(~signBeh & validUnits), bins);
    b = bar(bins, [n1',n2'], 'stacked');
    b(1).FaceColor = 'k';
    b(2).FaceColor = 'w';
    xlim([-1.05 1.05])
    title(sprintf('%s (laser off) (n = %d)', labels{m}, sum(validUnits)))
    xlabel('Diff index')
    ylabel('#Neurons')
    set(gca, 'box', 'off', 'XTick', -1:1)
    % cumulative distribution for laser off condition
    x = sort(diffOff(validUnits), 'ascend');
    y = (1:sum(validUnits)) ./ sum(validUnits);
    x = [-1; x; 1];
    y = [0 y 1];
    ind = validUnits & ~all(isnan(pseudoBehOff),2);
    nulls = pseudoBehOff(ind,:);
    xNull = sort(nulls, 1, 'ascend');
    xNull = sort(xNull, 2, 'ascend');
    limNull = prctile(xNull, [2.5 97.5], 2);
    limNull = [[-1 -1]; limNull; [1 1]];
    xNull = median(xNull, 2);
    xNull = [-1; xNull; 1];
    yNull = (1:length(xNull))' ./ length(xNull);
    figure
    hold on
    h = [0 0];
    fill([limNull(:,1);flip(limNull(:,2))], [yNull; flip(yNull)], ...
        'k', 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.2)
    h(1) = plot(x, y, 'k', 'LineWidth', 1);
    h(2) = plot(xNull, yNull, ':k', 'LineWidth', 1);
    exReal = diffOff(exIDs);
    [xUni, ind] = unique(x);
    yUni = y(ind);
    heights = interp1(xUni, yUni, exReal);
    for ex = 1:length(exIDs)
        plot(exReal(ex), heights(ex), 'o', 'MarkerEdgeColor', exCols(ex,:), ...
            'MarkerFaceColor', exCols(ex,:), 'LineWidth', 2)
    end
    title(['\Delta ' labels{m} ' (laser off)'])
    xlabel('Diff index')
    ylabel('#Neurons')
    set(gca, 'box', 'off', 'XTick', -1:1)
    % histogram for laser on condition
    figure
    n1 = hist(diffOn(signBehOn & validUnits), bins);
    n2 = hist(diffOn(~signBehOn & validUnits), bins);
    b = bar(bins, [n1',n2'], 'stacked');
    b(1).FaceColor = 'k';
    b(2).FaceColor = 'w';
    xlim([-1.05 1.05])
    title(sprintf('%s (laser on) (n = %d)', labels{m}, sum(validUnits)))
    xlabel('Diff index')
    ylabel('#Neurons')
    set(gca, 'box', 'off', 'XTick', -1:1)
    % cumulative distribution for laser on condition
    x = sort(diffOn(validUnits), 'ascend');
    y = (1:sum(validUnits)) ./ sum(validUnits);
    x = [-1; x; 1];
    y = [0 y 1];
    ind = validUnits & ~all(isnan(pseudoBehOn),2);
    nulls = pseudoBehOn(ind,:);
    xNull = sort(nulls, 1, 'ascend');
    xNull = sort(xNull, 2, 'ascend');
    limNull = prctile(xNull, [2.5 97.5], 2);
    limNull = [[-1 -1]; limNull; [1 1]];
    xNull = median(xNull, 2);
    xNull = [-1; xNull; 1];
    yNull = (1:length(xNull))' ./ length(xNull);
    figure
    hold on
    h = [0 0];
    fill([limNull(:,1);flip(limNull(:,2))], [yNull; flip(yNull)], ...
        'k', 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.2)
    h(1) = plot(x, y, 'k', 'LineWidth', 1);
    h(2) = plot(xNull, yNull, ':k', 'LineWidth', 1);
    exReal = diffOn(exIDs);
    [xUni, ind] = unique(x);
    yUni = y(ind);
    heights = interp1(xUni, yUni, exReal);
    for ex = 1:length(exIDs)
        plot(exReal(ex), heights(ex), 'o', 'MarkerEdgeColor', exCols(ex,:), ...
            'MarkerFaceColor', exCols(ex,:), 'LineWidth', 2)
    end
    title(['\Delta ' labels{m} ' (laser on)'])
    xlabel('Diff index')
    ylabel('#Neurons')
    set(gca, 'box', 'off', 'XTick', -1:1)
    
    % histogram for laser differences (diffOn - diffOff)
    figure
    n1 = hist(diffs(signLaser & validUnits), bins);
    n2 = hist(diffs(~signLaser & validUnits), bins);
    b = bar(bins, [n1',n2'], 'stacked');
    b(1).FaceColor = 'k';
    b(2).FaceColor = 'w';
    hold on
    plot(nanmedian(diffs(validUnits)), max(n1+n2), 'kv', ...
        'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'none')
    xlim([-1.05 1.05])
    title(sprintf('%s (n = %d)', labels{m}, sum(validUnits)))
    xlabel('laser on - laser off')
    ylabel('#Neurons')
    set(gca, 'box', 'off', 'XTick', -1:1)
    % cumulative distribution for laser differences
    x = sort(diffs(validUnits), 'ascend');
    y = (1:sum(validUnits)) ./ sum(validUnits);
    x = [-1; x; 1];
    y = [0 y 1];
    ind = validUnits;
    nulls = abs(pseudoLaserOn(ind,:)) - abs(pseudoLaserOff(ind,:));
    xNull = sort(nulls, 1, 'ascend');
    xNull = sort(xNull, 2, 'ascend');
    limNull = prctile(xNull, [2.5 97.5], 2);
    limNull = [[-1 -1]; limNull; [1 1]];
    xNull = nanmedian(xNull, 2);
    xNull = [-1; xNull; 1];
    yNull = (1:length(xNull))' ./ length(xNull);
    figure
    hold on
    h = [0 0];
    fill([limNull(:,1);flip(limNull(:,2))], [yNull; flip(yNull)], ...
        'k', 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.2)
    plot(x, y, 'k', 'LineWidth', 1);
    plot(xNull, yNull, ':k', 'LineWidth', 1);
    exReal = diffs(exIDs);
    [xUni, ind] = unique(x);
    yUni = y(ind);
    heights = interp1(xUni, yUni, exReal);
    for ex = 1:length(exIDs)
        plot(exReal(ex), heights(ex), 'o', 'MarkerEdgeColor', exCols(ex,:), ...
            'MarkerFaceColor', exCols(ex,:), 'LineWidth', 2)
    end
    title(['\Delta ' labels{m}])
    xlabel('laser on - laser off')
    ylabel('#Neurons')
    set(gca, 'box', 'off', 'XTick', -1:1)
end

%% Make plots (OLD)

measures = {minis, maxis, mods, stimMeans, maxis2};
nullBehMeasures = {nullBehMinis, nullBehMaxis, nullBehMods, ...
    nullBehStimMeans, nullBehMaxis2};
nullLaserMeasures = {nullLaserMinis, nullLaserMaxis, nullLaserMods, ...
    nullLaserStimMeans, nullLaserMaxis2};
labels = {'Minimum', 'Maximum', 'Modulation', 'Mean', 'Maximum (w/ untuned)'};
modFuns = @(a,b) (b-a)./(abs(a)+abs(b));
funLabels = 'Diff-index (large-small)/(large+small)';
binSizes = 0.1;
mini = -1;%-0.95;
maxi = 1;%0.95;
bins = mini : binSizes : maxi;

validUnits = all(~isnan(cat(2,maxis3{:})),2) & ...
    ~isnan(modFuns(measures{5}{1}, measures{5}{2})) & ...
    ~isnan(modFuns(measures{5}{3}, measures{5}{4}));

% Plot max. resp. and blank resp. when laser off versus laser on
measure = maxis3;
nullMeas = nullLaserMaxis3;
conds = {'Small','Large'};
for c = 1:2
    mini = min([measure{c}(validUnits); measure{c+2}(validUnits)]);
    maxi = max([measure{c}(validUnits); measure{c+2}(validUnits)]);
    diffs = modFuns(measure{c}, measure{c+2});
    pseudoDiffsLaser = modFuns(nullMeas{c}, nullMeas{c+2});
    confInt = prctile(pseudoDiffsLaser, [2.5 97.5], 2);
    sgnf = diffs<confInt(:,1) | diffs>confInt(:,2);
    h = [0 0];
    figure
    hold on
    h(1) = plot(measure{c}(validUnits & ~sgnf), measure{c+2}(validUnits & ~sgnf), ...
        'o','Color',[1 1 1].*.8,'MarkerFaceColor',[1 1 1].*.8);
    h(2) = plot(measure{c}(validUnits & sgnf), measure{c+2}(validUnits & sgnf), ...
        'ok','MarkerFaceColor','k');
    plot([mini maxi],[mini maxi],'r')
    axis image
    set(gca,'box','off')
    xlabel('Maximum (laser off)')
    ylabel('Maximum (laser on)')
    p = signrank(measure{c}(validUnits), measure{c+2}(validUnits));
    title(sprintf('%s pupil: (off - on)/off = %.2f (p = %.3f, n = %d)', ...
        conds{c}, nanmedian((measure{c}(validUnits) - ...
        measure{c+2}(validUnits)) ./ measure{c}(validUnits)), p, ...
        sum(validUnits)))
    legend(h,{'p>=0.05','p<0.05'},'Location','NorthWest')
    
    mini = min([blank{c}(validUnits); blank{c+2}(validUnits)]);
    maxi = max([blank{c}(validUnits); blank{c+2}(validUnits)]);
    figure
    hold on
    plot(blank{c}(validUnits),blank{c+2}(validUnits),'ok','MarkerFaceColor','k')
    plot([mini maxi],[mini maxi],'r')
    axis image
    set(gca,'box','off')
    xlabel('Blank (laser off)')
    ylabel('Blank (laser on)')
    p = signrank(blank{c}(validUnits), blank{c+2}(validUnits));
    title(sprintf('%s pupil: (off - on)/off = %.2f (p = %.3f, n = %d)', ...
        conds{c}, nanmedian((blank{c}(validUnits) - ...
        blank{c+2}(validUnits)) ./ blank{c}(validUnits)), p, ...
        sum(validUnits)))
    
    ratiosMax = modFuns(measure{c+2}(validUnits), measure{c}(validUnits));
    ratiosBlank = modFuns(blank{c+2}(validUnits), blank{c}(validUnits));
%     ratiosMax = (measure{c}(validUnits) - measure{c+2}(validUnits)) ...
%         ./ measure{c}(validUnits);
%     ratiosBlank = (blank{c}(validUnits) - blank{c+2}(validUnits)) ...
%         ./ blank{c}(validUnits);
    figure
    hold on
    plot(ratiosBlank,ratiosMax,'ok','MarkerFaceColor','k')
    axis square
    axis([-1 1 -1 1])
    set(gca,'box','off')
    xlabel('DI Blank')
    ylabel('DI Maximum')
    ind = ~isnan(ratiosMax) & ~isnan(ratiosBlank) & ~isinf(ratiosBlank);
    [rho,p] = corr(ratiosMax(ind),ratiosBlank(ind));
    title(sprintf('%s pupil; rho = %.2f (p = %.3f, n = %d)', ...
        conds{c}, rho, p, sum(ind)))
end

% Plot scatters and histograms: diff. index - laser off vs. on
for m = 5 % 1:length(measures)
    if m == 3
        validUnits = ...
            ~isnan(modFuns(measures{m}{1}, measures{m}{2})) & ...
            ~isnan(modFuns(measures{m}{3}, measures{m}{4}));
    end
    diffOff = modFuns(measures{m}{1}, measures{m}{2});
    diffOn = modFuns(measures{m}{3}, measures{m}{4});
    pseudoBehOff = modFuns(nullBehMeasures{m}{1}, nullBehMeasures{m}{2});
    pseudoBehOn = modFuns(nullBehMeasures{m}{3}, nullBehMeasures{m}{4});
    pseudoLaserOff = modFuns(nullLaserMeasures{m}{1}, nullLaserMeasures{m}{2});
    pseudoLaserOn = modFuns(nullLaserMeasures{m}{3}, nullLaserMeasures{m}{4});
    
    confIntBehOff = prctile(pseudoBehOff, [2.5 97.5], 2);
    signBeh = diffOff < confIntBehOff(:,1) | diffOff > confIntBehOff(:,2);
    pVals = sum(pseudoBehOff < diffOff, 2) ./ size(pseudoBehOff,2);
    ind = pVals > 0.5;
    pVals(ind) = 1 - pVals(ind);
    pVals = 2 .* pVals;
    pVals(pVals==0) = 1/size(pseudoBehOff,2);
    % Fisher's method to combine p-values
    chi_vals = -2.*log(pVals(validUnits));
    group_pval_behOff = 1 - chi2cdf(sum(chi_vals),2*length(chi_vals));
    
    confIntBehOn = prctile(pseudoBehOn, [2.5 97.5], 2);
    signBehOn = diffOn < confIntBehOn(:,1) | diffOn > confIntBehOn(:,2);
    pVals = sum(pseudoBehOn < diffOn, 2) ./ size(pseudoBehOn,2);
    ind = pVals > 0.5;
    pVals(ind) = 1 - pVals(ind);
    pVals = 2 .* pVals;
    pVals(pVals==0) = 1/size(pseudoBehOn,2);
    % Fisher's method to combine p-values
    chi_vals = -2.*log(pVals(validUnits));
    group_pval_behOn = 1 - chi2cdf(sum(chi_vals),2*length(chi_vals));
    
    % changed on 26.03.2019
    diffs = abs(diffOn) - abs(diffOff);
    pseudoDiffsLaser = abs(pseudoLaserOn) - abs(pseudoLaserOff);
    confIntLaser = prctile(pseudoDiffsLaser, [2.5 97.5], 2);
    signLaser = diffs < confIntLaser(:,1) | diffs > confIntLaser(:,2);
    pVals = sum(pseudoDiffsLaser < diffs, 2) ./ size(pseudoDiffsLaser,2);
    ind = pVals > 0.5;
    pVals(ind) = 1 - pVals(ind);
    pVals = 2 .* pVals;
    pVals(pVals==0) = 1/size(pseudoDiffsLaser,2);
    % Fisher's method to combine p-values
    chi_vals = -2.*log(pVals(validUnits));
    group_pval_laser = 1 - chi2cdf(sum(chi_vals),2*length(chi_vals));
%     diffs = diffOff - diffOn;
%     pseudoDiffsLaser = pseudoLaserOff - pseudoLaserOn;
%     confIntLaser = prctile(pseudoDiffsLaser, [2.5 97.5], 2);
%     signLaser = diffs < confIntLaser(:,1) | diffs > confIntLaser(:,2);
    
%     validUnits = ~isnan(diffOff) & ~isnan(diffOn);

%     n = [0 0 0 0];
%     h = [0 0 0 0];
%     figure('Position', [635 685 660 420])
%     hold on
%     fill([-1 1 1 -1],[-1 1 -1 1], 'k', 'EdgeColor', 'none', ...
%         'FaceColor', 'k', 'FaceAlpha', .1)
%     % (1) significant for laser (large dots)
%     % (1a) not significant for pupil when laser off (black dots)
%     ind = signLaser & ~signBeh;
%     n(1) = sum(ind);
%     h(1) = plot(diffOff(ind), diffOn(ind), 'o', 'MarkerSize', 7, ...
%         'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w');
%     % (1b) significant for pupil when laser off (red dots)
%     ind = signLaser & signBeh;
%     n(2) = sum(ind);
%     h(2) = plot(diffOff(ind), diffOn(ind), 'o', 'MarkerSize', 7, ...
%         'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'w');
%     % (2) not significant for laser (small dots)
%     % (2a) not significant for pupil when laser off (black dots)
%     ind = ~signLaser & ~signBeh;
%     n(3) = sum(ind);
%     h(3) = plot(diffOff(ind), diffOn(ind), 'o', 'MarkerSize', 5, ...
%         'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k');
%     % (2b) significant for pupil when laser off (red dots)
%     ind = ~signLaser & signBeh;
%     n(4) = sum(ind);
%     h(4) = plot(diffOff(ind), diffOn(ind), 'o', 'MarkerSize', 5, ...
%         'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r');
%     for ex = 1:length(exUnits)
%         plot(diffOff(exIDs(ex)), diffOn(exIDs(ex)), 'o', ...
%             'MarkerSize', 10, 'MarkerEdgeColor', exCols(ex,:), ...
%             'LineWidth', 2, 'MarkerFaceColor', 'none')
%     end
%     plot([-1 1], [-1 1], 'k:')
%     legend(h, ['p_{laser}<.05, p_{pupil}\geq.05 (n=' num2str(n(1)) ')'], ...
%         ['p_{laser}<.05, p_{pupil}<.05 (n=' num2str(n(2)) ')'], ...
%         ['p_{laser}\geq.05, p_{pupil}\geq.05 (n=' num2str(n(3)) ')'], ...
%         ['p_{laser}\geq.05, p_{pupil}<.05 (n=' num2str(n(4)) ')'], ...
%         'Location', 'northeastoutside')
%     axis([-1 1 -1 1])
%     axis square
%     xlabel('laser off')
%     ylabel('laser on')
%     title(['\Delta' sprintf('%s: |off| - |on| = %.3f (p = %.3f, n = %d)', labels{m}, ...
%         nanmedian(diffs), ...
%         signrank(abs(diffOff),abs(diffOn)), sum(~isnan(diffs)))])
%     set(gca, 'box', 'off', 'XTick', -1:1, 'YTick', -1:1)
    
    % histogram for laser off condition
    figure
    n1 = hist(diffOff(signBeh & validUnits), bins);
    n2 = hist(diffOff(~signBeh & validUnits), bins);
    b = bar(bins, [n1',n2'], 'stacked');
    b(1).FaceColor = 'k';
    b(2).FaceColor = 'w';
    xlim([-1.05 1.05])
    title(sprintf('%s (laser off) (n = %d)', labels{m}, sum(validUnits)))
    xlabel('Diff index')
    ylabel('#Neurons')
    set(gca, 'box', 'off', 'XTick', -1:1)
    % cumulative distribution for laser off condition
    x = sort(diffOff(validUnits), 'ascend');
    y = (1:sum(validUnits)) ./ sum(validUnits);
    x = [-1; x; 1];
    y = [0 y 1];
    ind = validUnits & ~all(isnan(pseudoBehOff),2);
    nulls = pseudoBehOff(ind,:);
    xNull = sort(nulls, 1, 'ascend');
    xNull = sort(xNull, 2, 'ascend');
    limNull = prctile(xNull, [2.5 97.5], 2);
    limNull = [[-1 -1]; limNull; [1 1]];
    xNull = median(xNull, 2);
    xNull = [-1; xNull; 1];
    yNull = (1:length(xNull))' ./ length(xNull);
    figure
    hold on
    h = [0 0];
    fill([limNull(:,1);flip(limNull(:,2))], [yNull; flip(yNull)], ...
        'k', 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.2)
    h(1) = plot(x, y, 'k', 'LineWidth', 1);
    h(2) = plot(xNull, yNull, ':k', 'LineWidth', 1);
    exReal = diffOff(exIDs);
    [xUni, ind] = unique(x);
    yUni = y(ind);
    heights = interp1(xUni, yUni, exReal);
    for ex = 1:length(exIDs)
        plot(exReal(ex), heights(ex), 'o', 'MarkerEdgeColor', exCols(ex,:), ...
            'MarkerFaceColor', exCols(ex,:), 'LineWidth', 2)
    end
    title(['\Delta ' labels{m} ' (laser off)'])
    xlabel('Diff index')
    ylabel('#Neurons')
    set(gca, 'box', 'off', 'XTick', -1:1)
    % histogram for laser on condition
    figure
    n1 = hist(diffOn(signBehOn & validUnits), bins);
    n2 = hist(diffOn(~signBehOn & validUnits), bins);
    b = bar(bins, [n1',n2'], 'stacked');
    b(1).FaceColor = 'k';
    b(2).FaceColor = 'w';
    xlim([-1.05 1.05])
    title(sprintf('%s (laser on) (n = %d)', labels{m}, sum(validUnits)))
    xlabel('Diff index')
    ylabel('#Neurons')
    set(gca, 'box', 'off', 'XTick', -1:1)
    % cumulative distribution for laser on condition
    x = sort(diffOn(validUnits), 'ascend');
    y = (1:sum(validUnits)) ./ sum(validUnits);
    x = [-1; x; 1];
    y = [0 y 1];
    ind = validUnits & ~all(isnan(pseudoBehOn),2);
    nulls = pseudoBehOn(ind,:);
    xNull = sort(nulls, 1, 'ascend');
    xNull = sort(xNull, 2, 'ascend');
    limNull = prctile(xNull, [2.5 97.5], 2);
    limNull = [[-1 -1]; limNull; [1 1]];
    xNull = median(xNull, 2);
    xNull = [-1; xNull; 1];
    yNull = (1:length(xNull))' ./ length(xNull);
    figure
    hold on
    h = [0 0];
    fill([limNull(:,1);flip(limNull(:,2))], [yNull; flip(yNull)], ...
        'k', 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.2)
    h(1) = plot(x, y, 'k', 'LineWidth', 1);
    h(2) = plot(xNull, yNull, ':k', 'LineWidth', 1);
    exReal = diffOn(exIDs);
    [xUni, ind] = unique(x);
    yUni = y(ind);
    heights = interp1(xUni, yUni, exReal);
    for ex = 1:length(exIDs)
        plot(exReal(ex), heights(ex), 'o', 'MarkerEdgeColor', exCols(ex,:), ...
            'MarkerFaceColor', exCols(ex,:), 'LineWidth', 2)
    end
    title(['\Delta ' labels{m} ' (laser on)'])
    xlabel('Diff index')
    ylabel('#Neurons')
    set(gca, 'box', 'off', 'XTick', -1:1)
    
    % histogram for laser differences (abs(diffOff) - abs(diffOn))
    figure
    n1 = hist(diffs(signLaser & validUnits), bins);
    n2 = hist(diffs(~signLaser & validUnits), bins);
    b = bar(bins, [n1',n2'], 'stacked');
    b(1).FaceColor = 'k';
    b(2).FaceColor = 'w';
    hold on
    plot(nanmedian(diffs(validUnits)), max(n1+n2), 'kv', ...
        'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'none')
    p = signtest(abs(diffOn(validUnits)),abs(diffOff(validUnits)));
    xlim([-1.05 1.05])
    title(sprintf('%s (n = %d)', labels{m}, sum(validUnits)))
    xlabel('|laser on| - |laser off|')
    ylabel('#Neurons')
    set(gca, 'box', 'off', 'XTick', -1:1)
    % cumulative distribution for laser differences
    x = sort(diffs(validUnits), 'ascend');
    y = (1:sum(validUnits)) ./ sum(validUnits);
    x = [-1; x; 1];
    y = [0 y 1];
    ind = validUnits;
    nulls = abs(pseudoLaserOn(ind,:)) - abs(pseudoLaserOff(ind,:));
    xNull = sort(nulls, 1, 'ascend');
    xNull = sort(xNull, 2, 'ascend');
    limNull = prctile(xNull, [2.5 97.5], 2);
    limNull = [[-1 -1]; limNull; [1 1]];
    xNull = nanmedian(xNull, 2);
    xNull = [-1; xNull; 1];
    yNull = (1:length(xNull))' ./ length(xNull);
    figure
    hold on
    h = [0 0];
    fill([limNull(:,1);flip(limNull(:,2))], [yNull; flip(yNull)], ...
        'k', 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.2)
    plot(x, y, 'k', 'LineWidth', 1);
    plot(xNull, yNull, ':k', 'LineWidth', 1);
    exReal = diffs(exIDs);
    [xUni, ind] = unique(x);
    yUni = y(ind);
    heights = interp1(xUni, yUni, exReal);
    for ex = 1:length(exIDs)
        plot(exReal(ex), heights(ex), 'o', 'MarkerEdgeColor', exCols(ex,:), ...
            'MarkerFaceColor', exCols(ex,:), 'LineWidth', 2)
    end
    title(['\Delta ' labels{m}])
    xlabel('|laser on| - |laser off|')
    ylabel('#Neurons')
    set(gca, 'box', 'off', 'XTick', -1:1)
end