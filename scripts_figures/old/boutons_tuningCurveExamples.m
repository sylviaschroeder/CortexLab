%% Define data

k = 4; % SS069, 2016-10-21, 1
examples = [1 5; 1 9; 1 15; 1 120]; % [plane, cellID]

folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\boutons\nonVisualEffects';

data = load(fullfile(folderResults, 'pupil', ...
    'tuning_prefDirSigmaDIFixed_lineFit.mat'));
tuning = data.tuning;
degrees = data.x;

data = load(fullfile(folderResults, 'kernelFit', 'results.mat'));
results = data.results;

%% Plot traces of 1st example cell aligned to stimulus
buffer = 2; % in sec (before and after stim period)
for ex = 1:size(examples,1)
    response = permute(squeeze(results(k).plane(examples(ex,1)).responses( ...
        examples(ex,2),:,:,:)), [3 2 1]); % [stimuli x trials x time]
    stimDur = results(k).plane(examples(ex,1)).stimDuration;
    responseTime = results(k).plane(examples(ex,1)).responseTime;
    traceInds = responseTime >= -buffer & responseTime <= stimDur + buffer;
    responseTime = responseTime(traceInds);
    response = response(:,:,traceInds);
    conditions = NaN(size(response,1), size(response,2));
    for c = 1:2
        conditions(~isnan(tuning(k).plane(examples(ex,1)).cond(c) ...
            .cell(examples(ex,2)).responses)) = c;
    end
    
    mini = min(response(:));
    maxi = max(response(:));
    rng = maxi - mini;
    % maxi = maxi + .01*rng;
    % mini = mini - 0.01*rng;
    xDist = .5;
    yDist = .1 * rng;
    traceDur = stimDur + 2*buffer;
    
    figure('Position',[3 640 1915 420])
    hold on
    y0 = 0;
    cols = [0 0 0; 1 0 0];
    marks = linspace(0, maxi, 5);
    marks(end) = [];
    for c = 1:2
        x0 = 0;
        for st = 1:size(response,1)
            fill([0 stimDur stimDur 0] + x0, ...
                [mini mini maxi maxi] + y0, 'k', 'FaceColor', cols(c,:), ...
                'FaceAlpha', 0.1, 'EdgeColor', 'none')
            ind = find(conditions(st,:) == c);
            for tr = 1:length(ind)
                plot(responseTime + x0, ...
                    squeeze(response(st,ind(tr),:)) + y0, ...
                    'Color', cols(c,:).*.3 + [1 1 1].*.7)
            end
            plot(responseTime + x0, ...
                squeeze(mean(response(st,ind,:),2)) + y0, 'Color', cols(c,:), ...
                'LineWidth', 2)
            plot(repmat([-buffer stimDur+buffer]',1,4) + x0, ...
                repmat(marks,2,1) + y0, 'k:')
            x0 = x0 + traceDur + xDist;
        end
        y0 = y0 + mini - maxi - yDist;
    end
    axis tight
    set(gca, 'XTick', [0 stimDur], 'YTick', [0 maxi])
    xlabel('Stimuli')
    ylabel('\DeltaF/F')
    title(sprintf('Plane %d, cell %d', tuning(k).planes(examples(ex,1)), ...
        examples(ex,2)))
end

%% Plot tuning curves of selected examples
lin = {'-','-'};
cols = {'k', 'r'};

% examples
for n = 1:size(examples,1)
%     normPar = abs(sum(tuning(k).plane(examples(n,1)).cond(1) ...
%         .cell(examples(n,2)).parameters([2 4])));
    normPar = 1;
    
    figure
    hold on
    h = [0 0];
    m = cell(1,2);
    s = cell(1,2);
    for c = 1:2
        h(c) = plot(degrees, tuning(k).plane(examples(n,1)).cond(c) ...
            .cell(examples(n,2)).curve ./ normPar, lin{c}, ...
            'Color', cols{c}, 'LineWidth',2);
        m{c} = nanmean(tuning(k).plane(examples(n,1)).cond(c) ...
            .cell(examples(n,2)).responses,2) ./ normPar;
        s{c} = nanstd(tuning(k).plane(examples(n,1)).cond(c) ...
            .cell(examples(n,2)).responses ./ normPar,0,2) ./ ...
            sqrt(sum(~isnan(tuning(k).plane(examples(n,1)).cond(c) ...
            .cell(examples(n,2)).responses),2));
        errorbar(tuning(k).plane(examples(n,1)).cond(c) ...
            .cell(examples(n,2)).directions, m{c}, s{c}, 'o', ...
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
    title(sprintf('Plane %d, cell %d', tuning(k).planes(examples(n,1)), ...
        examples(n,2)))
    xlim([-10 370])
    ylim([mini maxi])
    xlabel('Direction (in degrees)')
    ylabel('\DeltaF/F')
    
    figure('Position', [1250 680 560 420])
    hold on
    if tuning(k).plane(examples(n,1)).isSuppressed(examples(n,2)) == 1
        for c = 1:2
            m{c} = -m{c};
        end
        tmp = maxi;
        maxi = -mini;
        mini = -tmp;
    end
    plot(m{1}, m{2}, 'ko', 'MarkerFaceColor', 'k')
    intercept = tuning(k).plane(examples(n,1)).lineFit(examples(n,2)) ...
        .intercept / normPar;
    slope = tuning(k).plane(examples(n,1)).lineFit(examples(n,2)) ...
        .slope;
    if tuning(k).plane(examples(n,1)).isSuppressed(examples(n,2)) == 1
        intercept = -intercept;
    end
    plot([mini maxi],[mini maxi].*slope+intercept, ...
        'Color','k','LineWidth',2)
    plot([mini maxi],[mini maxi],'k:')
    axis([mini maxi mini maxi])
    axis square
    xlabel(tuning(k).plane(examples(n,1)).cond(1).name)
    ylabel(tuning(k).plane(examples(n,1)).cond(2).name)
    title(sprintf('Plane %d, cell %d', tuning(k).planes(examples(n,1)), ...
        examples(n,2)))
    set(gca,'box','off')
end

%% Plot difference index for minimum, maximum, ... of tuning curve for population
data = load(fullfile(folderResults, 'pupil', ...
    'nullTuning_prefDirSigmaDIFixed.mat'));
null = data.null;

minis = [];
maxis = [];
stimMeans =[];
nullMinis = [];
nullMaxis = [];
nullStimMeans = [];
exIndices = NaN(size(examples,1),1);
isSuppr = [];
for iExp = 1:length(tuning)
    fprintf('Dataset %d: %s %s exp.: %d\n', iExp, tuning(iExp).subject, ...
        tuning(iExp).date, tuning(iExp).exp);
    for iPlane = 1:length(tuning(iExp).plane)
        data = tuning(iExp).plane(iPlane);
        if isempty(data.cellIDs)
            continue
        end
        neurons = find(~cellfun(@isempty, {data.cond(1).cell.parameters}));
        draws = size(null(iExp).plane(iPlane).cond(1).cell(neurons(1)).parameters,1);
        
        for j = 1:length(neurons)
            iCell = neurons(j);
            means = NaN(1,2);
            prefs = NaN(1,2);
            nonprefs = NaN(1,2);
            nullMeans = NaN(1,2,draws);
            nullPrefs = NaN(1,2,draws);
            nullNonprefs = NaN(1,2,draws);
            for c = 1:2
                pars = data.cond(c).cell(iCell).parameters;
                curve = data.cond(c).cell(iCell).curve;
                if length(pars) == 1 % not tuned
                    means(c) = pars;
                else
                    means(c) = mean(curve);
                    oris = mod(pars(1) + [0 90 180], 360);
                    sr = orituneWrappedConditions(pars, oris);
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
                        sr(p,:) = orituneWrappedConditions(pars(p,:), oris);
                        curves(p,:) = orituneWrappedConditions(pars(p,:), degrees);
                    end
                    nullMeans(1,c,:) = mean(curves,2);
                    nullPrefs(1,c,:) = sr(:,1);
                    nullNonprefs(1,c,:) = sr(:,ind);
                end
            end
            
            minis(end+1,:) = nonprefs;
            maxis(end+1,:) = prefs;
            stimMeans(end+1,:) = means;
            nullMinis(end+1,:,:) = nullNonprefs;
            nullMaxis(end+1,:,:) = nullPrefs;
            nullStimMeans(end+1,:,:) = nullMeans;
            if k == iExp
                ex = find(ismember(examples(:,1),iPlane) & ...
                    ismember(examples(:,2),iCell));
                if ~isempty(ex)
                    exIndices(ex) = size(minis,1);
                end
            end
            
            isSuppr(end+1,:) = data.isSuppressed(iCell);
        end
    end
end

mods = abs(maxis - minis);
nullMods = abs(nullMaxis - nullMinis);

% invert sign of responses of suppressed cells
minis(isSuppr==1,:) = -minis(isSuppr==1,:);
maxis(isSuppr==1,:) = -maxis(isSuppr==1,:);
stimMeans(isSuppr==1,:) = -stimMeans(isSuppr==1,:);
nullMinis(isSuppr==1,:,:) = -nullMinis(isSuppr==1,:,:);
nullMaxis(isSuppr==1,:,:) = -nullMaxis(isSuppr==1,:,:);
nullStimMeans(isSuppr==1,:,:) = -nullStimMeans(isSuppr==1,:,:);

measures = {minis, maxis, mods, stimMeans};
nullMeasures = {nullMinis, nullMaxis, nullMods, nullStimMeans};
labels = {'Minimum', 'Maximum', 'Modulation depth', 'Mean'};
modFuns = @(a,b) (b-a)./(abs(a)+abs(b));
funLabels = 'Diff-index (large-small)/(large+small)';
binSizes = 0.05;

% Plot histograms and cumulative distributions of difference indices (+
% examples)
cols = lines(length(exIndices));
for m = 1:length(measures)
    figure
    real = modFuns(measures{m}(:,1), measures{m}(:,2));
    pseudo = modFuns(squeeze(nullMeasures{m}(:,1,:)), ...
        squeeze(nullMeasures{m}(:,2,:)));
    confInt = prctile(pseudo, [2.5 97.5], 2);
    isSignf = real < confInt(:,1) | real > confInt(:,2);
%     isSignf = true(size(real));
    mini = -.975; %round(min(real) / binSizes) * binSizes;
    maxi = .975; %round(max(real) / binSizes) * binSizes;
    bins = mini:binSizes:maxi;
    n1 = hist(real(isSignf), bins);
    n2 = hist(real(~isSignf), bins);
    bar(bins, [n1',n2'], 'stacked')
    colormap([0 0 0;1 1 1])
    xlim([mini-.5*binSizes maxi+.5*binSizes])
    title(sprintf('%s (n=%d)', labels{m}, sum(~isnan(real))))
    xlabel(funLabels)
    ylabel('#Neurons')
    
%     figure('Position', [1250 680 560 420])
%     ind = ~isnan(real);
%     x = sort(real(ind), 'ascend');
%     y = (1:sum(ind)) ./ sum(ind);
%     ind = ~isnan(pseudo);
%     xNull = sort(pseudo(ind), 'ascend');
%     yNull = (1:sum(ind(:))) ./ sum(ind(:));
%     xNull = [min(x); xNull; max(x)];
%     yNull = [0 yNull 1];
%     hold on
%     plot(x, y, 'k', 'LineWidth', 2)
%     plot(xNull, yNull, ':', 'Color', [.5 .5 .5], 'LineWidth', 2)
%     exReal = real(exIndices);
%     exSign = isSignf(exIndices);
%     [xUni, ind] = unique(x);
%     yUni = y(ind);
%     heights = interp1(xUni, yUni, exReal);
%     for j = 1:length(exIndices)
%         if exSign(j)
%             plot(exReal(j), heights(j), 'o', 'MarkerEdgeColor', cols(j,:), ...
%                 'MarkerFaceColor', cols(j,:), 'LineWidth', 2)
%         else
%             plot(exReal(j), heights(j), 'o', 'MarkerEdgeColor', cols(j,:), ...
%                 'LineWidth', 2)
%         end
%     end
%     xlim(x([1 end]))
%     xlabel(funLabels)
%     ylabel('Proportion of neurons')
%     title(sprintf('%s (n=%d)', labels{m}, sum(~isnan(real))))
%     legend('data','null distribution', 'Location', 'NorthWest')
%     legend('boxoff')
end
figure
hold on
for j = 1:length(exIndices)
    plot(j, 1, 'o', 'MarkerEdgeColor', cols(j,:), ...
        'MarkerFaceColor', cols(j,:), 'MarkerSize', 30)
end

% Plot scatter plots
for s = 1:2
    m1 = (s-1)*2 + 1;
    m2 = (s-1)*2 + 2;
    real1 = modFuns(measures{m1}(:,1), measures{m1}(:,2));
    real2 = modFuns(measures{m2}(:,1), measures{m2}(:,2));
    mini = -.975; %round(min(real) / binSizes) * binSizes;
    maxi = .975; %round(max(real) / binSizes) * binSizes;
    figure
    scatter(real1, real2, 20, 'k', 'filled')
    alpha(.3)
    hold on
    for j = 1:length(exIndices)
        plot(real1(exIndices(j)), real2(exIndices(j)), 'o', ...
            'MarkerEdgeColor', cols(j,:), 'LineWidth', 2)
    end
    xlim([-1 1])
    ylim([-1 1])
    title(sprintf('n = %d', sum(~isnan(real1)&~isnan(real2))))
    xlabel(labels{m1})
    ylabel(labels{m2})
    axis square
end