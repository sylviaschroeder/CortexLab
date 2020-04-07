folderROIData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\InfoStructs';
% folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\nonVisualEffects\modelGratingResp\';
folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\boutons\nonVisualEffects\';

%% Parameters
% smoothing (low-pass filter) of non-visual signal before
smoothStd = 0.25; %in sec

fields = {'expGratings','expGrayScreen','expDark'};
stimuli = {'gratings', 'grayScreen', 'dark'};
signals = {'pupil', 'running'};
% fields = {'expGratings','expDark'};
% stimuli = {'gratings', 'dark'};
% signals = {'running'};

%% Compare correlation with nonvisual signal between drifting grating and gray monitor
% Load database
db_driftingGratings_blank
% db_boutons_driftingGratings_blanks

draws = 200;
corrs = struct([]);

data = load(fullfile(folderResults, 'kernelFit', 'results.mat'));
results = data.results;

for s = 1:length(signals)
    nonVisualSignal = signals{s};
    nonVisual = cell(length(db), length(fields));
    for k = 1:length(db)
        fprintf('Dataset %d: %s %s\n', k, db(k).subject, ...
            db(k).date);
        corrs(k).subject = db(k).subject;
        corrs(k).date = db(k).date;
        corrs(k).exp = cell(1,length(fields));
        for f = 1:length(fields)
            corrs(k).exp{f} = db(k).(fields{f});
        end
        
        cellIDs = cell(length(fields), length(db(k).planes));
        calciumTraces = cell(length(fields), length(db(k).planes));
        nv = cell(length(fields), length(db(k).planes));
        for st = 1:length(fields) % analyse different stimuli
            if isempty(db(k).(fields{st}))
                continue
            end
            folder = fullfile(folderROIData, db(k).subject, ...
                db(k).date, num2str(db(k).(fields{st})));
            fileStart = [db(k).date '_' num2str(db(k).(fields{st})) '_' ...
                db(k).subject];
            file = [fileStart '_2P_plane%03d_ROI.mat'];
            for iPlane = 1:length(db(k).planes)
                % load meta
                data=load(fullfile(folder, sprintf(file,db(k).planes(iPlane))));
                meta = data.meta;
                if iPlane == 1
                    nonVisData = [];
                    % load ball or pupil data
                    switch nonVisualSignal
                        case 'running'
                            ballData = nonVis.getRunningSpeed(meta);
                            if ~isempty(ballData)
                                stdSamples = round(smoothStd / median(diff(ballData.t)));
                                convWindow = normpdf(-3*stdSamples:3*stdSamples, 0, stdSamples);
                                nonVisData = conv(ballData.total, convWindow, 'same');
                                nonVisData = nonVisData / median(diff(ballData.t)) / 53;
                                nonVisTime = ballData.t;
                            end
                            label = 'running speed';
                        case 'pupil'
                            [pupilData, nonVisTime] = nonVis.loadPupilData(meta);
                            if ~isempty(pupilData)
                                nonVisTime(length(pupilData.x)+1:end) = [];
                                nonVisData = nonVis.getPupilDiam(pupilData);
                            end
                            label = 'pupil diameter';
                    end
                end
%                 meta.folderTL = strrep(meta.folderTL, 'ioo.ucl.ac.uk', ...
%                     'cortexlab.net');
                frameTimes = ppbox.getFrameTimes(meta);
                ind = isnan(nonVisData);
                indInterp = hist(nonVisTime(ind), frameTimes) > 0;
                warning off
                nv{st,iPlane} = interp1(nonVisTime, nonVisData, frameTimes, 'pchip')';
                nv{st,iPlane}(indInterp) = NaN;
                % only consider ROIs that are unique and not switch-on
                ind = ~all(isnan(meta.F_final),1)';
                if isfield(meta.ROI, 'isDuplicate')
                    ind = ind & meta.ROI.isDuplicate == 0;
                end
                if isfield(meta.ROI, 'isSwitchOn')
                    ind = ind & meta.ROI.isSwitchOn == 0;
                end
                cellIDs{st,iPlane} = ind; 
                calciumTraces{st,iPlane} = meta.F_final;
            end
            nonVisual{k,st} = nv{st,1};
        end
        set = find(strcmp(db(k).subject, {results.subject}) & ...
            strcmp(db(k).date, {results.date}));
        for iPlane = 1:length(db(k).planes)
            ind = all(cat(2,cellIDs{:,iPlane}), 2);
            for st = 1:length(fields)
                if isempty(calciumTraces{st,iPlane})
                    corrs(k).plane(iPlane).(stimuli{st}).rhos = NaN(length(ind),1);
                    corrs(k).plane(iPlane).(stimuli{st}).pValues = NaN(length(ind),1);
                    corrs(k).plane(iPlane).(stimuli{st}).nullRhos = NaN(length(ind),draws);
                    continue
                end
                j = ~isnan(nv{st,iPlane}) & ...
                    all(~isnan(calciumTraces{st,iPlane}(:,ind)),2);
                [r,p] = corr(nv{st,iPlane}(j), calciumTraces{st,iPlane}(j,ind));
                rhos = NaN(size(calciumTraces{st,iPlane},2), 1);
                rhos(ind) = r;
                pValues = NaN(size(calciumTraces{st,iPlane},2), 1);
                pValues(ind) = p;
                corrs(k).plane(iPlane).(stimuli{st}).rhos = rhos;
                corrs(k).plane(iPlane).(stimuli{st}).pValues = pValues;
                
                % generate null distribution by shifting non-visual signal in time
                shifts = randi(length(nv{st,iPlane}), draws, 1);
                nonVisualShifted = NaN(length(nv{st,iPlane}), draws);
                r = NaN(sum(ind), draws);
                for j = 1:draws
                    nonVisualShifted(:,j) = circshift(nv{st,iPlane}, shifts(j));
                    t = ~isnan(nonVisualShifted(:,j)) & ...
                        all(~isnan(calciumTraces{st,iPlane}(:,ind)),2);
                    r(:,j) = corr(nonVisualShifted(t,j), ...
                        calciumTraces{st,iPlane}(t,ind));
                end
                rhos = NaN(size(calciumTraces{st,iPlane},2), draws);
                rhos(ind,:) = r;
                corrs(k).plane(iPlane).(stimuli{st}).nullRhos = rhos;
            end
            
            l = length(ind);
            ind = find(ind);
            if ~isempty(set)
                both = intersect(results(set).plane(iPlane).cellIDs, ind);
                j1 = ismember(both, results(set).plane(iPlane).cellIDs);
                j2 = ismember(both, ind);
                responsive = NaN(l, 1);
                if ~isempty(ind)
                    responsive(ind(j2)) = ~cellfun(@isempty, ...
                        {results(set).plane(iPlane).kernelFit(j1).kernel})';
                end
                corrs(k).plane(iPlane).responsive = responsive;
                if isfield(results(set).plane, 'isGad')
                    isGad = NaN(l, 1);
                    isGad(ind(j2)) = results(set).plane(iPlane).isGad(j1);
                    corrs(k).plane(iPlane).isGad = isGad;
                end
            end
        end
    end
    clear data
    data.corrs = corrs;
    data.nonVisual = nonVisual;
    save(fullfile(folderResults,nonVisualSignal, ...
        'corrsDuringGratingsAndGrayScreen.mat'), '-struct', 'data')
end

%% Plot population results
data = load(fullfile(folderResults,'running', ...
    'corrsDuringGratingsAndGrayScreen.mat'));
corrsRunning = data.corrs;
running = data.nonVisual;
% data = load(fullfile(folderResults,'pupil', ...
%     'corrsDuringGratingsAndGrayScreen.mat'));
% corrsPupil = data.corrs;
% pupil = data.nonVisual;

% example = 5; % neurons: 7, boutons: 2

% rhosPupil = cell(length(corrsPupil), length(stimuli));
% nullRhosPupil = cell(length(corrsPupil), length(stimuli));
rhosRunning = cell(length(corrsRunning), length(stimuli));
nullRhosRunning = cell(length(corrsRunning), length(stimuli));
% resp = cell(length(corrsPupil),1);
% isGad = cell(length(corrsPupil),1);
for k = 1:length(corrsRunning)
    for iPlane = 1:length(corrsRunning(k).plane)
        for st = 1:length(stimuli)
%             rhosPupil{k,st} = [rhosPupil{k,st}; ...
%                 corrsPupil(k).plane(iPlane).(stimuli{st}).rhos];
%             nullRhosPupil{k,st} = [nullRhosPupil{k,st}; ...
%                 corrsPupil(k).plane(iPlane).(stimuli{st}).nullRhos];
            rhosRunning{k,st} = [rhosRunning{k,st}; ...
                corrsRunning(k).plane(iPlane).(stimuli{st}).rhos];
            nullRhosRunning{k,st} = [nullRhosRunning{k,st}; ...
                corrsRunning(k).plane(iPlane).(stimuli{st}).nullRhos];
        end
%         resp{k} = [resp{k}; corrsPupil(k).plane(iPlane).responsive];
%         if isfield(corrsPupil(k).plane, 'isGad')
%             isGad{k} = [isGad{k}; corrsPupil(k).plane(iPlane).isGad];
%         end
    end 
end

% (1) Scatterplot: corr. during gratings vs. corr. during blank screen
figure
hold on
% pupil
rhosGratings = cat(1,rhosPupil{:,1});
rhosBlanks = cat(1,rhosPupil{:,2});
% running
% rhosGratings = cat(1,rhosRunning{:,1});
% rhosBlanks = cat(1,rhosRunning{:,2});
ind = ~any(isnan([rhosGratings, rhosBlanks]), 2);
rho = corr(cat(1,rhosGratings(ind)), cat(1,rhosBlanks(ind)));
n = sum(ind);
isResp = cat(1, resp{:});
cols = repmat([.5 1]',1,3);
conds = [0 1];
for c = 1:2
    plot(rhosGratings(isResp==conds(c)), rhosBlanks(isResp==conds(c)), 'k.', 'MarkerSize', 3)
end
xlabel('during drifting gratings')
ylabel('during gray screen')
% title(sprintf('Correlation with pupil diameter (n=%d, rho=%.3f)',n,rho))
title(sprintf('Correlation with running speed (n=%d, rho=%.3f)',n,rho))
axis([-.55 .85 -.55 .85]) 
axis square
% mark those neurons whose traces were plotted in figures from
% corrSingleCellTracesWithBehaviour.m (assuming both were sorted according
% to the same bahvioural data, i.e. running or pupil)
r1 = rhosPupil{example,1}; % gratings, pupil
r2 = rhosPupil{example,2}; % gray screen, pupil
% r1 = rhosRunning{example,1}; % gratings, pupil
% r2 = rhosRunning{example,2}; % gray screen, pupil

r3 = rhosPupil{example,2};
% r3 = rhosRunning{example,2}; % gray screen, running
% ind = ~isnan(r1);
ind = ~isnan(r3);
r1 = r1(ind);
r2 = r2(ind);
r3 = r3(ind);
% [~,largest] = sort(r1,'descend');
% [~,smallest] = sort(r1,'ascend');
[~,largest] = sort(r3,'descend');
[~,smallest] = sort(r3,'ascend');
numCells = 10;
hold on
plot(r1(largest(1:numCells)),r2(largest(1:numCells)),'r.', 'MarkerSize', 10)
plot(r1(smallest(1:numCells)),r2(smallest(1:numCells)),'b.', 'MarkerSize', 10)

% (2) Histograms: distribution of corr. coeffs. for gratings and blank screen
% pupil
% rhosGratings = cat(1,rhosPupil{:,1});
% rhosBlanks = cat(1,rhosPupil{:,2});
% nullGratings = cat(1,nullRhosPupil{:,1});
% nullBlanks = cat(1,nullRhosPupil{:,2});
% running
rhosGratings = cat(1,rhosRunning{:,1});
rhosBlanks = cat(1,rhosRunning{:,2});
nullGratings = cat(1,nullRhosRunning{:,1});
nullBlanks = cat(1,nullRhosRunning{:,2});
ind = isnan(rhosGratings);
rhosGratings(ind) = [];
nullGratings(ind,:) = [];
ind = isnan(rhosBlanks);
rhosBlanks(ind) = [];
nullBlanks(ind,:) = [];
prctlGratings = prctile(nullGratings, [2.5 97.5], 2);
prctlBlanks = prctile(nullBlanks, [2.5 97.5], 2);
isSgnfGratings = rhosGratings < prctlGratings(:,1) | ...
    rhosGratings > prctlGratings(:,2);
isSgnfBlanks = rhosBlanks < prctlBlanks(:,1) | rhosBlanks > prctlBlanks(:,2);
% gratings
figure
bins = -.5:.05:.8;
n1 = hist(rhosGratings(isSgnfGratings),bins)';
n2 = hist(rhosGratings(~isSgnfGratings),bins)';
bar(bins, [n1,n2], 'stacked')
colormap([0 0 0; 1 1 1])
xlabel('Correlation coeff.')
ylabel('#Neurons')
% title(sprintf('Corr. with pupil during gratings (n = %d)', length(rhosGratings)))
title(sprintf('Corr. with running during gratings (n = %d)', length(rhosGratings)))
xlim([-.55 .85])
ax = gca;
ax.Box = 'off';
legend('p < 0.05', 'p >= 0.05')
% blanks
figure
n1 = hist(rhosBlanks(isSgnfBlanks),bins)';
n2 = hist(rhosBlanks(~isSgnfBlanks),bins)';
bar(bins, [n1,n2], 'stacked')
colormap([0 0 0; 1 1 1])
xlabel('Correlation coeff.')
ylabel('#Neurons')
% title(sprintf('Corr. with pupil during gray screens (n = %d)', length(rhosBlanks)))
title(sprintf('Corr. with running during darkness (n = %d)', length(rhosBlanks)))
xlim([-.55 .85])
ax = gca;
ax.Box = 'off';
legend('p < 0.05', 'p >= 0.05')

% (3) Cumulative distribution plots: for pupil
% pupil
rhosGratings = cat(1,rhosPupil{:,1});
rhosBlanks = cat(1,rhosPupil{:,2});
nullGratings = cat(1,nullRhosPupil{:,1});
nullBlanks = cat(1,nullRhosPupil{:,2});
r1 = rhosPupil{example,1}; % gratings, pupil
r2 = rhosPupil{example,2}; % gray screen, pupil
r3 = rhosPupil{example,2}; % gray screen, pupil
% r3 = rhosRunning{example,2}; % gray screen, running
% running
% rhosGratings = cat(1,rhosRunning{:,1});
% rhosBlanks = cat(1,rhosRunning{:,2});
% nullGratings = cat(1,nullRhosRunning{:,1});
% nullBlanks = cat(1,nullRhosRunning{:,2});
% r1 = rhosRunning{7,1}; % gratings, running
% r2 = rhosRunning{7,2}; % gray screen, running
% r3 = rhosRunning{7,2}; % gray screen, running

ind = isnan(rhosGratings);
rhosGratings(ind) = [];
nullGratings(ind,:) = [];
ind = isnan(rhosBlanks);
rhosBlanks(ind) = [];
nullBlanks(ind,:) = [];
ind = ~isnan(r3);
r1 = r1(ind);
r2 = r2(ind);
r3 = r3(ind);
% gratings
x = sort(rhosGratings, 'ascend');
y = (1:length(x)) ./ length(x);
xNull = sort(nullGratings(:), 'ascend');
yNull = (1:length(xNull))' ./ (length(xNull));
xNull = [min(x); xNull; max(x)];
yNull = [0; yNull; 1];
figure
hold on
plot(x, y, 'k', 'LineWidth', 2)
plot(xNull, yNull, ':', 'Color', [.5 .5 .5], 'LineWidth', 2)
heights = interp1(x, y, r1(largest(1:numCells)));
plot(r1(largest(1:numCells)), heights, 'o', 'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', 'r')
heights = interp1(x, y, r1(smallest(1:numCells)));
plot(r1(smallest(1:numCells)), heights, 'o', 'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', 'b')
xlim([-.55 .85])
xlabel('Correlation with pupil')
ylabel('Proportion of neurons')
% title(sprintf('Corr. with pupil during gratings (n = %d)', length(rhosGratings)))
title(sprintf('Corr. with running during gratings (n = %d)', length(rhosGratings)))
legend('data', 'null distribution', 'Location', 'NorthWest')
% gray screens
x = sort(rhosBlanks, 'ascend');
y = (1:length(x)) ./ length(x);
xNull = sort(nullBlanks(:), 'ascend');
yNull = (1:length(xNull))' ./ (length(xNull));
xNull = [min(x); xNull; max(x)];
yNull = [0; yNull; 1];
figure
hold on
plot(x, y, 'k', 'LineWidth', 2)
plot(xNull, yNull, ':', 'Color', [.5 .5 .5], 'LineWidth', 2)
heights = interp1(x, y, r2(largest(1:numCells)));
plot(r2(largest(1:numCells)), heights, 'o', 'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', 'r')
heights = interp1(x, y, r2(smallest(1:numCells)));
plot(r2(smallest(1:numCells)), heights, 'o', 'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', 'b')
xlim([-.55 .85])
xlabel('Correlation with pupil')
ylabel('Proportion of neurons')
% title(sprintf('Corr. with pupil during gray screens (n = %d)', length(rhosBlanks)))
title(sprintf('Corr. with running during gray screens (n = %d)', length(rhosBlanks)))
legend('data', 'null distribution', 'Location', 'NorthWest')

% (4) Cumulative distribution plots: for pupil, for inh. and exc. neurons
% pupil
% rhosGratings = cat(1,rhosPupil{:,1});
% rhosBlanks = cat(1,rhosPupil{:,2});
% nullGratings = cat(1,nullRhosPupil{:,1});
% nullBlanks = cat(1,nullRhosPupil{:,2});
% running
rhosGratings = cat(1,rhosRunning{:,1});
rhosBlanks = cat(1,rhosRunning{:,2});
nullGratings = cat(1,nullRhosRunning{:,1});
nullBlanks = cat(1,nullRhosRunning{:,2});
isGadGratings = cat(1,isGad{:});
isGadBlanks = cat(1,isGad{:});
ind = isnan(rhosGratings);
rhosGratings(ind) = [];
nullGratings(ind,:) = [];
isGadGratings(ind,:) = [];
ind = isnan(rhosBlanks);
rhosBlanks(ind) = [];
nullBlanks(ind,:) = [];
isGadBlanks(ind) = [];
% gratings
xInh = sort(rhosGratings(isGadGratings == 1), 'ascend');
yInh = (1:length(xInh)) ./ length(xInh);
xNullInh = sort(reshape(nullGratings(isGadGratings == 1,:),[],1), 'ascend');
yNullInh = (1:length(xNullInh))' ./ (length(xNullInh));
xNullInh = [min(xInh); xNullInh; max(xInh)];
yNullInh = [0; yNullInh; 1];
xExc = sort(rhosGratings(isGadGratings == -1), 'ascend');
yExc = (1:length(xExc)) ./ length(xExc);
xNullExc = sort(reshape(nullGratings(isGadGratings == -1,:),[],1), 'ascend');
yNullExc = (1:length(xNullExc))' ./ (length(xNullExc));
xNullExc = [min(xExc); xNullExc; max(xExc)];
yNullExc = [0; yNullExc; 1];
cols = lines(2);
h = zeros(1,4);
figure
hold on
plot([0 0], [0 1], 'k:')
h(1) = plot(xInh, yInh, 'Color', cols(2,:), 'LineWidth', 2);
h(2) = plot(xNullInh, yNullInh, ':', 'Color', cols(2,:), 'LineWidth', 2);
h(3) = plot(xExc, yExc, 'Color', cols(1,:), 'LineWidth', 2);
h(4) = plot(xNullExc, yNullExc, ':', 'Color', cols(1,:), 'LineWidth', 2);
xlim([-.55 .85])
xlabel('Correlation with pupil')
ylabel('Proportion of neurons')
[~,p] = kstest2(xInh, xExc);
% title(sprintf('Corr. with pupil during gratings (ks-test: p = %.3f)', p))
title(sprintf('Corr. with running during gratings (ks-test: p = %.3f)', p))
legend(h, sprintf('inhibitory (n=%d)',length(xInh)), 'null of inh.', ...
    sprintf('excitatory (n=%d)',length(xExc)), 'null of exc.', ...
    'Location', 'NorthWest')
% gray screens
xInh = sort(rhosBlanks(isGadBlanks == 1), 'ascend');
yInh = (1:length(xInh)) ./ length(xInh);
xNullInh = sort(reshape(nullBlanks(isGadBlanks == 1,:),[],1), 'ascend');
yNullInh = (1:length(xNullInh))' ./ (length(xNullInh));
xNullInh = [min(xInh); xNullInh; max(xInh)];
yNullInh = [0; yNullInh; 1];
xExc = sort(rhosBlanks(isGadBlanks == -1), 'ascend');
yExc = (1:length(xExc)) ./ length(xExc);
xNullExc = sort(reshape(nullBlanks(isGadBlanks == -1,:),[],1), 'ascend');
yNullExc = (1:length(xNullExc))' ./ (length(xNullExc));
xNullExc = [min(xExc); xNullExc; max(xExc)];
yNullExc = [0; yNullExc; 1];
cols = lines(2);
h = zeros(1,4);
figure
hold on
plot([0 0], [0 1], 'k:')
h(1) = plot(xInh, yInh, 'Color', cols(2,:), 'LineWidth', 2);
h(2) = plot(xNullInh, yNullInh, ':', 'Color', cols(2,:), 'LineWidth', 2);
h(3) = plot(xExc, yExc, 'Color', cols(1,:), 'LineWidth', 2);
h(4) = plot(xNullExc, yNullExc, ':', 'Color', cols(1,:), 'LineWidth', 2);
xlim([-.55 .85])
xlabel('Correlation with pupil')
ylabel('Proportion of neurons')
[~,p] = kstest2(xInh, xExc);
% title(sprintf('Corr. with pupil during gray screens (ks-test: p = %.3f)', p))
title(sprintf('Corr. with running during gray screens (ks-test: p = %.3f)', p))
legend(h, sprintf('inhibitory (n=%d)',length(xInh)), 'null of inh.', ...
    sprintf('excitatory (n=%d)',length(xExc)), 'null of exc.', ...
    'Location', 'NorthWest')

% Compare running with pupil
corrBeh = NaN(size(pupil,1),2);
for k=1:size(pupil,1)
    for st = 1:2
        if isempty(running{k,st}) || isempty(pupil{k,st})
            continue
        end
        ind = ~isnan(pupil{k,st});
        corrBeh(k,st) = corr(running{k,st}(ind),pupil{k,st}(ind));
    end
end
% plot running vs. pupil for all datasets
p = [];
r = [];
for k = 1:length(pupil)
    tmp = cat(1, pupil{k,:});
    p = [p; (tmp - min(tmp)) ./ (max(tmp) - min(tmp))];
    tmp = cat(1, running{k,:});
    r = [r; (tmp - min(tmp)) ./ (max(tmp) - min(tmp))];
end
ind = ~isnan(p) & ~isnan(r);
p = p(ind);
r = r(ind);
figure
plot(p,r,'k.','MarkerSize',1)
ind = ~isnan(p);
rho = corr(r(ind),p(ind));
axis([0 1 0 1])
axis square
xlabel('Pupil diameter')
ylabel('Running speed')
title(sprintf('Corr. coeff.: %.3f (average across datasets: %.3f)', rho, ...
    nanmean(corrBeh(:))))
% make density plot
figure
[N,C] = hist3([r,p],[30 30]);
contourf(C{2},C{1},log10(N+1),20,'LineColor','none')
colormap gray
l = [0 10 100 1000];
c = colorbar('Ticks',log10(l+1),'TickLabels',l);
c.Label.String = '# Samples';
axis square
xlabel('Running speed')
ylabel('Pupil diameter')
title(sprintf('Corr. coeff.: %.3f (average across datasets: %.3f)', rho, ...
    nanmean(corrBeh(:))))
% make conditional probability plot
pBins = prctile(p,0:10:100);
pBins(end) = pBins(end) + 1;
rBins = unique(prctile(r,0:10:100));
prob = NaN(length(rBins)-1, length(pBins)-1);
for b = 1:length(pBins)-1
    ind = p>=pBins(b) & p<pBins(b+1);
    n = histcounts(r(ind), rBins);
    prob(:,b) = n ./ sum(ind);
end
figure
imagesc([10 100], [20 100], prob)
h = colorbar;
h.Label.String = 'P(running|pupil)';
colormap gray
ax = gca;
ax.YDir = 'normal';
axis square
xlabel('Percentile of pupil diameter')
ylabel('Percentile of running speed')

% Compare corr. of neural traces with running vs. with pupil 
% (2) during gray screens
r = cat(1,rhosRunning{:,2});
p = cat(1,rhosPupil{:,2});
figure
plot(r,p,'k.', 'MarkerSize', 3)
ind = ~any(isnan([r,p]),2);
rho = corr(r(ind),p(ind));
xlabel('Corr. with running')
ylabel('Corr. with pupil diameter')
title(sprintf('Correlation during gray screen (n = %d, rho = %.3f)', sum(ind), rho))
axis square
hold on
rb1 = rhosPupil{7,2};
rb2 = rhosRunning{7,2};
ind = isnan(rb2);
rb1(ind) = [];
rb2(ind) = [];

[~,maxis]=sort(rb2,'descend');
[~,minis]=sort(rb2,'ascend');

plot(rb2(maxis(1:numCells)),rb1(maxis(1:numCells)), '.r', 'MarkerSize', 10)
plot(rb2(minis(1:numCells)),rb1(minis(1:numCells)), '.b', 'MarkerSize', 10)
axis([-.55 .85 -.55 .85])
ax = gca;
ax.Box = 'off';

% (1) during gratings
r = cat(1,rhosRunning{:,1}); % during gratings
p = cat(1,rhosPupil{:,1});
figure
plot(r,p,'k.', 'MarkerSize', 3)
ind = ~any(isnan([r,p]),2);
rho = corr(r(ind),p(ind));
xlabel('Corr. with running')
ylabel('Corr. with pupil diameter')
title(sprintf('Correlation during gratings (n = %d, rho = %.3f)', sum(ind), rho))
axis square
hold on
rg1 = rhosPupil{7,1};
rg2 = rhosRunning{7,1};
rb2 = rhosRunning{7,2};
% ind = isnan(rg1);
ind = isnan(rb2);
rg1(ind) = [];
rg2(ind) = [];
% [~,maxis]=sort(rg1,'descend');
% [~,minis]=sort(rg1,'ascend');
plot(rg2(maxis(1:numCells)),rg1(maxis(1:numCells)), '.r', 'MarkerSize', 10)
plot(rg2(minis(1:numCells)),rg1(minis(1:numCells)), '.b', 'MarkerSize', 10)
axis([-.55 .85 -.55 .85])
ax = gca;
ax.Box = 'off';

%% OLD: Compare correlation with behaviour across different visual responsiveness
% data = load(fullfile(folderResults, nonVisualSignal, ...
%     'linFit_nonvis_baseline_fullAdditiveModelOnly.mat'));
% results = data.results;
% data = load(fullfile(folderResults, nonVisualSignal, ...
%     'linFit_nonvis_baseline.mat'));
% resultsExpVar = data.results;
% rhoActivity = [];
% pActivity = [];
% rhoBaseline = [];
% pBaseline = [];
% explVarsFullAdditive = [];
% explVarsStimulusOnly = [];
% for iExp = 1:length(results)
%     fprintf('Dataset %d: %s %s exp.: %d\n', iExp, results(iExp).subject, ...
%         results(iExp).date, results(iExp).exp);
%     folder = [folderROIData filesep results(iExp).subject filesep ...
%         results(iExp).date filesep num2str(results(iExp).exp)];
%     fileStart = [results(iExp).date '_' num2str(results(iExp).exp) '_' ...
%         results(iExp).subject];
%     file = [fileStart '_2P_plane%03d_ROI.mat'];
%     for iPlane = 1:length(results(iExp).plane)
%         data=load(fullfile(folder, sprintf(file,results(iExp).planes(iPlane))));
%         meta = data.meta;
%         if iPlane == 1
%             nonVisData = [];
%             % load ball or pupil data
%             switch nonVisualSignal
%                 case 'running'
%                     ballData = nonVis.getRunningSpeed(meta);
%                     if ~isempty(ballData)
%                         filtWindow = ceil(smoothing / median(diff(ballData.t)));
%                         if mod(filtWindow,2) == 0
%                             filtWindow = filtWindow-1;
%                         end
%                         nonVisData = sgolayfilt(ballData.total, filtPoly, filtWindow);
%                         nonVisTime = ballData.t;
%                     end
%                     label = 'running speed';
%                 case 'pupil'
%                     [pupilData, nonVisTime] = nonVis.loadPupilData(meta);
%                     if ~isempty(pupilData)
%                         nonVisTime(length(pupilData.x)+1:end) = [];
%                         nonVisData = nonVis.getPupilDiam(pupilData);
%                     end
%                     label = 'pupil diameter';
%             end
%         end
%         frameTimes = ppbox.getFrameTimes(meta);
%         nv = interp1(nonVisTime, nonVisData, frameTimes, 'pchip')';
%         ids = results(iExp).plane(iPlane).cellIDs;
%         traces = bsxfun(@rdivide, meta.Fcorr(:,ids) - meta.F0(:,ids), ...
%             max(1, mean(meta.F0(:,ids),1)));
%         for iCell = 1:length(ids)
%             [r,p] = corr(nv,traces(:,iCell));
%             rhoActivity(end+1) = r;
%             pActivity(end+1) = p;
%             if ~isempty(results(iExp).plane(iPlane).modelPars(iCell).baseline)
%                 [r,p] = corr(nv, results(iExp).plane(iPlane) ...
%                     .modelPars(iCell).baseline);
%                 explVarsFullAdditive(end+1) = resultsExpVar(iExp).plane(iPlane) ...
%                     .crossVal(iCell).explVars(4);
%                 explVarsStimulusOnly(end+1) = resultsExpVar(iExp).plane(iPlane) ...
%                     .crossVal(iCell).explVars(12);
%             else
%                 r = NaN;
%                 p = NaN;
%                 explVarsFullAdditive(end+1) = NaN;
%                 explVarsStimulusOnly(end+1) = NaN;
%             end
%             rhoBaseline(end+1) = r;
%             pBaseline(end+1) = p;
%         end
%     end
% end
% 
% figure('Position',[680 680 1130 420])
% subplot(1,2,1)
% hold on
% subplot(1,2,2)
% hold on
% mini = floor(min(rhoActivity)*10)/10;
% maxi = ceil(max(rhoActivity)*10)/10;
% bins = mini+.05:.1:maxi;
% edges = mini : .1 : maxi+.05;
% cols = lines(3);
% groups = NaN(size(explVarsStimulusOnly));
% groups(isnan(explVarsStimulusOnly)) = 1;
% groups(explVarsStimulusOnly<=0) = 2;
% groups(explVarsStimulusOnly>0) = 3;
% names = {'No visual responses','No reliable visual response',...
%     'Reliable visual response'};
% for k = 1:3
%     n = hist(rhoActivity(groups==k), bins);
%     n = n ./ sum(groups==k);
%     subplot(1,2,1)
%     plot(bins, n, 'Color', cols(k,:), 'LineWidth', 2)
%     subplot(1,2,2)
%     plot(bins, cumsum(n), 'Color', cols(k,:), 'LineWidth', 2)
% %     stairs(edges, n([1:end end]), 'Color', cols(k,:), 'LineWidth', 2)
% %     subplot(3,1,k)
% %     hist(rhoActivity(groups==k), bins)
% %     h = findobj(gca,'Type','patch');
% %     h.FaceColor = cols(k,:);
% %     xlim([mini maxi])
% %     title(names{k})
% %     ylabel('# Neurons')
% end
% subplot(1,2,1)
% ylimit = get(gca,'YLim');
% plot([0 0],ylimit,'k:')
% xlim(edges([1 end]))
% ylim(ylimit)
% xlabel('Corr. coeff. with nonvisual signal')
% ylabel('Proportions of neurons')
% title(sprintf('Correlation with %s',label))
% subplot(1,2,2)
% ylimit = get(gca,'YLim');
% plot([0 0],ylimit,'k:')
% legend(names,'Location','southeast')
% xlim(edges([1 end]))
% ylim(ylimit)
% xlabel('Corr. coeff. with nonvisual signal')
% ylabel('Cum. sum')
% 
% figure
% scatter(rhoActivity(explVarsFullAdditive>0), rhoBaseline(explVarsFullAdditive>0), 'r')
% hold on
% scatter(rhoActivity(explVarsFullAdditive<=0), rhoBaseline(explVarsFullAdditive<=0), 'k')
% legend('Good model fits','Bad model fits')
% mini = floor(min([rhoActivity, rhoBaseline])*10)/10;
% maxi = ceil(max([rhoActivity, rhoBaseline])*10)/10;
% axis([mini maxi mini maxi])
% title('Correlation with X and nonvisual signal')
% xlabel('X = Original activity')
% ylabel('X = Baseline')
% axis square
% 
% figure
% scatter(explVarsStimulusOnly, rhoActivity, 'k')
% xlabel('Explained var. by stimulus only')
% ylabel('Corr. of activity with nonvisual signal')
