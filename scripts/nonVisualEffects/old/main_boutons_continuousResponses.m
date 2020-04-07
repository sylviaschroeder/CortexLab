folderROIData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\InfoStructs';
folderResultsNeurons = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\nonVisualEffects\modelGratingResp\';
folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\boutons\nonVisualEffects\';

db_boutons_driftingGratings_blanks

%% Parameters
% fields = {'expGratings','expGrayScreen','expDark'};
% stimuli = {'gratings', 'grayScreen', 'dark'};
fields = {'expGratings','expDark'};
stimuli = {'gratings','dark'};

%% Plot traces of most and least correlated neurons
% num neurons for each condition
n = 10;
cols = [.7 0 0; 0 0 .7];
minSamples = 30;
runningThreshold = 5;
reference = 'running';

data = load(fullfile(folderResults,'running', ...
    'corrsDuringGratingsAndGrayScreen.mat'));
corrsRunning = data.corrs;
running = data.nonVisual;
% data = load(fullfile(folderResults,'pupil', ...
%     'corrsDuringGratingsAndGrayScreen.mat'));
% corrs = data.corrs;
% pupil = data.nonVisual;

expType = 2; % see 'stimuli'

if ~exist(fullfile(folderResults, ['plots_' reference], ...
        'continuousTraceCorrWithBehaviour'), 'dir')
    mkdir(fullfile(folderResults, ['plots_' reference], ...
        'continuousTraceCorrWithBehaviour'));
end

for k = 1:length(db)
    exp = db(k).(fields{expType});
    if isempty(exp)
        continue
    end
    folder = fullfile(folderROIData, db(k).subject, ...
        db(k).date, num2str(exp));
    fileStart = [db(k).date '_' num2str(exp) '_' ...
        db(k).subject];
    file = [fileStart '_2P_plane%03d_ROI.mat'];
    calciumTraces = cell(1, length(db(k).planes));
    time = cell(1, length(db(k).planes));
    rhos = [];
    cellIDs = [];
    for iPlane = 1:length(db(k).planes)
        % load meta
        data=load(fullfile(folder, sprintf(file,db(k).planes(iPlane))));
        meta = data.meta;
        time{iPlane} = ppbox.getFrameTimes(meta);
        calciumTraces{iPlane} = medfilt1(meta.F_final,11,[],1);
        
        if strcmp(reference, 'running')
            r = corrsRunning(k).plane(iPlane).(stimuli{expType}).rhos;
        else
            r = corrs(k).plane(iPlane).(stimuli{expType}).rhos;
        end
        rhos = [rhos; r];
        cellIDs = [cellIDs; [ones(length(r),1).*iPlane (1:length(r))']];
    end
    valid = find(~isnan(rhos));
    [~, sorted] = sort(rhos(valid));
    set{1} = valid(sorted(end-n+1:end));
    set{2} = valid(sorted(1:n));
    
    if strcmp(reference, 'running')
        beh = running{k,expType};
        ind = beh' > runningThreshold;
    else
        beh = pupil{k,expType};
        ind = beh' > nanmedian(pupil{k,expType});
    end
    sta = find(diff(ind)==1);
    sto = find(diff(ind)==-1);
    if sta(1)>sto(1)
        sta = [1, sta];
    end
    if sta(end)>sto(end)
        sto(end+1) = length(beh);
    end
    ind = (sto - sta) >= minSamples;
    starts = sta(ind);
    stops = sto(ind);
    
    dist = 5;
    ax = [];
    figure('Position', [1 41 1920 1083])
    
    subplot(10,1,1)
%     hold on
%     mini = min(pupil{k,expType});
%     maxi = max(pupil{k,expType});
%     fill(time{1}([starts;stops;stops;starts]), ...
%         [mini mini maxi maxi]', 'k', 'EdgeColor', 'none', 'FaceColor', ...
%         [1 1 1].*.9)
%     plot(time{1},pupil{k,expType},'k')
%     ylabel({'Pupil';'diam.'})
%     axis tight
%     a = gca;
%     a.Box = 'off';
%     a.XTickLabel = [];
%     a.XTick = [];
%     ax(end+1) = a;
    title(sprintf('Correlation with behaviour - %s (reference: %s, %s)', ...
        stimuli{expType}, reference, stimuli{expType}))
    
    subplot(10,1,2)
    hold on
    mini = min(running{k,expType});
    maxi = max(running{k,expType});
    fill(time{1}([starts;stops;stops;starts]), ...
        [mini mini maxi maxi]', 'k', 'EdgeColor', 'none', 'FaceColor', ...
        [1 1 1].*.9)
    plot(time{1},running{k,expType},'k')
    ylabel({'Running';'(cm/s)'})
    axis tight
    a = gca;
    a.Box = 'off';
    a.XTickLabel = [];
    a.XTick = [];
    ax(end+1) = a;
    
    subplot(10,1,3:10)
    hold on
    mini = -(2*n + .7) * dist;
    maxi = 2.5*dist;
    fill(time{1}([starts;stops;stops;starts]), ...
        [mini mini maxi maxi]', 'k', 'EdgeColor', 'none', 'FaceColor', ...
        [1 1 1].*.9)
    y = 0;
    for s = 1:2
        for j = 1:n
            tr = calciumTraces{cellIDs(set{s}(j),1)}(:,cellIDs(set{s}(j),2));
            tr = (tr-nanmean(tr))/nanstd(tr);
            plot(time{cellIDs(set{s}(j),1)},tr+y, 'Color', cols(s,:))
            y = y-dist;
        end
        y = y - dist;
    end
    axis tight
    ylim([mini-.2 maxi])
    xlabel('Time (s)')
    ylabel('\DeltaF/F (normalised)')
    a = gca;
    a.Box = 'off';
    a.YTick = [];
    ax(end+1) = a;
    
    linkaxes(ax, 'x')
    
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(fullfile(folderResults, ['plots_' reference], 'continuousTraceCorrWithBehaviour', ...
        sprintf('traces_%s_%s_%s.jpg', db(k).subject, db(k).date, stimuli{expType})), ...
        '-djpeg','-r0')
    close gcf
%     
%     disp(['Cell IDs of plotted traces for ' stimuli{expType}])
%     disp('Most correlated neurons')
%     disp(cellsMostCorr)
%     disp('Least correlated neurons')
%     disp(cellsLeastCorr)
end

%% Plot population results (against those of neurons)
behaviour = 'running';

data = load(fullfile(folderResults, behaviour, ...
    'corrsDuringGratingsAndGrayScreen.mat'));
corrs = data.corrs;
data = load(fullfile(folderResultsNeurons, behaviour, ...
    'corrsDuringGratingsAndGrayScreen.mat'));
corrsNeurons = data.corrs;

example = 5;
numCells = 10;

rhos = cell(length(corrs), 3);
nullRhos = cell(length(corrs), 3);
rhosNeurons = cell(length(corrsNeurons), 2);
nullRhosNeurons = cell(length(corrsNeurons), 2);

for k = 1:length(corrs)
    for iPlane = 1:length(corrs(k).plane)
        for st = 1:3
            if k == 2 && st == 3 % animal not running here
                rhos{k,st} = [rhos{k,st}; ...
                    NaN(size(corrs(k).plane(iPlane).(stimuli{st}).rhos))];
                nullRhos{k,st} = [nullRhos{k,st}; ...
                    NaN(size(corrs(k).plane(iPlane).(stimuli{st}).nullRhos))];
                continue
            end
            rhos{k,st} = [rhos{k,st}; ...
                corrs(k).plane(iPlane).(stimuli{st}).rhos];
            nullRhos{k,st} = [nullRhos{k,st}; ...
                corrs(k).plane(iPlane).(stimuli{st}).nullRhos];            
        end
    end 
end
for k = 1:length(corrsNeurons)
    for iPlane = 1:length(corrsNeurons(k).plane)
        for st = 1:2
            rhosNeurons{k,st} = [rhosNeurons{k,st}; ...
                corrsNeurons(k).plane(iPlane).(stimuli{st}).rhos];
            nullRhosNeurons{k,st} = [nullRhosNeurons{k,st}; ...
                corrsNeurons(k).plane(iPlane).(stimuli{st}).nullRhos];            
        end
    end 
end

% (1) Histograms: distribution of corr. coeffs. for gratings and blank screen
rhosGratings = cat(1,rhos{:,1});
rhosBlanks = cat(1,rhos{:,2});
rhosDark = cat(1,rhos{:,3});
nullGratings = cat(1,nullRhos{:,1});
nullBlanks = cat(1,nullRhos{:,2});
nullDark = cat(1,nullRhos{:,3});

ind = isnan(rhosGratings);
rhosGratings(ind) = [];
nullGratings(ind,:) = [];
ind = isnan(rhosBlanks);
rhosBlanks(ind) = [];
nullBlanks(ind,:) = [];
ind = isnan(rhosDark);
rhosDark(ind) = [];
nullDark(ind,:) = [];
prctlGratings = prctile(nullGratings, [2.5 97.5], 2);
prctlBlanks = prctile(nullBlanks, [2.5 97.5], 2);
prctlDark = prctile(nullDark, [2.5 97.5], 2);
isSgnfGratings = rhosGratings < prctlGratings(:,1) | ...
    rhosGratings > prctlGratings(:,2);
isSgnfBlanks = rhosBlanks < prctlBlanks(:,1) | rhosBlanks > prctlBlanks(:,2);
isSgnfDark = rhosDark < prctlDark(:,1) | rhosDark > prctlDark(:,2);

rhosGratingsNeurons = cat(1,rhosNeurons{:,1});
rhosBlanksNeurons = cat(1,rhosNeurons{:,2});
nullGratingsNeurons = cat(1,nullRhosNeurons{:,1});
nullBlanksNeurons = cat(1,nullRhosNeurons{:,2});

ind = isnan(rhosGratingsNeurons);
rhosGratingsNeurons(ind) = [];
nullGratingsNeurons(ind,:) = [];
ind = isnan(rhosBlanksNeurons);
rhosBlanksNeurons(ind) = [];
nullBlanksNeurons(ind,:) = [];
prctlGratingsNeurons = prctile(nullGratingsNeurons, [2.5 97.5], 2);
prctlBlanksNeurons = prctile(nullBlanksNeurons, [2.5 97.5], 2);
isSgnfGratingsNeurons = rhosGratingsNeurons < prctlGratingsNeurons(:,1) | ...
    rhosGratingsNeurons > prctlGratingsNeurons(:,2);
isSgnfBlanksNeurons = rhosBlanksNeurons < prctlBlanksNeurons(:,1) | ...
    rhosBlanksNeurons > prctlBlanksNeurons(:,2);

% gratings
figure
bins = -.5:.05:.8;
n1 = hist(rhosGratings(isSgnfGratings),bins)';
n2 = hist(rhosGratings(~isSgnfGratings),bins)';
bar(bins, [n1,n2] ./ sum([n1;n2]), 'stacked')
colormap([0 0 0; 1 1 1])
hold on
n1 = hist(rhosGratingsNeurons(isSgnfGratingsNeurons),bins)';
n2 = hist(rhosGratingsNeurons(~isSgnfGratingsNeurons),bins)';
stairs(bins-.025, n1 ./ sum([n1;n2]), 'Color', [0 .5 0], 'LineWidth', 2)
stairs(bins-.025, (n1+n2) ./ sum([n1;n2]), ':', 'Color', [0 .5 0], 'LineWidth', 2)
xlabel('Correlation coeff.')
ylabel('Proportion of boutons (neurons)')
title(sprintf('Corr. with %s during gratings (n = %d)', behaviour, length(rhosGratings)))
xlim([-.55 .85])
ax = gca;
ax.Box = 'off';
legend('p < 0.05', 'p >= 0.05', 'p < 0.05 (neurons)', 'p >= 0.05 (neurons)')
legend('boxoff')
% blanks
figure
n1 = hist(rhosBlanks(isSgnfBlanks),bins)';
n2 = hist(rhosBlanks(~isSgnfBlanks),bins)';
bar(bins, [n1,n2] ./ sum([n1;n2]), 'stacked')
colormap([0 0 0; 1 1 1])
hold on
n1 = hist(rhosGratingsNeurons(isSgnfGratingsNeurons),bins)';
n2 = hist(rhosGratingsNeurons(~isSgnfGratingsNeurons),bins)';
stairs(bins-.025, n1 ./ sum([n1;n2]), 'Color', [0 .5 0], 'LineWidth', 2)
stairs(bins-.025, (n1+n2) ./ sum([n1;n2]), ':', 'Color', [0 .5 0], 'LineWidth', 2)
xlabel('Correlation coeff.')
ylabel('Proportion of boutons (neurons)')
title(sprintf('Corr. with %s during gray screens (n = %d)', behaviour, length(rhosBlanks)))
xlim([-.55 .85])
ax = gca;
ax.Box = 'off';
legend('p < 0.05', 'p >= 0.05', 'p < 0.05 (neurons)', 'p >= 0.05 (neurons)')
legend('boxoff')
% dark
figure
n1 = hist(rhosDark(isSgnfDark),bins)';
n2 = hist(rhosDark(~isSgnfDark),bins)';
bar(bins, [n1,n2] ./ sum([n1;n2]), 'stacked')
colormap([0 0 0; 1 1 1])
xlabel('Correlation coeff.')
ylabel('Proportion of boutons')
title(sprintf('Corr. with %s during darkness (n = %d)', behaviour, length(rhosDark)))
xlim([-.55 .85])
ylim([0 .2])
ax = gca;
ax.Box = 'off';
legend('p < 0.05', 'p >= 0.05')
legend('boxoff')

% (2) Cumulative distribution plots
rhosGratings = cat(1,rhos{:,1});
rhosBlanks = cat(1,rhos{:,2});
nullGratings = cat(1,nullRhos{:,1});
nullBlanks = cat(1,nullRhos{:,2});
r1 = rhos{example,1}; % gratings, pupil
r2 = rhos{example,2}; % gray screen, pupil
r3 = rhos{example,2}; % gray screen, pupil

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
[~,largest] = sort(r3,'descend');
[~,smallest] = sort(r3,'ascend');

rhosGratingsNeurons = cat(1,rhosNeurons{:,1});
rhosBlanksNeurons = cat(1,rhosNeurons{:,2});
nullGratingsNeurons = cat(1,nullRhosNeurons{:,1});
nullBlanksNeurons = cat(1,nullRhosNeurons{:,2});

ind = isnan(rhosGratingsNeurons);
rhosGratingsNeurons(ind) = [];
nullGratingsNeurons(ind,:) = [];
ind = isnan(rhosBlanksNeurons);
rhosBlanksNeurons(ind) = [];
nullBlanksNeurons(ind,:) = [];

% gratings
x = sort(rhosGratings, 'ascend');
y = (1:length(x)) ./ length(x);
xNull = sort(nullGratings(:), 'ascend');
yNull = (1:length(xNull))' ./ (length(xNull));
xNull = [min(x); xNull; max(x)];
yNull = [0; yNull; 1];
h = [0 0 0 0];
figure
hold on
h(1) = plot(x, y, 'k', 'LineWidth', 2);
h(2) = plot(xNull, yNull, ':', 'Color', [.5 .5 .5], 'LineWidth', 2);
heights = interp1(x, y, r1(largest(1:numCells)));
plot(r1(largest(1:numCells)), heights, 'o', 'MarkerEdgeColor', 'none', ...
    'MarkerFaceColor', 'r')
heights = interp1(x, y, r1(smallest(1:numCells)));
plot(r1(smallest(1:numCells)), heights, 'o', 'MarkerEdgeColor', 'none', ...
    'MarkerFaceColor', 'b')

x = sort(rhosGratingsNeurons, 'ascend');
y = (1:length(x)) ./ length(x);
xNull = sort(nullGratingsNeurons(:), 'ascend');
yNull = (1:length(xNull))' ./ (length(xNull));
xNull = [min(x); xNull; max(x)];
yNull = [0; yNull; 1];
h(3) = plot(x, y, 'Color', [0 .5 0], 'LineWidth', 2);
h(4) = plot(xNull, yNull, ':', 'Color', [0 .5 0], 'LineWidth', 2);

xlim([-.55 .85])
xlabel('Correlation with pupil')
ylabel('Proportion of boutons (neurons)')
title(sprintf('Corr. with pupil during gratings (n = %d)', length(rhosGratings)))
legend(h, 'data', 'null', 'data (neurons)', 'null (neurons)', 'Location', 'NorthWest')
legend('boxoff')

% gray screens
x = sort(rhosBlanks, 'ascend');
y = (1:length(x)) ./ length(x);
xNull = sort(nullBlanks(:), 'ascend');
yNull = (1:length(xNull))' ./ (length(xNull));
xNull = [min(x); xNull; max(x)];
yNull = [0; yNull; 1];
h = [0 0 0 0];
figure
hold on
h(1) = plot(x, y, 'k', 'LineWidth', 2);
h(2) = plot(xNull, yNull, ':', 'Color', [.5 .5 .5], 'LineWidth', 2);
heights = interp1(x, y, r2(largest(1:numCells)));
plot(r2(largest(1:numCells)), heights, 'o', 'MarkerEdgeColor', 'none', ...
    'MarkerFaceColor', 'r')
heights = interp1(x, y, r2(smallest(1:numCells)));
plot(r2(smallest(1:numCells)), heights, 'o', 'MarkerEdgeColor', 'none', ...
    'MarkerFaceColor', 'b')

x = sort(rhosBlanksNeurons, 'ascend');
y = (1:length(x)) ./ length(x);
xNull = sort(nullBlanksNeurons(:), 'ascend');
yNull = (1:length(xNull))' ./ (length(xNull));
xNull = [min(x); xNull; max(x)];
yNull = [0; yNull; 1];
h(3) = plot(x, y, 'Color', [0 .5 0], 'LineWidth', 2);
h(4) = plot(xNull, yNull, ':', 'Color', [0 .5 0], 'LineWidth', 2);

xlim([-.55 .85])
xlabel('Correlation with pupil')
ylabel('Proportion of boutons (neurons)')
title(sprintf('Corr. with pupil during gray screens (n = %d)', length(rhosBlanks)))
legend(h, 'data', 'null', 'data (neurons)', 'null (neurons)', 'Location', 'NorthWest')
legend('boxoff')