folderROIData = 'C:\DATA\InfoStructs';
folderResults = 'C:\RESULTS\nonVisualEffects\modelGratingResp\';

%% Parameters
minR2 = 0.3;
nonVisualSignal = 'pupil';

doPlot = 0;
doPause = 0;

data = load(fullfile(folderResults, nonVisualSignal, ...
    'slopeFits_lowVsHighNonVisualSignal_onlySignifParamsFit.mat'));
fits = data.fits;
data = load(fullfile(folderResults, nonVisualSignal, ...
    'tuning_prefDirSigmaFixedAcrossConditions_simultaneous.mat'));
tuning = data.tuning;

%% Check features of contrast-suppressed cells: exc./inh., tuning width
gadGroups = [-1 0 1];
isGad = cell(1,length(tuning));
isSuppressed = cell(1,length(tuning));
isTuned = cell(1,length(tuning));
tuningWidths = cell(1,length(tuning));
% isResponsive = [];
for iExp = 1:length(tuning)
    fprintf('Dataset %d: %s %s exp.: %d\n', iExp, tuning(iExp).subject, ...
        tuning(iExp).date, tuning(iExp).exp);
    for iPlane = 1:length(tuning(iExp).plane)
        isGad{iExp} = [isGad{iExp}; tuning(iExp).plane(iPlane).isGad];
        isSuppressed{iExp} = [isSuppressed{iExp}; tuning(iExp).plane(iPlane).isSuppressed];
        tuningPars = {tuning(iExp).plane(iPlane).cond(1).cell.parameters};
        widths = NaN(length(tuningPars),1);
        ind = ~cellfun(@isempty, tuningPars);
        pars = cell2mat(tuningPars(ind));
        widths(ind) = cat(1, pars(5,:));
        tuningWidths{iExp} = [tuningWidths{iExp}; widths];
        pValStim = {fits(iExp).plane(iPlane).cell.ANOVApStim};
        ind = ~cellfun(@isempty, pValStim);
        tuned = NaN(size(widths));
        tuned(ind) = [pValStim{ind}];
        isTuned{iExp} = [isTuned{iExp}; tuned];
    end
end

%% Plot results
% (1) Histgram of tuning width for suppressed and non-suppressed cells
isSupp = cat(1, isSuppressed{:});
tuned = cat(1, isTuned{:});
widths = cat(1, tuningWidths{:});
c = [-1 1];
str = {'Excited-by-contrast','Suppressed-by-contrast'};
for k = 1:2
    figure
    hold on
    bins = 30:5:180;
    ind = isSupp==c(k);
    n = hist(widths(ind & tuned<0.05), bins) ./ sum(ind) .* 100;
    bar(bins, n, 0.9, 'k')
    bar(200, nansum(ind & tuned>=0.05) ./ sum(ind) .* 100, 5, 'w')
    xlabel('Tuning width')
    ylabel('% Neurons')
    title(sprintf('%s (n = %d)', str{k}, sum(ind)))
end

% (2) Percentage of exc. and inh. cells with classs of excited- and
% suppressed-by-contrast cells
supprExc = NaN(length(tuning),1);
supprInh = NaN(length(tuning),1);
enhExc = NaN(length(tuning),1);
enhInh = NaN(length(tuning),1);
exc = NaN(length(tuning),1);
inh = NaN(length(tuning),1);
for k = 1:length(tuning)
    totalSuppr = nansum(isSuppressed{k}==1 & isGad{k}~=0);
    supprExc(k) = nansum(isSuppressed{k}==1 & isGad{k}==-1) / totalSuppr .* 100;
    supprInh(k) = nansum(isSuppressed{k}==1 & isGad{k}==1) / totalSuppr .* 100;
    totalEnh = nansum(isSuppressed{k}==-1 & isGad{k}~=0);
    enhExc(k) = nansum(isSuppressed{k}==-1 & isGad{k}==-1) / totalEnh .* 100;
    enhInh(k) = nansum(isSuppressed{k}==-1 & isGad{k}==1) / totalEnh .* 100;
    total = sum(~isnan(isSuppressed{k}) & isGad{k}~=0);
    exc(k) = nansum(~isnan(isSuppressed{k}) & isGad{k}==-1) / total .* 100;
    inh(k) = nansum(~isnan(isSuppressed{k}) & isGad{k}==1) / total .* 100;
end
excMed = median(exc);
inhMed = median(inh);
figure
hold on
plot([0 6],[1 1].*excMed,'b')
plot([0 6],[1 1].*inhMed,'r')
boxplot([enhExc, enhInh, NaN(size(enhExc)), ...
    supprExc, supprInh])
% boxplot([enhExc./excMed, enhInh./inhMed, NaN(size(enhExc)), ...
%     supprExc./excMed, supprInh./inhMed])
set(gca,'XTick',[1 2 4 5],'XTickLabel',{'Exc.', 'Inh.'},'box','off');
title('Enhanced- and suppressed-by-contrast cells')
ylabel('% Neurons within class')
ylim([0 100])


bars = zeros(2,3);
for g=1:3
    bars(1,g) = sum(isSuppressed==-1 & isGad==gadGroups(g));
    bars(2,g) = sum(isSuppressed==1 & isGad==gadGroups(g));
end
figure
bar(bsxfun(@rdivide, bars, sum(bars,1)) .* 100)
set(gca, 'XTick', [1 2], 'XTickLabel', {'not suppressed', 'suppressed'})
ylabel('Neurons (%)')
title('Cell type of contrast-suppressed cells')