function [results, figureHandles] = ...
    getCorrToNonVisData(calciumTraces, calciumTime, nonVisData, ...
    nonVisTime, nonVisName, neuronIDs, rhoThreshold, doPlot, stimOnsets, ...
    isGad)

signifP = 0.05;
colors = [.87 .49 0; 0 0 0; 0 .5 0];

if nargin < 6 || isempty(neuronIDs)
    neuronIDs = (1:size(calciumTraces,2))';
end
if nargin < 7 || isempty(rhoThreshold)
    rhoThreshold = 0;
end
if nargin < 8 || isempty(doPlot)
    doPlot = true;
end
if nargin < 9
    stimOnsets = [];
end
if nargin < 10
    isGad = [];
end
figureHandles = [];

zCal = bsxfun(@rdivide, bsxfun(@minus, calciumTraces, nanmean(calciumTraces,1)), ...
    nanstd(calciumTraces));
zNonVis = bsxfun(@rdivide, bsxfun(@minus, nonVisData, nanmean(nonVisData)), ...
    nanstd(nonVisData));

ind = nonVisTime < calciumTime(1) | nonVisTime > calciumTime(end);
nonVisTime(ind) = [];
zNonVis(ind) = [];

% interpolate non-visual data so that samples match sample times of
% calcium data
% zNonVis = interp1(nonVisTime, zNonVis, calciumTime, 'pchip');

% correlate calcium traces with running and pupil diameter
ind = ~isnan(zNonVis);
[rho, p] = corr(zCal(ind,:), zNonVis(ind)');

[rho, order] = sort(rho, 'descend');
zCal = zCal(:,order);
p = p(order);
neuronOrder = neuronIDs(order,:);
if ~isempty(isGad)
    isGad = isGad(order);
end

results.rho = rho;
results.p = p;
results.order = order;
results.orderedNeurons = neuronOrder;

if ~doPlot
    return
end

nonVisInds = find(p < signifP & rho >= rhoThreshold);
sets = ceil(length(nonVisInds) / 40);
setSize = ceil(length(nonVisInds) / sets);
handles = zeros(1, sets);
for iSet = 1:sets
    handles(iSet) = figure('Position', [1115 40 775 1075]);
    hold on
    inds = ((iSet-1)*setSize+1):min(iSet*setSize, length(nonVisInds));
    if ~isempty(stimOnsets)
        plot(repmat(reshape(stimOnsets,1,[]),2,1), ...
            repmat([-5*length(inds)-5;5],1,length(stimOnsets)), ...
            'Color', [0.8 0.8 0.8], 'LineWidth', 0.2)
    end
    plot(calciumTime, zNonVis, 'Color', [.3 .75 .93], 'LineWidth', 1.5)
    for iCell = 1:length(inds)
        c = 'k';
        if ~isempty(isGad)
            switch isGad(nonVisInds(inds(iCell)))
                case 1, c = colors(1,:);
                case 0, c = colors(2,:);
                case -1, c = colors(3,:);
            end
        end
        plot(calciumTime, zCal(:,nonVisInds(inds(iCell))) - 5*iCell, ...
            'Color', c)
        text(calciumTime(end), -5*iCell+2, sprintf('%.3f',rho(nonVisInds(inds(iCell)))))
    end
    xlim([0 calciumTime(end)*1.1])
    ylim([-5*length(inds)-5 5])
    set(gca, 'YTick', -5*length(inds):5:0, ...
        'YTickLabel', [flip(cellstr(num2str(neuronOrder(inds,:)))); {nonVisName}])
    xlabel('Time (in s)')
    ylabel('Neuron ID')
    title([nonVisName ' and neural responses - set ' num2str(iSet)])
end
figureHandles.traces = handles;
figureHandles.rho = figure;
bins = floor(min(rho)*50)/50:0.02:ceil(max(rho)*50)/50;
if isempty(isGad)
    n = hist(rho, bins);
    bar(bins, n, 'k')
else
    n1 = hist(rho(isGad == 1), bins);
    n2 = hist(rho(isGad == 0), bins);
    n3 = hist(rho(isGad == -1), bins);
    b = bar(bins, [n1' n2' n3'], 'stacked');
    for k = 1:3
        b(k).FaceColor = colors(k,:);
        b(k).EdgeColor = 'none';
    end
    legend('GAD+', 'GAD?', 'GAD-')
end
title(['Correlation coeff. ' nonVisName ' and neural response'])
xlabel('Rho')
ylabel('#Neurons')