function plotRFDistsVsBrainDists(RFdists, neuralDists)

clusterDist = 25;


RFd = [];
neuralD = [];
for set = 1:length(RFdists)
    numNeurons = size(RFdists{set},1);
    ind = ones(numNeurons);
    ind = full(spdiags(ind, 1:numNeurons, numNeurons, numNeurons));
    RFd = [RFd; RFdists{set}(ind == 1)];
    neuralD = [neuralD; neuralDists{set}(ind == 1)];
end

maxDist = max(neuralD);
clusterEdges = 0:clusterDist:maxDist+clusterDist;
meanDiffs = NaN(length(clusterEdges)-1, 1);
semDiffs = NaN(length(clusterEdges)-1, 1);
for clust = 1:length(meanDiffs)
    d = RFd(neuralD > clusterEdges(clust) & neuralD < clusterEdges(clust+1));
    meanDiffs(clust) = nanmean(d);
    semDiffs(clust) = nanstd(d) / sqrt(sum(~isnan(d)));
end

figure
plot(neuralD, RFd, 'k.')
hold on
means = reshape([meanDiffs, meanDiffs]', [], 1);
sems = reshape([semDiffs, semDiffs]', [], 1);
clusterEdges = reshape([clusterEdges; clusterEdges], [], 1);
clusterEdges([1 end]) = [];
m = means;
m(isnan(m)) = 0;
s = sems;
s(isnan(s)) = 0;
fill([clusterEdges; flip(clusterEdges)], [m + s; flip(m - s)], ...
    'k', 'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
plot(clusterEdges, means, 'r', 'LineWidth', 2)

xlim([0 clusterEdges(find(~isnan(means),1,'last'))])
xlabel('Distance between neurons')
ylabel('Distance between RFs')