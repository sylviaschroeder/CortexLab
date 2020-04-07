function plotCorrDiffsVsDistance(correlationDifferences, distances, ...
    permutedDifferences)

clusterDist = 25;
if nargin < 3
    permutedDifferences = [];
end

if ~iscell(correlationDifferences) && ~iscell(distances) && ...
        ~isvector(correlationDifferences) && ~isvector(distances)
    correlationDifferences = {correlationDifferences};
    distances = {distances};
elseif (~iscell(correlationDifferences) && iscell(distances)) || ...
        (iscell(correlationDifferences) && ~iscell(distances))
    display('correlationDifferences and distances have to have the same format!')
    return
end

if iscell(correlationDifferences)
    diffs = [];
    dist = [];
    for set = 1:length(correlationDifferences)
        numNeurons = size(correlationDifferences{set},1);
        ind = ones(numNeurons);
        ind = full(spdiags(ind, 1:numNeurons, numNeurons, numNeurons));
        diffs = [diffs; correlationDifferences{set}(ind == 1)];
        dist = [dist; distances{set}(ind == 1)];
    end
else
    diffs = correlationDifferences;
    dist = distances;
end

maxDist = max(dist);
clusterEdges = 0:clusterDist:maxDist+clusterDist;
meanDiffs = NaN(length(clusterEdges)-1, 1);
semDiffs = NaN(length(clusterEdges)-1, 1);
permMeanDiffs = NaN(length(clusterEdges)-1, size(permutedDifferences,2));
for clust = 1:length(meanDiffs)
    ind = dist > clusterEdges(clust) & dist < clusterEdges(clust+1);
    d = diffs(ind);
    meanDiffs(clust) = nanmean(d);
    semDiffs(clust) = nanstd(d) / sqrt(sum(~isnan(d)));
    permMeanDiffs(clust,:) = nanmean(permutedDifferences(ind,:),1);
end
permPercentiles = prctile(permMeanDiffs,[2.5 50 97.5],2);

figure
plot(dist, diffs, 'k.')
hold on
means = reshape([meanDiffs, meanDiffs]', [], 1);
permLow = reshape([permPercentiles(:,1),permPercentiles(:,1)]', [], 1);
permHigh = reshape([permPercentiles(:,3),permPercentiles(:,3)]', [], 1);
permMed = reshape([permPercentiles(:,2),permPercentiles(:,2)]', [], 1);
clusterEdges = reshape([clusterEdges; clusterEdges], [], 1);
clusterEdges([1 end]) = [];
m = means;
m(isnan(m)) = 0;
plot(clusterEdges, means, 'r', 'LineWidth', 1)

pl = permLow;
pl(isnan(pl)) = 0;
ph = permHigh;
ph(isnan(ph)) = 0;
fill([clusterEdges; flip(clusterEdges)], [pl; flip(ph)], ...
    'k', 'FaceColor', 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
plot(clusterEdges, permMed, 'b', 'LineWidth', 1)

xlim([0 clusterEdges(find(~isnan(means),1,'last'))])
xlabel('Distance between neurons')
ylabel('Difference in corrrlation strength')