function plotOriDiffsVsDistance(orientationDifferences, distances, ...
    titleLine, permutedDifferences)

clusterDist = 25;

if nargin < 3 || isempty(titleLine)
    titleLine = 'orientation';
end
if nargin < 4
    permutedDifferences = [];
end

if ~iscell(orientationDifferences) && ~iscell(distances) && ...
        ~isvector(orientationDifferences) && ~isvector(distances)
    orientationDifferences = {orientationDifferences};
    distances = {distances};
elseif (~iscell(orientationDifferences) && iscell(distances)) || ...
        (iscell(orientationDifferences) && ~iscell(distances))
    display('orientationDifferences and distances have to have the same format!')
    return
end

if iscell(orientationDifferences)
    diffs = [];
    dist = [];
    for set = 1:length(orientationDifferences)
        numNeurons = size(orientationDifferences{set},1);
        ind = ones(numNeurons);
        ind = full(spdiags(ind, 1:numNeurons, numNeurons, numNeurons));
        diffs = [diffs; orientationDifferences{set}(ind == 1)];
        dist = [dist; distances{set}(ind == 1)];
    end
else
    diffs = orientationDifferences;
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
sems = reshape([semDiffs, semDiffs]', [], 1);
permLow = reshape([permPercentiles(:,1),permPercentiles(:,1)]', [], 1);
permHigh = reshape([permPercentiles(:,3),permPercentiles(:,3)]', [], 1);
permMed = reshape([permPercentiles(:,2),permPercentiles(:,2)]', [], 1);
clusterEdges = reshape([clusterEdges; clusterEdges], [], 1);
clusterEdges([1 end]) = [];
m = means;
m(isnan(m)) = 0;
s = sems;
s(isnan(s)) = 0;
% fill([clusterEdges; flip(clusterEdges)], [m + s; flip(m - s)], ...
%     'k', 'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
plot(clusterEdges, means, 'r', 'LineWidth', 1)
plot(clusterEdges, means+sems, 'r:')
plot(clusterEdges, means-sems, 'r:')

pl = permLow;
pl(isnan(pl)) = 0;
ph = permHigh;
ph(isnan(ph)) = 0;
fill([clusterEdges; flip(clusterEdges)], [pl; ...
    flip(ph)], ...
    'k', 'FaceColor', 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
plot(clusterEdges, permMed, 'b', 'LineWidth', 1)

xlim([0 clusterEdges(find(~isnan(means),1,'last'))])
xlabel('Distance between neurons')
ylabel(sprintf('Difference in preferred %s', titleLine))