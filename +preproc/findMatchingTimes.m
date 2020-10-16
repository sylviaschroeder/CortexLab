function [indicesA, indicesB, timeDiffs] = findMatchingTimes(timesA, timesB)

timesA = timesA(:);
timesB = timesB(:);

differences = abs(timesA - timesB');
[minDiffs, bestMatches] = min(differences, [], 2);

indicesB = unique(bestMatches);
indicesA = NaN(size(indicesB));
for j = 1:length(indicesB)
    ind = find(bestMatches == indicesB(j));
    [~, best] = min(minDiffs(ind));
    indicesA(j) = ind(best);
end
timeDiffs = timesA(indicesA) - timesB(indicesB);