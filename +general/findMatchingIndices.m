function ind = findMatchingIndices(origIndices, origTime, newTime)

origSR = median(diff(origTime));
newSR = median(diff(newTime));
if islogical(origIndices)
    origIndices = find(origIndices);
end

if origSR <= newSR
    n = hist(origTime(origIndices), newTime);
    ind = find(n > 0);
else
    n = cumsum(hist(newTime, origTime));
    m = [1 n(1:end-1)+1];
    ind = [];
    for k = origIndices(:)'
        ind = [ind, m(k):n(k)];
    end
end