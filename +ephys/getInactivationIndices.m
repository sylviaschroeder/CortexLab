function modulationIndices = getInactivationIndices(spikeTimes, clusters, ...
    eventTimes, window, stimIDs, laserOn)

clusterIDs = unique(clusters);
modulationIndices = NaN(size(clusterIDs));

laserOnTrials = ismember(stimIDs, find(laserOn==1));
laserOffTrials = ismember(stimIDs, find(laserOn==0));

for iClus = 1:length(clusterIDs)
    st = spikeTimes(clusters == clusterIDs(iClus));
    eventCounts = zeros(size(eventTimes));
    for iEvent = 1:length(eventTimes)
        eventCounts(iEvent) = sum(st >= eventTimes(iEvent)+window(1) & ...
            st <= eventTimes(iEvent)+window(2));
    end
    laserOnCount = sum(eventCounts(laserOnTrials));
    laserOffCount = sum(eventCounts(laserOffTrials));
    modulationIndices(iClus) = log(laserOnCount / laserOffCount);
end