function writeEphysToNPY(mouseName, date, tags)
% spikes.times (needs time alignment with master)
% spikes.clusters (from manual; don't bother about noise)
% spikes.depths (from drift map)
% spikes.amps
% spikes.templates
% clusters.ids
% cluster.groups (1 MUA, 2 good, 3 unsorted)
% templates.waveforms ([template x sample x subchan])
% templates.waveformsChannels ([template x subchan])
%   channels.localCoordinates ([chan x 2], x-coord, y-coord.

useDriftmap = true;
numWfChans = 10;

if nargin < 3
    tags = getEphysTags(mouseName, date);
end
root = getRootDir(mouseName, date);

destDir = fullfile(root, 'alf');

sp = loadAllKsDir(mouseName, date);

for tg = 1:length(tags)
    thisDest = fullfile(destDir, tags{tg});
    if ~exist(thisDest, 'dir'); mkdir(thisDest); end
    
    writeNPY(sp(tg).st, fullfile(thisDest, 'spikes.times.npy'));
    writeNPY(sp(tg).clu, fullfile(thisDest, 'spikes.clusters.npy'));
    writeNPY(sp(tg).spikeAmps, fullfile(thisDest, 'spikes.amps.npy'));
    writeNPY(sp(tg).spikeTemplates, fullfile(thisDest, 'spikes.templates.npy'));
    
    if ~useDriftmap
        writeNPY(sp(tg).spikeDepths, fullfile(thisDest, 'spikes.depths.npy'));
    end
    
    % include unsorted clusters into clusters.groups
    cids = sp(tg).cids(:);
    cgs = sp(tg).cgs(:);
    allCID = unique(sp(tg).clu);
    allCG = 3*ones(size(allCID));
    for c = 1:length(cids)
        allCG(allCID==cids(c))=cgs(c);
    end
    writeNPY(allCID, fullfile(thisDest, 'clusters.ids.npy'));
    writeNPY(allCG, fullfile(thisDest, 'clusters.groups.npy'));
    
    allTemps = unique(sp(tg).spikeTemplates);
    wf = NaN(length(allTemps), size(sp(tg).tempsUnW,2), 2*numWfChans+1);
    wfChans = NaN(length(allTemps), 2*numWfChans+1);
    [~, maxChans] = max(max(abs(sp(tg).temps(allTemps+1,:,:)),[],2),[],3);
    for temp = 1:length(allTemps)
        indChans = (-numWfChans:numWfChans) + maxChans(temp);
        valid = indChans > 0 & indChans <= size(sp(tg).temps,3);
        wf(temp,:,valid) = sp(tg).tempsUnW(allTemps(temp)+1,:,indChans(valid));
        wfChans(temp,valid) = indChans(valid);
    end
    writeNPY(wf, fullfile(thisDest, 'templates.waveforms.npy'));
    writeNPY(wfChans, fullfile(thisDest, 'templates.waveformsChannels.npy'));
    writeNPY([sp(tg).xcoords, sp(tg).ycoords], ...
        fullfile(thisDest, 'channels.localCoordinates.npy'));
end
clear sp wf wfChans

if useDriftmap
    for tg = 1:length(tags)
        thisDest = fullfile(destDir, tags{tg});
        ksDir = getKSdir(mouseName, date, tags{tg});
        [~,~, sd] = ksDriftmap(ksDir);
        writeNPY(sd, fullfile(thisDest, 'spikes.depths.npy'));
    end
end