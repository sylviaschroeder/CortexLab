%% Define datasets
db_ephys_darkScreens;

%% Parameters
spikeRes = 0.005;
timeBeforeAfterStimuli = 10; %in sec

%% Folders
ksRaw = '\\basket.cortexlab.net\data\nick\';
ksSorted = '\\basket.cortexlab.net\data\sylvia\';
expInfo = '\\ZSERVER.cortexlab.net\Data\expInfo';
resultFolder = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\electrophys';

for k = 1:length(db)
    fprintf('Dataset %d of %d\n', k, length(db))
    expRoot = fileparts(dat.expPath(db(k).subject, db(k).date, 1, 'main', 'master'));
    rawDir = fullfile(expRoot, ['ephys_' db(k).tag]);
    alignDir = fullfile(expRoot, 'alignments');
    
    %% Load sync information
    if ~strcmp(db(k).tag, db(k).masterTimebase)
        bEphysToMaster = readNPY(fullfile(alignDir, ...
            sprintf('correct_ephys_%s_to_ephys_%s.npy', db(k).tag, db(k).masterTimebase)));
    else % this one is master, so use a dummy conversion
        bEphysToMaster = [1; 0];
    end
    
    bTLtoMaster = readNPY(fullfile(alignDir, ...
        sprintf('correct_timeline_%d_to_ephys_%s.npy', db(k).tlExp, db(k).masterTimebase)));
    
    %% Load protocol information
    time = (round(db(k).time(1)/spikeRes)*spikeRes : ...
        spikeRes : db(k).time(2))';
    
    %% Load spikes
    try
        s = loadKSdir(fullfile(ksRaw, db(k).subject, db(k).date, ['ephys_' db(k).tag]));
    catch
        s = loadKSdir(fullfile(ksRaw, db(k).subject, db(k).date));
    end
    s.st = applyCorrection(s.st, bEphysToMaster);
    
    try
        s.clu = readNPY(fullfile(ksSorted, db(k).subject, db(k).date, ...
            ['ephys_' db(k).tag], 'spike_clusters.npy'));
        [s.cids, s.cgs] = readClusterGroupsCSV(fullfile(ksSorted, db(k).subject, ...
            db(k).date, ['ephys_' db(k).tag], 'cluster_group.tsv'));
    catch
        s.clu = readNPY(fullfile(ksSorted, db(k).subject, db(k).date, ...
            'spike_clusters.npy'));
        [s.cids, s.cgs] = readClusterGroupsCSV(fullfile(ksSorted, db(k).subject, ...
            db(k).date, 'cluster_groups.csv'));
    end
    if length(s.clu) > length(s.st) % noise clusters were removed in cids but not in clu
        noiseSpikes = ismember(s.clu, s.cids(s.cgs==0));
        s.clu(noiseSpikes) = [];
    end
    
    [spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] = ...
        templatePositionsAmplitudes(s.temps, s.winv, s.ycoords, s.spikeTemplates, s.tempScalingAmps);
    
    goodClusters = s.cids(s.cgs == 2);
    
    spikeCounts = NaN(length(time), length(goodClusters));
    depths = NaN(1, length(goodClusters));
    spikeWidths = NaN(1, length(goodClusters));
    waveforms = cell(1, length(goodClusters));
    nSpikes = cell(1, length(goodClusters));
    notSC = [];
    for iCell = 1:length(goodClusters)
        spikeInds = s.clu==goodClusters(iCell) & s.st>=time(1) & s.st<=time(end);
        depths(iCell) = mean(spikeDepths(spikeInds));
        tmplt = unique(s.spikeTemplates(spikeInds));
        wdth = 0;
        waveform = zeros(size(s.temps,2),length(tmplt));
        n = zeros(1,length(tmplt));
        for t = 1:length(tmplt)
            wdth = wdth + templateDuration(tmplt(t)+1) * ...
                (sum(s.spikeTemplates == tmplt(t) & spikeInds) / sum(spikeInds));
            tempUnW = squeeze(s.temps(tmplt(t)+1,:,:))*s.winv;
            [~,ind] = max(max(tempUnW, [], 1), [], 2);
            waveform(:,t) = tempUnW(:,ind);
            n(t) = sum(s.spikeTemplates == tmplt(t) & spikeInds);
        end
        spikeWidths(iCell) = wdth ./ s.sample_rate * 1000;
        waveforms{iCell} = waveform;
        nSpikes{iCell} = n;
        if sum(spikeInds)==0 || depths(iCell)<db(k).yPosRange(1) || ...
                depths(iCell)>db(k).yPosRange(2)
            notSC(end+1) = iCell;
            continue
        end
        sc = hist(s.st(spikeInds), time);
        spikeCounts(:, iCell) = sc;
    end
    cellIDs = goodClusters(setdiff(1:length(goodClusters), notSC));
    spikeCounts(:,notSC) = [];
    depths(notSC) = [];
    spikeWidths(notSC) = [];
    waveforms(notSC) = [];
    nSpikes(notSC) = [];
    
%     r = load(fullfile(resultFolder, db(k).subject, db(k).date, ...
%         sprintf('%02d_data.mat', db(k).exp)));
    
    r.sampleRate = s.sample_rate;
    r.time = time;
    r.spikeCounts = spikeCounts;
    r.cellIDs = cellIDs;
    r.depths = depths;
    r.spikeWidths = spikeWidths;
    r.waveforms = waveforms;
    r.nSpikes = nSpikes;
    r.timelineToEphys = bTLtoMaster;
    folder = fullfile(resultFolder, db(k).subject, db(k).date);
    if ~isdir(folder)
        mkdir(folder);
    end
    save(fullfile(folder, 'dark_data.mat'), '-struct', 'r');
end