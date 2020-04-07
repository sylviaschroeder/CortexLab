%% Define datasets
db_ephys_driftingGratings;

%% Parameters
spikeRes = 0.005;
timeBeforeAfterStimuli = 10; %in sec

%% Folders
resultFolder = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\electrophys';

%% Constants
tags.stimDur = {'dur'};
tags.stimOnset = {'ton'};
tags.stimOffset = {'toff'};
tags.laserOn = {'tstart1'};
tags.laserOff = {'tend1'};
tags.ori = {'ori1','ori'};
tags.contrast = {'c1', 'cg'};
tags.amplitude = {'amp1'};

for k = 1:length(db)
    fprintf('Dataset %d of %d\n', k, length(db))
    expRoot = fileparts(dat.expPath(db(k).subject, db(k).date, 1, 'main', 'master'));
    if isempty(db(k).tag)
        rawDir = fullfile(expRoot, 'ephys');
    else
        rawDir = fullfile(expRoot, ['ephys_' db(k).tag]);
    end
    tlFile = dat.expFilePath(db(k).subject, db(k).date, db(k).tlExp, 'timeline', 'master');
%     tlFile = fullfile(expInfo, db(k).subject, date, db(k).tlExp, ...
%         sprintf('%s_1_%s_Timeline.dat', date, db(k).subject));
    protocolFile = dat.expFilePath(db(k).subject, db(k).date, db(k).exp, 'parameters', 'master');
    
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
    
    stimOnTL = readNPY(fullfile(alignDir, ...
        sprintf('mpep_%d_onsets_in_timeline_%d.npy', db(k).exp, db(k).tlExp)));
    
    stimOffTL = readNPY(fullfile(alignDir, ...
        sprintf('mpep_%d_offsets_in_timeline_%d.npy', db(k).exp, db(k).tlExp)));
    
    stimOn = applyCorrection(stimOnTL, bTLtoMaster);
    stimOff = applyCorrection(stimOffTL, bTLtoMaster);
    
    %% Load protocol information
    data = load(protocolFile);
    if isfield(data, 'Protocol')
        Protocol = data.Protocol;
    else
        Protocol = data.parameters.Protocol;
    end
    
    stimDur = round(Protocol.pars(strcmp(Protocol.parnames, ...
        tags.stimDur{1}),:)' ./ 10 .* 60) ./ 60; % 60 Hz monitor frame rate
    stimOn_rel = round(Protocol.pars(strcmp(Protocol.parnames, ...
        tags.stimOnset{1}),:)' ./ 1000 .* 60) ./ 60; 
    stimOff_rel = round(Protocol.pars(strcmp(Protocol.parnames, ...
        tags.stimOffset{1}),:)' ./ 1000 .* 60) ./ 60 - stimDur;
    laserOn_rel = round(Protocol.pars(strcmp(Protocol.parnames, ...
        tags.laserOn{1}),:)' ./ 1000 .* 60) ./ 60; 
    laserOff_rel = round(Protocol.pars(strcmp(Protocol.parnames, ...
        tags.laserOff{1}),:)' ./ 1000 .* 60) ./ 60 - stimDur;

    stimIDs = zeros(1, numel(Protocol.seqnums));
    for q = 1:size(Protocol.seqnums,1)
        stimIDs(Protocol.seqnums(q,:)) = q;
    end
    stimSequence.seq = stimIDs(1:length(stimOn))';
    stimSequence.labels = 1:Protocol.npars;
    
    stimOn_rel = stimOn_rel(stimSequence.seq);
    stimOff_rel = stimOff_rel(stimSequence.seq);
    laserOn_rel = laserOn_rel(stimSequence.seq) - stimOn_rel;
    laserOff_rel = laserOff_rel(stimSequence.seq) - stimOff_rel;
    laser.onset = laserOn_rel;
    laser.offset = laserOff_rel;
    
    time = (round((stimOn(1)-timeBeforeAfterStimuli)/spikeRes)*spikeRes : ...
        spikeRes : (stimOff(end)+timeBeforeAfterStimuli))';
    stimTimes.onset = stimOn(:) + stimOn_rel;
    stimTimes.offset = stimOff(:) + stimOff_rel;
    stimMatrix = ppbox.buildStimMatrix(stimSequence, stimTimes, time);
    tag = [];
    for t = 1:length(tags.contrast)
        tag = strcmp(tags.contrast{t}, Protocol.parnames);
        if any(tag)
            tag = tags.contrast{t};
            break
        end
    end
    if isempty(tag)
        disp('ERROR: stimulus parameter for contrast not found.')
        return
    end
    blanks = find(Protocol.pars(strcmp(Protocol.parnames, tag),:)' == 0);
    tag = [];
    for t = 1:length(tags.ori)
        tag = strcmp(tags.ori{t}, Protocol.parnames);
        if any(tag)
            tag = tags.ori{t};
            break
        end
    end
    if isempty(tag)
        disp('ERROR: stimulus parameter for orientation not found.')
        return
    end
    directions = Protocol.pars(strcmp(Protocol.parnames, tag),:)';
    directions(blanks) = NaN;
    tag = [];
    for t = 1:length(tags.amplitude)
        tag = strcmp(tags.amplitude{t}, Protocol.parnames);
        if any(tag)
            tag = tags.amplitude{t};
            break
        end
    end
    if isempty(tag)
        disp('ERROR: stimulus parameter for laser amplitude not found.')
        return
    end
    laserOn = Protocol.pars(strcmp(Protocol.parnames, 'amp1'),:)' > 0;
    
    %% Load spikes
    s = loadAllKsDir(db(k).subject, db(k).date);
    goodClusters = s.cids(s.cgs == 2);
    ksDir = getKSdir(db(k).subject, db(k).date, s.name);
    [~,~,sd] = ksDriftmap(ksDir);
    
    spikeCounts = NaN(length(time), length(goodClusters));
    depths = NaN(1, length(goodClusters));
    spikeWidths = NaN(1, length(goodClusters));
    waveforms = cell(1, length(goodClusters));
    nSpikes = cell(1, length(goodClusters));
    notSC = [];
    for iCell = 1:length(goodClusters)
        spikeInds = s.clu==goodClusters(iCell) & s.st>=time(1) & s.st<=time(end);
        depths(iCell) = mean(sd(spikeInds));
        tmplt = unique(s.spikeTemplates(spikeInds));
        wdth = 0;
        waveform = zeros(size(s.temps,2),length(tmplt));
        n = zeros(1,length(tmplt));
        for t = 1:length(tmplt)
            wdth = wdth + s.tempDur(tmplt(t)+1) * ...
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
    r.stimMatrix = stimMatrix;
    r.directions = directions;
    r.blanks = blanks;
    r.stimTimes = stimTimes;
    r.laserOn = laserOn;
    r.timelineToEphys = bTLtoMaster;
    folder = fullfile(resultFolder, db(k).subject, db(k).date);
    if ~isfolder(folder)
        mkdir(folder);
    end
    save(fullfile(folder, sprintf('%02d_data.mat', db(k).exp)), '-struct', 'r');
end