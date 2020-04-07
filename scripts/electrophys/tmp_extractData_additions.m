%% Define datasets
db_ephys_driftingGratings;

%% Folders
expInfo = '\\ZSERVER.cortexlab.net\Data\expInfo';
resultFolder = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\electrophys';
expInfoFolder = '\\ZSERVER.cortexlab.net\Data\expInfo';

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
    
    % load data (time, spikeCounts, cellIDs, dephts, stimMatrix,
    % directions, blanks, stimTimes, laserOn, laserTimesRelative,
    % nSpikes, sampleRate, spikeWidths, stimSequence, timelineToEphys,
    % waveforms)
    load(fullfile(resultFolder, db(k).subject, db(k).date, ...
        sprintf('%02d_data_old.mat', db(k).exp)));
    
    if strcmp(db(k).subject, 'SS061')
        protocolFile = '\\ZSERVER.cortexlab.net\Data\trodes\M160426_SS061\20160511\1\Protocol.mat';
    else
        protocolFile = fullfile(expInfoFolder, db(k).subject, db(k).date, ...
            num2str(db(k).exp), sprintf('%s_%d_%s_parameters.mat', db(k).date, ...
            db(k).exp, db(k).subject));
    end
    
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
    stimSequence.seq = stimIDs(1:length(stimTimes.onset))';
    stimSequence.labels = 1:Protocol.npars;
    
    stimOn_rel = stimOn_rel(stimSequence.seq);
    stimOff_rel = stimOff_rel(stimSequence.seq);
    laserOn_rel = laserOn_rel(stimSequence.seq) - stimOn_rel;
    laserOff_rel = laserOff_rel(stimSequence.seq) - stimOff_rel;
    laser.onset = laserOn_rel;
    laser.offset = laserOff_rel;
    
    stimTimes.onset = stimTimes.onset(:) + stimOn_rel;
    stimTimes.offset = stimTimes.offset(:) + stimOff_rel;
    stimMatrix = ppbox.buildStimMatrix(stimSequence, stimTimes, time);
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

    %% Save data
    r.sampleRate = sampleRate;
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
    r.stimSequence = stimSequence;
    r.laserOn = laserOn;
    r.laserTimesRelative = laser;
    r.timelineToEphys = timelineToEphys;
    folder = fullfile(resultFolder, db(k).subject, db(k).date);
    if ~isdir(folder)
        mkdir(folder);
    end
    save(fullfile(folder, sprintf('%02d_data.mat', db(k).exp)), '-struct', 'r');
end