db_ephys_V1

%% Folders
folderBase = '\\ZUBJECTS.cortexlab.net\Subjects';
folderSave = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\DataToPublish\v1 neurons ephys';

%% Constants
tags.stimDur = {'dur'};
tags.stimOnset = {'ton'};
tags.stimOffset = {'toff'};
tags.laserOn = {'tstart1'};
tags.laserOff = {'tend1'};
tags.ori = {'ori1','ori'};
tags.contrast = {'c1', 'cg'};
tags.amplitude = {'amp1'};

%% Collect data
for k = 1:2
    fprintf('Dataset %d of %d\n', k, length(db))
    data = [];
    subject = db(k).subject;
    folderSession = fullfile(folderSave, subject, db(k).date, '001');
    if ~isfolder(folderSession)
        mkdir(folderSession)
    end
    alignDir = fullfile(folderBase, db(k).subject, db(k).date, 'alignments');
    
    [expNums, ~, ~, ~, ~, tl, hasTimeline] = ...
        dat.whichExpNums(db(k).subject, db(k).date);
    TLexp = expNums(hasTimeline);
    TLexp = TLexp(end);
    tl = tl{end};
    bTLtoMaster = readNPY(fullfile(alignDir, ...
        sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, db(k).masterTimebase)));
    tlTime = applyCorrection(tl.rawDAQTimestamps, bTLtoMaster);
    
    % grating parameters
    d = load(fullfile(folderBase, db(k).subject, db(k).date, ...
        num2str(db(k).exp), sprintf('%s_%d_%s_parameters.mat', ...
        db(k).date, db(k).exp, db(k).subject)));
    parsGratings = d.parameters.Protocol;
    seqnums = parsGratings.seqnums;
    seqnums(:,db(k).excludeReps) = [];
    validTrials = sort(reshape(seqnums,[],1));
    stimOnTL = readNPY(fullfile(alignDir, ...
        sprintf('mpep_%d_onsets_in_timeline_%d.npy', db(k).exp, TLexp)));
    stimOffTL = readNPY(fullfile(alignDir, ...
        sprintf('mpep_%d_offsets_in_timeline_%d.npy', db(k).exp, TLexp)));
    stimOn = applyCorrection(stimOnTL(validTrials), bTLtoMaster);
    stimOff = applyCorrection(stimOffTL(validTrials), bTLtoMaster);
    stimIDs = repmat((1:parsGratings.npfilestimuli)',parsGratings.nrepeats,1);
    stimIDs = stimIDs(validTrials);
    [~,order] = sort(seqnums(:));
    stimSeq = stimIDs(order);
    stimDur = round(parsGratings.pars(strcmp(parsGratings.parnames, ...
        tags.stimDur{1}),:)' ./ 10 .* 60) ./ 60; % 60 Hz monitor frame rate
    stimOn_rel = round(parsGratings.pars(strcmp(parsGratings.parnames, ...
        tags.stimOnset{1}),:)' ./ 1000 .* 60) ./ 60;
    stimOff_rel = round(parsGratings.pars(strcmp(parsGratings.parnames, ...
        tags.stimOffset{1}),:)' ./ 1000 .* 60) ./ 60 - stimDur;
    laserOn_rel = round(parsGratings.pars(strcmp(parsGratings.parnames, ...
        tags.laserOn{1}),:)' ./ 1000 .* 60) ./ 60; 
    laserOff_rel = round(parsGratings.pars(strcmp(parsGratings.parnames, ...
        tags.laserOff{1}),:)' ./ 1000 .* 60) ./ 60 - stimDur;
    laserOn_rel = laserOn_rel - stimOn_rel;
    laserOff_rel = laserOff_rel - stimOff_rel;
    stimOn_rel = stimOn_rel(stimSeq);
    stimOff_rel = stimOff_rel(stimSeq);
    stimOn = stimOn + stimOn_rel;
    stimOff = stimOff + stimOff_rel;
    stimDur = mean(stimOff - stimOn);
    tag = tags.contrast{ismember(tags.contrast, parsGratings.parnames)};
    blank = parsGratings.pars(strcmp(parsGratings.parnames,tag),:) == 0;
    tag = tags.ori{ismember(tags.ori, parsGratings.parnames)};
    directions = parsGratings.pars(strcmp(parsGratings.parnames,tag),:);
    directions(blank) = NaN;
    laserOn = parsGratings.pars(strcmp(parsGratings.parnames,tags.amplitude),:) > 0;
    
    writeNPY([stimOn stimOff], fullfile(folderSession, ...
        '_ss_grating.intervals.npy'));
    writeNPY(stimSeq, fullfile(folderSession, ...
        '_ss_grating._ss_gratingID.npy'));
    writeNPY(directions', fullfile(folderSession, ...
        '_ss_gratingID.directions.npy'));
    writeNPY(laserOn', fullfile(folderSession, ...
        '_ss_gratingID.laserOn.npy'));
    writeNPY(laserOn_rel, fullfile(folderSession, ...
        '_ss_gratingID.laserOnTime.npy'));
    writeNPY(laserOff_rel, fullfile(folderSession, ...
        '_ss_gratingID.laserOffTime.npy'));
    
    % spike data
    sp = loadAllKsDir(db(k).subject, db(k).date);
    if ~isempty(db(k).tag)
        probe = find(strcmp({sp.name}, db(k).tag));
        name = ['ephys_' sp(probe).name];
    else
        probe = 1;
        name = 'ephys';
    end
    rawF = fullfile(folderBase, db(k).subject, db(k).date, ...
        name, 'sorting');
    if isfolder(rawF)
        [~, ~, spikeDepths] = ksDriftmap(rawF);
    end
    templates = findTempForEachClu(sp(probe).clu, sp(probe).spikeTemplates);
    
    units = sp(probe).cids(:);
    depths = NaN(length(units),1);
    for iCell = 1:length(units)
        depths(iCell) = nanmean(spikeDepths(sp(probe).clu == units(iCell)));
    end
    included = depths>=db(k).yPosRange(1) & depths<=db(k).yPosRange(2);
    units = units(included);
    depths = depths(included);
    ind = ismember(sp(probe).clu, units);
    waveforms = NaN(length(units), size(sp(probe).tempsUnW,2), ...
        size(sp(probe).tempsUnW,3));
    for iCell = 1:length(units)
        waveforms(iCell,:,:) = squeeze(sp(probe). ...
            tempsUnW(templates(units(iCell)+1)+1,:,:));
    end
    cluster_id = units;
    group = cell(length(units),1);
    group(sp(probe).cgs(included)==1) = {'mua'};
    group(sp(probe).cgs(included)==2) = {'good'};
    groups = table(cluster_id, group);
    
    writetable(groups, fullfile(folderSession, ...
        'cluster.groups.csv'))
    writeNPY(sp(probe).st(ind), fullfile(folderSession, ...
        'spike.times.npy'));
    writeNPY(sp(probe).spikeAmps(ind), fullfile(folderSession, ...
        'spike.amps.npy'));
    writeNPY(sp(probe).clu(ind), fullfile(folderSession, ...
        'spike.clusters.npy'));
    writeNPY(spikeDepths(ind), fullfile(folderSession, ...
        'spike.depths.npy'));
    writematrix(units, fullfile(folderSession, 'clusters.uuids.csv'));
    writeNPY(waveforms, fullfile(folderSession, ...
        'clusters.waveforms.npy'));
    writeNPY([sp(probe).xcoords sp(probe).ycoords], fullfile(folderSession, ...
        'channels.localCoordinates.npy'));
    writeNPY(sp(probe).sample_rate, fullfile(folderSession, ...
        'probe._ss_sampleRate.npy'));
    writeNPY(db(k).yPosRange, fullfile(folderSession, ...
        'probe.v1Depth.npy'));
end