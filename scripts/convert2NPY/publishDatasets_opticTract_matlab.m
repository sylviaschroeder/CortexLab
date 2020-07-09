db_ephys_opticTract

%% Folders
folderBase = '\\ZUBJECTS.cortexlab.net\Subjects';
folderSave = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\DataToPublish\opticTract';

%% Collect data

for k = 1:length(db)
    fprintf('Dataset %d of %d\n', k, length(db))
    data = [];
    subject = db(k).subject;
    ind = strfind(subject, 'SS');
    subject = subject(ind:end);
    data.subject = subject;
    data.date = db(k).date;
    
    alignDir = fullfile(folderBase, db(k).subject, db(k).date, 'alignments');
    sp = loadAllKsDir(db(k).subject, db(k).date);
    probe = db(k).OTprobe;
    rawF = fullfile(folderBase, db(k).subject, db(k).date, ...
        ['ephys_' sp(probe).name], 'sorting');
    if isfolder(rawF)
        [~, ~, spikeDepths] = ksDriftmap(rawF);
    end
    templates = findTempForEachClu(sp(probe).clu, sp(probe).spikeTemplates);
    
    units = db(k).OTunits{1}(db(k).OTgood{1}==1);
    ind = ismember(sp(probe).clu, units);
    
    data.sampleRate = sp(probe).sample_rate;
    data.spikeTimes = sp(probe).st(ind);
    data.clusters = sp(probe).clu(ind);
    data.spikeAmps = sp(probe).spikeAmps(ind);
    data.spikeDepths = spikeDepths(ind);
    data.xcoords = sp(probe).xcoords;
    data.ycoords = sp(probe).ycoords;
    
    data.units.cids = units;
    data.units.waveform = NaN(size(sp(probe).tempsUnW,2), ...
        size(sp(probe).tempsUnW,3), length(units));
    data.units.goodTimes = db(k).OTtimes{1}(db(k).OTgood{1}==1);
    for iCell = 1:length(units)
        data.units.waveform(:,:,iCell) = squeeze(sp(probe). ...
            tempsUnW(templates(units(iCell)+1)+1,:,:));
    end
    
    [expNums, ~, ~, ~, ~, tl, hasTimeline] = ...
        dat.whichExpNums(db(k).subject, db(k).date);
    TLexp = expNums(hasTimeline);
    TLexp = TLexp(end);
    tl = tl{end};
    bTLtoMaster = readNPY(fullfile(alignDir, ...
        sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, db(k).probeNames{1})));
    tlTime = applyCorrection(tl.rawDAQTimestamps, bTLtoMaster);
    
    rotary = double(tl.rawDAQData(:, strcmp({tl.hw.inputs.name}, ...
        'rotaryEncoder')));
    runSpeed = nonVis.getRunningSpeed_wheel(rotary, tlTime, 0);
    runTime = runSpeed.t;
    runSpeed = runSpeed.total;
    cmPerUnit = 2*pi * 8.75 / (4 * 1024);
    runSpeed = runSpeed * cmPerUnit;
    data.runningSpeed = runSpeed;
    data.time_runningSpeed = runTime;
    
    data.dark.time = db(k).darkTime;
    
    if ~isempty(db(k).expFlicker)
        d = load(fullfile(folderBase, db(k).subject, db(k).date, ...
            num2str(db(k).expFlicker), sprintf('%s_%d_%s_parameters.mat', ...
            db(k).date, db(k).expFlicker, db(k).subject)));
        parsFlicker = d.parameters.Protocol;
        d = load(fullfile(folderBase, db(k).subject, db(k).date, ...
            num2str(db(k).expFlicker), sprintf('%s_%d_%s_hardwareInfo.mat', ...
            db(k).date, db(k).expFlicker, db(k).subject)));
        monitorRefreshRate = d.myScreenInfo.FrameRate;
        flickerDurs = parsFlicker.pars(strcmp(parsFlicker.parnames, ...
            'nfr'),:) / monitorRefreshRate;
        
        stimOnTL = readNPY(fullfile(alignDir, ...
            sprintf('mpep_%d_onsets_in_timeline_%d.npy', db(k).expFlicker, TLexp)));
        stimOffTL = readNPY(fullfile(alignDir, ...
            sprintf('mpep_%d_offsets_in_timeline_%d.npy', db(k).expFlicker, TLexp)));
        flickerWhite = applyCorrection(stimOnTL, bTLtoMaster);
        flickerBlack = applyCorrection(stimOffTL, bTLtoMaster);
        flickerTimes = cell(size(parsFlicker.seqnums));
        if length(flickerWhite) == size(parsFlicker.pars,2)*parsFlicker.nrepeats % only stim onset was recorded, not every flick
            for stim = 1:size(flickerTimes,1)
                for trial = 1:size(flickerTimes,2)
                    ind = parsFlicker.seqnums(stim,trial);
                    flickerTimes{stim,trial} = (flickerWhite(ind) : ...
                        flickerDurs(stim) : flickerBlack(ind))' + ...
                        [0 flickerDurs(stim)/2];
                    if flickerBlack(ind) < flickerTimes{stim,trial}(end,2)
                        flickerTimes{stim,trial}(end,:) = [];
                    end
                end
            end
            flickerAll = [flickerWhite(1) flickerBlack(end)];
        else % each flicker was recorded
            flickerAll = reshape([flickerWhite, flickerBlack]', [], 1);
            longestFlicker = max(parsFlicker.pars( ...
                strcmp(parsFlicker.parnames, 'nfr'),:)) / monitorRefreshRate;
            durs = diff(flickerAll);
            flickerStimStartInds = [1; find(durs > ...
                longestFlicker * 1.5) + 1; length(flickerAll)+1];
            for stim = 1:size(flickerTimes,1)
                for trial = 1:size(flickerTimes,2)
                    ind = parsFlicker.seqnums(stim,trial);
                    flickerTimes{stim,trial} = reshape(flickerAll( ...
                        flickerStimStartInds(ind) : ...
                        flickerStimStartInds(ind+1)-1), 2, [])';
                end
            end
        end
        data.flicker.times = flickerTimes;
    end
    if ~isempty(db(k).expNoise)
        d = load(fullfile(folderBase, db(k).subject, db(k).date, ...
            num2str(db(k).expNoise), sprintf('%s_%d_%s_parameters.mat', ...
            db(k).date, db(k).expNoise, db(k).subject)));
        parsNoise = d.parameters.Protocol;
        stimFile = str2func(strtok(parsNoise.xfile, '.'));
        % load myScreenInfo
        load(fullfile(folderBase, db(k).subject, db(k).date, ...
            num2str(db(k).expNoise), sprintf('%s_%d_%s_hardwareInfo.mat', ...
            db(k).date, db(k).expNoise, db(k).subject)));
        myScreenInfo.windowPtr = NaN;
        stimOnTL = readNPY(fullfile(alignDir, ...
            sprintf('mpep_%d_onsets_in_timeline_%d.npy', db(k).expNoise, TLexp)));
        noiseOn = applyCorrection(stimOnTL, bTLtoMaster);
        % call x-file to create stimuli
        SS = stimFile(myScreenInfo, parsNoise.pars);
        stimFrames = cat(3, SS.ImageTextures{:});
        framesPerImage = parsNoise.pars(6,1);
        frameTimes = (0 : size(stimFrames, 3)-1) * framesPerImage / myScreenInfo.FrameRate;
        noiseFrameTimes = (frameTimes + noiseOn)';
        noisePosition = parsNoise.pars(2:5)./10;
        
        data.noise.onTimes = noiseFrameTimes;
        data.noise.position = noisePosition;
        data.noise.sequence = repmat((1:size(stimFrames,3))', length(noiseOn), 1);
        data.noise.frames = stimFrames;
    end
    if ~isempty(db(k).expOri)
        d = load(fullfile(folderBase, db(k).subject, db(k).date, ...
            num2str(db(k).expOri), sprintf('%s_%d_%s_parameters.mat', ...
            db(k).date, db(k).expOri, db(k).subject)));
        parsGratings = d.parameters.Protocol;
        stimOnTL = readNPY(fullfile(alignDir, ...
            sprintf('mpep_%d_onsets_in_timeline_%d.npy', db(k).expOri, TLexp)));
        stimOffTL = readNPY(fullfile(alignDir, ...
            sprintf('mpep_%d_offsets_in_timeline_%d.npy', db(k).expOri, TLexp)));
        stimOn = applyCorrection(stimOnTL, bTLtoMaster);
        stimOff = applyCorrection(stimOffTL, bTLtoMaster);
        stimIDs = repmat((1:parsGratings.npfilestimuli)',parsGratings.nrepeats,1);
        [~,order] = sort(parsGratings.seqnums(:));
        stimSeq = stimIDs(order);
        blank = parsGratings.pars(15,:) == 0;
        directions = parsGratings.pars(6,:);
        directions(blank) = NaN;
        
        data.gratings.onTimes = stimOn;
        data.gratings.offTimes = stimOff;
        data.gratings.sequence = stimSeq;
        data.gratings.stimuli.directions = directions;
    end
    
    save(fullfile(folderSave, sprintf('%s_%s.mat', subject, ...
        db(k).date)), 'data')
end