db_ephys_opticTract

%% Folders
folderBase = '\\ZUBJECTS.cortexlab.net\Subjects';
folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS';
folderSave = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\DataToPublish\opticTract';

%% Parameters
% to determine whether animal was running or not at spike times during
% darkness
runThreshold = 1;
runningSigma = 0.25;
sigma = 1;
binSizeRun = 1/7.5;
binEdges = db(k).darkTime(1) : binSizeRun : db(k).darkTime(2);
timeBins = binEdges(1:end-1) + binSizeRun/2;
sig = round(sigma / binSizeRun);
win = normpdf(-5*sig : 5*sig, 0, sig);

%% Collect data
d = load(fullfile(folderResults, 'receptiveFields\axons\receptiveFields.mat'));
RFs = d.RFs;
d = load(fullfile(folderResults, 'OpticTract\runningCorrelation_dark\correlations_runningFiltered.mat'));
corrs = d.corrs;

for k = 1:length(db)
    fprintf('Dataset %d of %d\n', k, length(db))
    subject = db(k).subject;
    ind = strfind(subject, 'SS');
    subject = subject(ind:end);
    folderSession = fullfile(folderSave, subject, db(k).date, '001');
    if ~isfolder(folderSession)
        mkdir(folderSession)
    end
    
    % spike data
    alignDir = fullfile(folderBase, db(k).subject, db(k).date, 'alignments');
    sp = loadAllKsDir(db(k).subject, db(k).date);
    probe = db(k).OTprobe;
    rawF = fullfile(folderBase, db(k).subject, db(k).date, ...
        ['ephys_' sp(probe).name], 'sorting');
    if isfolder(rawF)
        [~, ~, spikeDepths] = ksDriftmap(rawF);
    end
    
    units = db(k).OTunits{1}(db(k).OTgood{1}==1);
    ind = ismember(sp(probe).clu, units);
    writeNPY(sp(probe).st(ind), fullfile(folderSession, ...
        'spike.times.npy'));
    writeNPY(sp(probe).spikeAmps(ind), fullfile(folderSession, ...
        'spike.amps.npy'));
    writeNPY(sp(probe).clu(ind), fullfile(folderSession, ...
        'spike.clusters.npy'));
    writeNPY(spikeDepths(ind), fullfile(folderSession, ...
        'spike.depths.npy'));
    
    writematrix(units, fullfile(folderSession, 'clusters.uuids.csv'));
    
    goodIntervals = [];
    goodClusters = [];
    for iCell = 1:length(units)
        gi = db(k).OTtimes{1}{db(k).OTunits{1}==units(iCell)};
        goodIntervals = [goodIntervals; gi];
        goodClusters = [goodClusters; ones(size(gi,1),1).*units(iCell)];
    end
    
    writeNPY(goodClusters, fullfile(folderSession, ...
        '_ss_validTimes.clusters.npy'));
    writeNPY(goodIntervals, fullfile(folderSession, ...
        '_ss_validTimes.intervals.npy'));
    
    writeNPY([sp(probe).xcoords sp(probe).ycoords], fullfile(folderSession, ...
        'channels.localCoordinates.npy'));
    
    writeNPY(sp(probe).sample_rate, fullfile(folderSession, ...
        'probe._ss_sampleRate.npy'));
    
    % running data
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
    
    writeNPY(runSpeed, fullfile(folderSession, ...
        '_ss_running.speed.npy'));
    writeNPY(runTime, fullfile(folderSession, ...
        '_ss_running.timestamps.npy'));
    
    % spike waveforms
    if ~strcmp(sp(probe).name, db(k).probeNames{1})
        bEphysToMaster = readNPY(fullfile(alignDir, ...
            sprintf('correct_ephys_%s_to_ephys_%s.npy', sp(probe).name, ...
            db(k).probeNames{1})));
    else % this one is master, so use a dummy conversion
        bEphysToMaster = [1; 0];
    end
    run = interp1(runTime, runSpeed, timeBins, 'pchip'); % run speed during darkness only
    run = [ones(1,length(win)) .* mean(run(1:min(length(run),length(win)))), ...
        run, ones(1,length(win)) .* mean(run(end-min(length(run),length(win))+1:end))];
    run = conv(run, win, 'same');
    run = run(length(win)+1 : end-length(win));
    isRunning = run > runThreshold;
    params.dataDir = fullfile(folderBase, db(k).subject, db(k).date, ...
        ['ephys_' sp(probe).name]);
    params.fileName = sprintf('%s_%s_%s_g0_t0.imec.ap_CAR.bin', ...
        db(k).subject, db(k).date, sp(probe).name);
    params.dataType = 'int16';
    params.nCh = 385;
    params.wfWin = round([-0.001 .002] .* sp(probe).sample_rate);
    params.nWf = 10000;
    waveforms = NaN(length(units), size(sp(probe).tempsUnW,3), ...
        sum(abs(params.wfWin))+1, 2); % [unit x channel x time x (running/not running)]
    for iCell = 1:length(units)
        st_dark = sp(probe).st(sp(probe).clu==units(iCell));
        st_dark(st_dark<binEdges(1) | st_dark>binEdges(end)) = [];
        goodInt = find(goodClusters == units(iCell));
        if ~isempty(goodInt)
            st = [];
            for g = 1:length(goodInt)
                st = [st; st_dark(st_dark >= goodIntervals(g,1) & ...
                    st_dark <= goodIntervals(g,2))];
            end
            st_dark = st;
        end
        [~,~,b] = histcounts(st_dark, binEdges);
        % for spikeTimes: first convert spike times (in master time) to time of
        % probe (if this probe was not the "master" probe, then multiply by
        % sampling rate to get spike times in samples (not seconds)
        params.spikeTimes = round((st_dark - bEphysToMaster(2)) ./ ...
            bEphysToMaster(1) .* sp(probe).sample_rate); % spikes of this unit during darkness
        clusters = isRunning(b)'+1; % whether animal was running or not, at each spike time
        params.spikeClusters = clusters;
        wf = getWaveForms(params); % [(running/not running) x channel x time]
        waveforms(iCell,:,:,:) = permute(wf.waveFormsMean, [2 3 1]);
    end
    writeNPY(waveforms, fullfile(folderSession, ...
        'clusters.waveforms.npy'));
    
    % stimulus data
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
        flickerAll = reshape([flickerWhite, flickerBlack]', [], 1);
        longestFlicker = max(parsFlicker.pars( ...
            strcmp(parsFlicker.parnames, 'nfr'),:)) / monitorRefreshRate;
        durs = diff(flickerAll);
        flickerStimStartInds = [1; find(durs > ...
            longestFlicker * 1.5) + 1; length(flickerAll)+1];
        flickerFreqs = NaN(length(flickerAll),1);
        flickerColor = NaN(length(flickerAll),1);
        flickerRep = NaN(length(flickerAll),1);
        for stim = 1:size(parsFlicker.seqnums,1)
            freq = 1/(flickerDurs(stim)*2);
            for trial = 1:size(parsFlicker.seqnums,2)
                count = parsFlicker.seqnums(stim,trial);
                ind = flickerStimStartInds(count) : ...
                    flickerStimStartInds(count+1)-1;
                flickerFreqs(ind) = freq;
                flickerColor(ind) = repmat([1;-1], length(ind)/2, 1);
                flickerRep(ind) = trial;
            end
        end
        
        writeNPY(flickerAll, fullfile(folderSession, ...
            '_ss_flicker.times.npy'));
        writeNPY(flickerFreqs, fullfile(folderSession, ...
            '_ss_flicker.frequencies.npy'));
        writeNPY(flickerColor, fullfile(folderSession, ...
            '_ss_flicker.color.npy'));
        writeNPY(flickerRep, fullfile(folderSession, ...
            '_ss_flicker.repetition.npy'));
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
        stimFrames = permute(stimFrames, [3 1 2]);
        framesPerImage = parsNoise.pars(6,1);
        frameTimes = (0 : size(stimFrames, 1)-1) * framesPerImage / myScreenInfo.FrameRate;
        noiseFrameTimes = (frameTimes + noiseOn)';
        noisePosition = parsNoise.pars(2:5)./10;
        
        writeNPY(noiseFrameTimes(:), fullfile(folderSession, ...
            '_ss_sparseNoise.times.npy'));
        writeNPY(repmat((1:length(frameTimes))', length(noiseOn), 1), ...
            fullfile(folderSession, '_ss_sparseNoise._ss_sparseNoiseID.npy'));
        
        writeNPY(noisePosition', fullfile(folderSession, ...
            '_ss_sparseNoiseArea.edges.npy'));
        
        writeNPY(stimFrames, fullfile(folderSession, ...
            '_ss_sparseNoiseID.map.npy'));
        
        % Fitted RFs
        ind = true(length(RFs(k).explainedVariances),1);
        if k == 5
            ind(6) = false; % exclude unit 66 because there was no running during darkness
        end
        writeNPY(permute(RFs(k).receptiveFields(:,:,:,:,ind), [5 1 2 3 4]), ...
            fullfile(folderSession, '_ss_rf.maps.npy'));
        writeNPY(RFs(k).explainedVariances(ind), ...
            fullfile(folderSession, '_ss_rf.explVars.npy'));
        writeNPY(RFs(k).explainedVariances_runOnly(ind), ...
            fullfile(folderSession, '_ss_rf.explVarsRunning.npy'));
        writeNPY(RFs(k).explainedVariances_stimOnly(ind), ...
            fullfile(folderSession, '_ss_rf.explVarsStim.npy'));
        writeNPY(RFs(k).lambdasRun(ind)', ...
            fullfile(folderSession, '_ss_rf.lambdasRunning.npy'));
        writeNPY(RFs(k).lambdasStim(ind)', ...
            fullfile(folderSession, '_ss_rf.lambdasStim.npy'));
        writeNPY(RFs(k).pVal_RFonly(ind), ...
            fullfile(folderSession, '_ss_rf.pValues.npy'));
        writeNPY(RFs(k).RFTimes, ...
            fullfile(folderSession, '_ss_rfDescr.timestamps.npy'));
        writeNPY(RFs(k).stimPosition', ...
            fullfile(folderSession, '_ss_rfDescr.edges.npy'));
        writeNPY(RFs(k).runningKernels(:,ind), ...
            fullfile(folderSession, '_ss_rfRunningKernels.rate.npy'));
        writeNPY(RFs(k).runWindow', ...
            fullfile(folderSession, '_ss_rfRunningKernels.timestamps.npy'));
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
        
        writeNPY([stimOn stimOff], fullfile(folderSession, ...
            '_ss_grating.intervals.npy'));
        writeNPY(stimSeq, fullfile(folderSession, ...
            '_ss_grating._ss_gratingID.npy'));
        
        writeNPY(directions', fullfile(folderSession, ...
            '_ss_gratingID.directions.npy'));
    end
    
    writeNPY(db(k).darkTime, ...
        fullfile(folderSession, '_ss_darkness.intervals.npy'));
    
    ind = true(length(corrs(k).units),1);
    if k == 5
        ind(6) = false; % exclude unit 66 because there was no running during darkness
    end
    writeNPY(cat(1, corrs(k).units(ind).crosscorr)', ...
        fullfile(folderSession, '_ss_crossCorrs.values.npy'));
    writeNPY(permute(cat(3, corrs(k).units(ind).nullCrossCorr), [2 1 3]), ...
        fullfile(folderSession, '_ss_crossCorrs.nullValues.npy'));
    writeNPY(corrs(k).lags', ...
        fullfile(folderSession, '_ss_crossCorrs.timestamps.npy'));
end