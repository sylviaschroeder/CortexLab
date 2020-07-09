db_ephys_driftingGratings

%% Folders
folderBase = '\\ZUBJECTS.cortexlab.net\Subjects';
folderTuning = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\Electrophys_SC\nonVisualEffects\pupil';
folderEye = '\\zserver.cortexlab.net\Data\EyeCamera';
folderSave = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\DataToPublish\sc neurons ephys';

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
d = load(fullfile(folderTuning, 'tuning_prefDirSigmaDIFixed_lineFit.mat'));
tuning = d.tuning;

for k = 1:length(db)
    fprintf('Dataset %d of %d\n', k, length(db))
    data = [];
    subject = db(k).subject;
    folderSession = fullfile(folderSave, subject, db(k).date, '001');
    if ~isfolder(folderSession)
        mkdir(folderSession)
    end
%     
%     alignDir = fullfile(folderBase, db(k).subject, db(k).date, 'alignments');
%     sp = loadAllKsDir(db(k).subject, db(k).date);
%     if ~isempty(db(k).tag)
%         probe = find(strcmp({sp.name}, db(k).tag));
%         name = ['ephys_' sp(probe).name];
%     else
%         probe = 1;
%         name = 'ephys';
%     end
%     rawF = fullfile(folderBase, db(k).subject, db(k).date, ...
%         name, 'sorting');
%     if isfolder(rawF)
%         [~, ~, spikeDepths] = ksDriftmap(rawF);
%     end
%     templates = findTempForEachClu(sp(probe).clu, sp(probe).spikeTemplates);
%     
%     units = sp(probe).cids(sp(probe).cgs == 2); % only single units
%     depths = NaN(length(units),1);
%     for iCell = 1:length(units)
%         depths(iCell) = nanmean(spikeDepths(sp(probe).clu == units(iCell)));
%     end
%     included = depths > db(k).yPosRange(1) & depths < db(k).yPosRange(2);
%     units = units(included);
%     ind = ismember(sp(probe).clu, units);
%     waveforms = NaN(length(units), size(sp(probe).tempsUnW,2), ...
%         size(sp(probe).tempsUnW,3));
%     for iCell = 1:length(units)
%         waveforms(iCell,:,:) = squeeze(sp(probe). ...
%             tempsUnW(templates(units(iCell)+1)+1,:,:));
%     end
%     
%     writeNPY(sp(probe).st(ind), fullfile(folderSession, ...
%         'spike.times.npy'));
%     writeNPY(sp(probe).spikeAmps(ind), fullfile(folderSession, ...
%         'spike.amps.npy'));
%     writeNPY(sp(probe).clu(ind), fullfile(folderSession, ...
%         'spike.clusters.npy'));
%     writeNPY(spikeDepths(ind), fullfile(folderSession, ...
%         'spike.depths.npy'));
%     
%     writematrix(units, fullfile(folderSession, 'clusters.uuids.csv'));
%     writeNPY(waveforms, fullfile(folderSession, ...
%         'clusters.waveforms.npy'));
%     
%     writeNPY([sp(probe).xcoords sp(probe).ycoords], fullfile(folderSession, ...
%         'channels.localCoordinates.npy'));
%     
%     writeNPY(sp(probe).sample_rate, fullfile(folderSession, ...
%         'probe._ss_sampleRate.npy'));
%     
%     [expNums, ~, ~, ~, ~, tl, hasTimeline] = ...
%         dat.whichExpNums(db(k).subject, db(k).date);
%     TLexp = expNums(hasTimeline);
%     TLexp = TLexp(end);
%     tl = tl{end};
%     bTLtoMaster = readNPY(fullfile(alignDir, ...
%         sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, db(k).masterTimebase)));
%     tlTime = applyCorrection(tl.rawDAQTimestamps, bTLtoMaster);
%     
%     rotary = double(tl.rawDAQData(:, strcmp({tl.hw.inputs.name}, ...
%         'rotaryEncoder')));
%     runSpeed = nonVis.getRunningSpeed_wheel(rotary, tlTime, 0);
%     runTime = runSpeed.t;
%     runSpeed = runSpeed.total;
%     cmPerUnit = 2*pi * 8.75 / (4 * 1024);
%     runSpeed = runSpeed * cmPerUnit;
%     
%     writeNPY(runSpeed, fullfile(folderSession, ...
%         '_ss_running.speed.npy'));
%     writeNPY(runTime, fullfile(folderSession, ...
%         '_ss_running.timestamps.npy'));
%     
%     if isfile(fullfile(folderEye, db(k).subject, db(k).date, ...
%             sprintf('%02d_eye_processed.mat', db(k).exp)))
%         data = load(fullfile(folderEye, db(k).subject, db(k).date, ...
%             sprintf('%02d_eye_processed.mat', db(k).exp)), 'results');
%         pupilSize = nonVis.getPupilDiam(data.results);
%         pupilXY = [data.results.x, data.results.y];
%         pupilXY(isnan(data.results.area) | data.results.blink | ...
%             ~data.results.goodFit,:) = NaN;
%         data = load(fullfile(folderEye, db(k).subject, db(k).date, ...
%             sprintf('%02d_eyeTime.mat', db(k).exp)));
%         time_pupilSize = data.eyeTime';
%         
%         writeNPY(pupilSize, fullfile(folderSession, 'eye.diameter.npy'));
%         writeNPY(pupilXY, fullfile(folderSession, 'eye.xyPos.npy'));
%         writeNPY(time_pupilSize, fullfile(folderSession, ...
%             'eye.timestamps.npy'));
%     end
%     
%     d = load(fullfile(folderBase, db(k).subject, db(k).date, ...
%         num2str(db(k).exp), sprintf('%s_%d_%s_parameters.mat', ...
%         db(k).date, db(k).exp, db(k).subject)));
%     parsGratings = d.parameters.Protocol;
%     stimOnTL = readNPY(fullfile(alignDir, ...
%         sprintf('mpep_%d_onsets_in_timeline_%d.npy', db(k).exp, TLexp)));
%     stimOffTL = readNPY(fullfile(alignDir, ...
%         sprintf('mpep_%d_offsets_in_timeline_%d.npy', db(k).exp, TLexp)));
%     stimOn = applyCorrection(stimOnTL, bTLtoMaster);
%     stimOff = applyCorrection(stimOffTL, bTLtoMaster);
%     stimIDs = repmat((1:parsGratings.npfilestimuli)',parsGratings.nrepeats,1);
%     [~,order] = sort(parsGratings.seqnums(:));
%     stimSeq = stimIDs(order);
%     stimDur = round(parsGratings.pars(strcmp(parsGratings.parnames, ...
%         tags.stimDur{1}),:)' ./ 10 .* 60) ./ 60; % 60 Hz monitor frame rate
%     stimOn_rel = round(parsGratings.pars(strcmp(parsGratings.parnames, ...
%         tags.stimOnset{1}),:)' ./ 1000 .* 60) ./ 60;
%     stimOff_rel = round(parsGratings.pars(strcmp(parsGratings.parnames, ...
%         tags.stimOffset{1}),:)' ./ 1000 .* 60) ./ 60 - stimDur;
%     laserOn_rel = round(parsGratings.pars(strcmp(parsGratings.parnames, ...
%         tags.laserOn{1}),:)' ./ 1000 .* 60) ./ 60; 
%     laserOff_rel = round(parsGratings.pars(strcmp(parsGratings.parnames, ...
%         tags.laserOff{1}),:)' ./ 1000 .* 60) ./ 60 - stimDur;
%     laserOn_rel = laserOn_rel - stimOn_rel;
%     laserOff_rel = laserOff_rel - stimOff_rel;
%     stimOn_rel = stimOn_rel(stimSeq);
%     stimOff_rel = stimOff_rel(stimSeq);
%     stimOn = stimOn + stimOn_rel;
%     stimOff = stimOff + stimOff_rel;
%     stimDur = mean(stimOff - stimOn);
%     tag = tags.contrast{ismember(tags.contrast, parsGratings.parnames)};
%     blank = parsGratings.pars(strcmp(parsGratings.parnames,tag),:) == 0;
%     tag = tags.ori{ismember(tags.ori, parsGratings.parnames)};
%     directions = parsGratings.pars(strcmp(parsGratings.parnames,tag),:);
%     directions(blank) = NaN;
%     laserOn = parsGratings.pars(strcmp(parsGratings.parnames,tags.amplitude),:) > 0;
%     
%     writeNPY([stimOn stimOff], fullfile(folderSession, ...
%         '_ss_grating.intervals.npy'));
%     writeNPY(stimSeq, fullfile(folderSession, ...
%         '_ss_grating._ss_gratingID.npy'));
%     
%     writeNPY(directions', fullfile(folderSession, ...
%         '_ss_gratingID.directions.npy'));
%     writeNPY(laserOn', fullfile(folderSession, ...
%         '_ss_gratingID.laserOn.npy'));
%     writeNPY(laserOn_rel, fullfile(folderSession, ...
%         '_ss_gratingID.laserOnTime.npy'));
%     writeNPY(laserOff_rel, fullfile(folderSession, ...
%         '_ss_gratingID.laserOffTime.npy'));
    
    % Fitted tuning curves
    kTun = find(strcmp(db(k).subject, {tuning.subject}) & ...
        strcmp(db(k).date, {tuning.date}) & db(k).exp==[tuning.exp]);
    numCells = length(tuning(kTun).cellIDs);
    parsSLoff = NaN(numCells, 5);
    parsLLoff = NaN(numCells, 5);
    parsSLon = NaN(numCells, 5);
    parsLLon = NaN(numCells, 5);
    curvSLoff = NaN(numCells, 360);
    curvLLoff = NaN(numCells, 360);
    curvSLon = NaN(numCells, 360);
    curvLLon = NaN(numCells, 360);
    for iCell = 1:numCells
        p = tuning(kTun).cond(1).cell(iCell).parameters;
        parsSLoff(iCell,1:length(p)) = p;
        p = tuning(kTun).cond(2).cell(iCell).parameters;
        parsLLoff(iCell,1:length(p)) = p;
        p = tuning(kTun).cond(3).cell(iCell).parameters;
        parsSLon(iCell,1:length(p)) = p;
        p = tuning(kTun).cond(4).cell(iCell).parameters;
        parsLLon(iCell,1:length(p)) = p;
        curvSLoff(iCell,:) = tuning(kTun).cond(1).cell(iCell).curve;
        curvLLoff(iCell,:) = tuning(kTun).cond(2).cell(iCell).curve;
        curvSLon(iCell,:) = tuning(kTun).cond(3).cell(iCell).curve;
        curvLLon(iCell,:) = tuning(kTun).cond(4).cell(iCell).curve;
    end
    writeNPY(tuning(kTun).isSuppressed, fullfile(folderSession, ...
        'clusters.isSuppressed.npy'));
    writeNPY(parsSLoff, fullfile(folderSession, ...
        'clusters.tuningParametersSmallLaserOff.npy'));
    writeNPY(parsLLoff, fullfile(folderSession, ...
        'clusters.tuningParametersLargeLaserOff.npy'));
    writeNPY(parsSLon, fullfile(folderSession, ...
        'clusters.tuningParametersSmallLaserOn.npy'));
    writeNPY(parsLLon, fullfile(folderSession, ...
        'clusters.tuningParametersLargeLaserOn.npy'));
    writeNPY(curvSLoff, fullfile(folderSession, ...
        'clusters.tuningCurvesSmallLaserOff.npy'));
    writeNPY(curvLLoff, fullfile(folderSession, ...
        'clusters.tuningCurvesLargeLaserOff.npy'));
    writeNPY(curvSLon, fullfile(folderSession, ...
        'clusters.tuningCurvesSmallLaserOn.npy'));
    writeNPY(curvLLon, fullfile(folderSession, ...
        'clusters.tuningCurvesLargeLaserOn.npy'));
    writeNPY(tuning(kTun).nonVisual==2, fullfile(folderSession, ...
        '_ss_gratingID.largePupil.npy'));
end