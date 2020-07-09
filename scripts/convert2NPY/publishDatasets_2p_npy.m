% units = 'boutons';
units = 'neurons';

timeGap = 600; %in s, gap between experiments

%% Load database
if strcmp(units, 'boutons')
    db_boutons_driftingGratings_blanks
else
    db_driftingGratings_blank
end

%% Folders
folderBase = 'C:\STORAGE\OneDrive - University College London\Lab';
folderROIData = fullfile(folderBase, 'DATA\InfoStructs');
if strcmp(units, 'boutons')
    folderTuning = fullfile(folderBase, 'RESULTS\boutons\nonVisualEffects');
    folderRF = fullfile(folderBase, 'RESULTS\receptiveFields\boutons');
    d = load(fullfile(folderROIData, 'corrections_boutons.mat'));
    corrections = d.corrections;
    doCorrect = d.doCorrect;
    folderSave = fullfile(folderBase, 'DATA\DataToPublish\boutons');
else
    folderTuning = fullfile(folderBase, 'RESULTS\nonVisualEffects\modelGratingResp');
    folderRF = fullfile(folderBase, 'RESULTS\receptiveFields\SC neurons');
    corrections = [];
    folderSave = fullfile(folderBase, 'DATA\DataToPublish\sc neurons 2p');
end

%% Collect data
experiments = {'expGratings', 'expGrayScreen', 'expDark', 'expNoise'};
expNames = {'gratings','grayScreen','darkness','sparseNoise'};

d = load(fullfile(folderTuning, 'pupil', 'corrsDuringGratingsAndGrayScreen_sigma1.00.mat'));
corrsPupil = d.corrs;
d = load(fullfile(folderTuning, 'running', 'corrsDuringGratingsAndGrayScreen_sigma1.00.mat'));
corrsRunning = d.corrs;
d = load(fullfile(folderTuning, 'kernelFit', 'results.mat'));
results = d.results;
d = load(fullfile(folderTuning, 'pupil', 'tuning_prefDirSigmaDIFixed_lineFit.mat'));
tuning = d.tuning;
d = load(fullfile(folderTuning, 'pupil', 'nullTuning_prefDirSigmaDIFixed.mat'));
null = d.null;
d = load(fullfile(folderRF, 'receptiveFields.mat'));
RFs = d.RFs;

for k = 1:length(db)
    fprintf('Dataset %d of %d\n', k, length(db))
    subject = db(k).subject;
    ind = strfind(subject, 'SS');
    subject = subject(ind:end);
    folderSession = fullfile(folderSave, subject, db(k).date, '001');
    if ~isfolder(folderSession)
        mkdir(folderSession)
    end
    
    valid = [];
    xyzPos = [];
    planePerUnit = [];
    cellIDs = [];
    isGad = [];
    infoDone = false;
    F_all = [];
    time = [];
    
    runningSpeed = [];
    time_runningSpeed = [];
    pupilSize = [];
    pupilXY = [];
    time_pupilSize = [];
    
    rhosPupilGratings = [];
    nullRhosPupilGratings = [];
    rhosRunningGratings = [];
    nullRhosRunningGratings = [];
    rhosPupilGray = [];
    nullRhosPupilGray = [];
    rhosRunningGray = [];
    nullRhosRunningGray = [];
    rhosRunningDarkness = [];
    nullRhosRunningDarkness = [];
    
    amps_all = [];
    kernels_all = [];
    preds_all = [];
    
    largePupil = [];
    isSuppressed = [];
    tunExplVar = [];
    paramsLarge = [];
    paramsSmall = [];
    curvesLarge = [];
    curvesSmall = [];
    pupilCondDone = false;
    nullParamsLarge = [];
    nullParamsSmall = [];
    
    receptiveFields = [];
    runningKernels = [];
    rfExplVar = [];
    rfEvRun = [];
    rfEvStim = [];
    rfLambdaRun = [];
    rfLambdaStim = [];
    rfPval = [];
    rfTime = [];
    rfEdges = [];
    runKernelTime = [];
    t0 = 0;
    for exp = 1:length(experiments)
        if ~isfield(db, experiments{exp}) || isempty(db(k).(experiments{exp}))
            continue
        end
        
        folder = fullfile(folderROIData, db(k).subject, ...
            db(k).date, num2str(db(k).(experiments{exp})));
        file = [sprintf('%s_%d_%s', db(k).date, db(k).(experiments{exp}), ...
            db(k).subject) '_2P_plane%03d_ROI.mat'];
        valid_exp = [];
        F_exp = [];
        for iPlane = 1:length(db(k).planes)
            % load data
            d = load(fullfile(folder, sprintf(file, db(k).planes(iPlane))));
            meta = d.meta;
            meta.folderTL = strrep(meta.folderTL, 'ioo.ucl.ac.uk', ...
                'cortexlab.net');
            numCells = size(meta.F_final,2);
            
            if ~infoDone
                xyzPos = [xyzPos; cat(1, meta.ROI.CellXYZMicrons{:})];
                planePerUnit = [planePerUnit; ones(numCells,1) .* iPlane];
                cellIDs = [cellIDs; (1:numCells)'];
                if strcmp(units, 'boutons')
                    g = NaN(numCells,1);
                else
                    g = meta.ROI.isGad;
                end
                isGad = [isGad; g];
            end
            
            if iPlane == 1
                ballData = nonVis.getRunningSpeed(meta);
                if ~isempty(ballData)
                    runningSpeed = [runningSpeed; ...
                        ballData.total / median(diff(ballData.t)) / 53];
                    time_runningSpeed = [time_runningSpeed; ballData.t + t0];
                end
                [pupilData, tp] = nonVis.loadPupilData(meta);
                if ~isempty(pupilData) && ~strcmp(expNames{exp}, 'darkness')
                    pupilSize = [pupilSize; nonVis.getPupilDiam(pupilData)];
                    time_pupilSize = [time_pupilSize; tp(1:length(pupilData.x)) + t0];
                    xy = [pupilData.x, pupilData.y];
                    xy(isnan(pupilData.area) | pupilData.blink | ~pupilData.goodFit,:) = NaN;
                    pupilXY = [pupilXY; xy];
                end
                frameTimes = ppbox.getFrameTimes(meta)';
                frameDur = median(diff(frameTimes));
                time = [time; frameTimes + t0];
                timeShifts = (db(k).planes - db(k).planes(1)) .* ...
                    (frameDur / meta.nPlanes);
                
                switch expNames{exp}
                    case 'gratings'
                        [stimTimes, stimSeq, stimMatrix] = ...
                            ssLocal.getStimulusResponseInfo(meta);
                        [directions, blank] = gratings.getOrientations(stimSeq);
                        dirs = directions(:,1);
                        dirs(blank) = NaN;
                        
                        writeNPY([stimTimes.onset, stimTimes.offset] + t0, ...
                            fullfile(folderSession, '_ss_grating.intervals.npy'));
                        writeNPY(stimSeq.seq, ...
                            fullfile(folderSession, '_ss_grating._ss_gratingID.npy'));
                        writeNPY(dirs, ...
                            fullfile(folderSession, '_ss_gratingID.directions.npy'));
                        writeNPY(results(k).plane(1).kernelTime', ...
                            fullfile(folderSession, '_ss_gratingKernels.timestamps.npy'));
                        writeNPY(frameTimes + t0, ...
                            fullfile(folderSession, '_ss_gratingPredictions.timestamps.npy'))
                    case 'sparseNoise'
                        stimTimes = ppbox.getStimTimes(meta);
                        [stimFrames, stimPosition] = whiteNoise.getStimulusFrames(meta);
                        stimFrames = permute(stimFrames, [3 1 2]);
                        stimFrameDur = mean(stimTimes.offset - stimTimes.onset) / size(stimFrames,1);
                        stimFrameTimes = ((1:size(stimFrames,1))-1) .* stimFrameDur;
                        allStimFrameTimes = reshape((stimTimes.onset + stimFrameTimes)', [], 1);
                        
                        writeNPY(allStimFrameTimes + t0, ...
                            fullfile(folderSession, '_ss_sparseNoise.times.npy'));
                        writeNPY(repmat((1:size(stimFrames,1))', length(stimTimes.onset), 1), ...
                            fullfile(folderSession, '_ss_sparseNoise._ss_sparseNoiseID.npy'));
                        writeNPY(stimPosition, ...
                            fullfile(folderSession, '_ss_sparseNoiseArea.edges.npy'));
                        writeNPY(stimFrames, ...
                            fullfile(folderSession, '_ss_sparseNoiseID.map.npy'));
                end
            end
            
            F = meta.F_final;
            if ~isempty(corrections)
                a = corrections(k).plane(iPlane).a{db(k).(experiments{exp})};
                b = corrections(k).plane(iPlane).b{db(k).(experiments{exp})};
                F = doCorrect(a,b,F);
            end
            if size(F_exp,1) > size(F,1)
                F = padarray(F, size(F_exp,1)-size(F,1), NaN, 'post');
            elseif size(F_exp,1)>0 && size(F_exp,1)<size(F,1)
                F(size(F_exp,1)+1:end,:) = [];
            end
            F_exp = [F_exp, F];
            
            if strcmp(expNames{exp}, 'gratings')
                % correlations
                rhosPupilGratings = [rhosPupilGratings; ...
                    corrsPupil(k).plane(iPlane).gratings.rhos];
                nullRhosPupilGratings = [nullRhosPupilGratings; ...
                    corrsPupil(k).plane(iPlane).gratings.nullRhos];
                rhosRunningGratings = [rhosRunningGratings; ...
                    corrsRunning(k).plane(iPlane).gratings.rhos];
                nullRhosRunningGratings = [nullRhosRunningGratings; ...
                    corrsRunning(k).plane(iPlane).gratings.nullRhos];
                
                % kernel fit for gratings
                kernels = NaN(length(results(k).plane(iPlane).kernelTime), numCells);
                preds = NaN(length(frameTimes), numCells); 
                for iCell = 1:length(results(k).plane(iPlane).cellIDs)
                    a = results(k).plane(iPlane).kernelFit(iCell).alphaEachTrial;
                    if isempty(a)
                        continue
                    end
                    kernels(:,results(k).plane(iPlane).cellIDs(iCell)) = ...
                        results(k).plane(iPlane).kernelFit(iCell).kernel;
                    p = results(k).plane(iPlane).kernelFit(iCell).prediction;
                    preds(1:length(p),results(k).plane(iPlane).cellIDs(iCell)) = p;
                end
                if ~isempty(corrections)
                    a = corrections(k).plane(iPlane).a{db(k).(experiments{exp})}';
                    b = corrections(k).plane(iPlane).b{db(k).(experiments{exp})}';
                    preds = doCorrect(permute(a,[2 1]),permute(b,[2 1]),preds);
                end
                if size(preds_all,1) > size(preds,1)
                    preds = padarray(preds, size(preds_all,1)-size(preds,1), ...
                        NaN, 'post');
                elseif size(preds_all,1)>0 && size(preds_all,1)<size(preds,1)
                    preds(size(preds_all,1)+1:end,:) = [];
                end
                kernels_all = [kernels_all, kernels];
                preds_all = [preds_all, preds];
                
                % tuning curve fits
                kTun = find(strcmp(db(k).subject, {tuning.subject}) & ...
                    strcmp(db(k).date, {tuning.date}));
                amps = NaN(size(directions,1), ...
                    size(results(k).plane(iPlane).responses,3), numCells);
                sup = NaN(numCells,1);
                parL = NaN(numCells,5);
                parS = NaN(numCells,5);
                nParsL = NaN(numCells,5,200);
                nParsS = NaN(numCells,5,200);
                cL = NaN(numCells,360);
                cS = NaN(numCells,360);
                ev = NaN(numCells,1);
                for iCell = 1:length(tuning(kTun).plane(iPlane).cellIDs)
                    if isempty(tuning(kTun).plane(iPlane).cond(1).cell(iCell).responses)
                        continue
                    end
                    id = tuning(kTun).plane(iPlane).cellIDs(iCell);
                    a1 = tuning(kTun).plane(iPlane).cond(1).cell(iCell).responses;
                    a2 = tuning(kTun).plane(iPlane).cond(2).cell(iCell).responses;
                    amps(:,:,id) = nanmean(cat(3, a1, a2), 3);
                    sup(id) = tuning(kTun).plane(iPlane).isSuppressed(iCell);
                    p = tuning(kTun).plane(iPlane).cond(1).cell(iCell).parameters;
                    parS(id,1:length(p)) = p;
                    p = tuning(kTun).plane(iPlane).cond(2).cell(iCell).parameters;
                    parL(id,1:length(p)) = p;
                    p = null(kTun).plane(iPlane).cond(1).cell(iCell).parameters;
                    nParsS(id,1:size(p,2),:) = p';
                    p = null(kTun).plane(iPlane).cond(2).cell(iCell).parameters;
                    nParsL(id,1:size(p,2),:) = p';
                    cS(id,:) = tuning(kTun).plane(iPlane).cond(1).cell(iCell).curve;
                    cL(id,:) = tuning(kTun).plane(iPlane).cond(2).cell(iCell).curve;
                    ev(id) = tuning(kTun).plane(iPlane).crossValExplVar(iCell);
                    
                    if ~pupilCondDone
                        rS = tuning(kTun).plane(iPlane).cond(1).cell(iCell).responses;
                        largePupil = isnan(rS);
                        pupilCondDone = true;
                    end
                end
                if ~isempty(corrections)
                    a = corrections(k).plane(iPlane).a{db(k).(experiments{exp})}';
                    b = corrections(k).plane(iPlane).b{db(k).(experiments{exp})}';
                    amps = doCorrect(permute(a,[2 3 1]),permute(b,[2 3 1]),amps);
                    cL = doCorrect(a,b,cL);
                    cS = doCorrect(a,b,cS);
                    parL(:,[2 4]) = doCorrect(a,b,parL(:,[2 4]));
                    parS(:,[2 4]) = doCorrect(a,b,parS(:,[2 4]));
                    nParsL(:,[2 4],:) = doCorrect(a,b,nParsL(:,[2 4],:));
                    nParsS(:,[2 4],:) = doCorrect(a,b,nParsS(:,[2 4],:));
                end
                amps_all = cat(3, amps_all, amps);
                isSuppressed = [isSuppressed; sup];
                paramsLarge = [paramsLarge; parL];
                paramsSmall = [paramsSmall; parS];
                nullParamsLarge = [nullParamsLarge; nParsL];
                nullParamsSmall = [nullParamsSmall; nParsS];
                curvesLarge = [curvesLarge; cL];
                curvesSmall = [curvesSmall; cS];
                tunExplVar = [tunExplVar; ev];
                
            elseif strcmp(expNames{exp}, 'grayScreen')
                % correlations
                rhosPupilGray = [rhosPupilGray; ...
                    corrsPupil(k).plane(iPlane).grayScreen.rhos];
                nullRhosPupilGray = [nullRhosPupilGray; ...
                    corrsPupil(k).plane(iPlane).grayScreen.nullRhos];
                rhosRunningGray = [rhosRunningGray; ...
                    corrsRunning(k).plane(iPlane).grayScreen.rhos];
                nullRhosRunningGray = [nullRhosRunningGray; ...
                    corrsRunning(k).plane(iPlane).grayScreen.nullRhos];
                
            elseif strcmp(expNames{exp}, 'sparseNoise')
                kRF = find(strcmp(db(k).subject, {RFs.subject}) & ...
                    strcmp(db(k).date, {RFs.date}));                
                receptiveFields = [receptiveFields; ...
                    permute(RFs(kRF).plane(iPlane).receptiveFields, [5 1 2 3 4])];
                runningKernels = [runningKernels, ...
                    RFs(kRF).plane(iPlane).runningKernels];
                rfExplVar = [rfExplVar; RFs(kRF).plane(iPlane).explainedVariances];
                rfEvRun = [rfEvRun; RFs(kRF).plane(iPlane).explainedVariances_runOnly];
                rfEvStim = [rfEvStim; RFs(kRF).plane(iPlane).explainedVariances_stimOnly];
                rfLambdaRun = [rfLambdaRun; RFs(kRF).plane(iPlane).lambdasRun'];
                rfLambdaStim = [rfLambdaStim; RFs(kRF).plane(iPlane).lambdasStim'];
                rfPval = [rfPval; RFs(kRF).plane(iPlane).pVal_RFonly];
                if iPlane == 1
                    rfTime = RFs(kRF).RFTimes;
                    rfEdges = RFs(kRF).stimPosition;
                    runKernelTime = RFs(kRF).runWindow';
                end
            elseif strcmp(expNames{exp}, 'darkness') && strcmp(units, 'boutons')
                % correlations
                rhosRunningDarkness = [rhosRunningDarkness; ...
                    corrsRunning(k).plane(iPlane).dark.rhos];
                nullRhosRunningDarkness = [nullRhosRunningDarkness; ...
                    corrsRunning(k).plane(iPlane).dark.nullRhos];
            end
            
            % valid ROIs: unique and no "switch-on" responses
            if isfield(meta.ROI, 'isDuplicate')
                valid_exp = [valid_exp; meta.ROI.isDuplicate == 0 & ...
                    meta.ROI.isSwitchOn == 0 & ~all(isnan(meta.F_final),1)'];
            else
                valid_exp = [valid_exp; ~all(isnan(meta.F_final),1)'];
            end
        end
        infoDone = true;
        valid = [valid, valid_exp];
        F_all = [F_all; F_exp];
        
        writeNPY(frameTimes([1 end])' + t0, fullfile(folderSession, ...
            sprintf('_ss_recordings.%s_intervals.npy', expNames{exp})));
        
        t0 = time(end) + timeGap;
    end
    valid = all(valid,2);
    
    writeNPY(planePerUnit(valid), ...
        fullfile(folderSession, '_ss_2pRois._ss_2pPlanes.npy'));
    writeNPY(cellIDs(valid), ...
        fullfile(folderSession, '_ss_2pRois.ids.npy'));
    writeNPY(xyzPos(valid,:), ...
        fullfile(folderSession, '_ss_2pRois.xyz.npy'));
    writeNPY(isGad(valid), ...
        fullfile(folderSession, '_ss_2pRois.isGad.npy'));
    writeNPY(timeShifts, fullfile(folderSession, '_ss_2pPlanes.delay.npy'));
    writeNPY(time, fullfile(folderSession, '_ss_2pCalcium.timestamps.npy'));
    writeNPY(F_all(:,valid), ...
        fullfile(folderSession, '_ss_2pCalcium.dff.npy'));
    
    if ~isempty(db(k).expGratings)
        writeNPY(rhosPupilGratings(valid), ...
            fullfile(folderSession, '_ss_corrsPupil.rhosGratings.npy'));
        writeNPY(nullRhosPupilGratings(valid,:), ...
            fullfile(folderSession, '_ss_corrsPupil.nullRhosGratings.npy'));
        writeNPY(rhosRunningGratings(valid), ...
            fullfile(folderSession, '_ss_corrsRunning.rhosGratings.npy'));
        writeNPY(nullRhosRunningGratings(valid,:), ...
            fullfile(folderSession, '_ss_corrsRunning.nullRhosGratings.npy'));
        writeNPY(largePupil, ...
            fullfile(folderSession, '_ss_gratingTrials.largePupil.npy'));
        writeNPY(amps_all(:,:,valid), ...
            fullfile(folderSession, '_ss_gratingTrials.amplitudes.npy'));
        writeNPY(kernels_all(:,valid), ...
            fullfile(folderSession, '_ss_gratingKernels.dff.npy'));
        writeNPY(preds_all(:,valid), ...
            fullfile(folderSession, '_ss_gratingPredictions.dff.npy'));
        writeNPY(isSuppressed(valid), ...
            fullfile(folderSession, '_ss_tuning.isSuppressed.npy'));
        writeNPY(paramsLarge(valid,:), ...
            fullfile(folderSession, '_ss_tuning.parametersLarge.npy'));
        writeNPY(paramsSmall(valid,:), ...
            fullfile(folderSession, '_ss_tuning.parametersSmall.npy'));
        writeNPY(nullParamsLarge(valid,:,:), ...
            fullfile(folderSession, '_ss_tuning.nullParametersLarge.npy'));
        writeNPY(nullParamsSmall(valid,:,:), ...
            fullfile(folderSession, '_ss_tuning.nullParametersSmall.npy'));
        writeNPY(curvesLarge(valid,:), ...
            fullfile(folderSession, '_ss_tuning.curvesLarge.npy'));
        writeNPY(curvesSmall(valid,:), ...
            fullfile(folderSession, '_ss_tuning.curvesSmall.npy'));
        writeNPY(tunExplVar(valid), ...
            fullfile(folderSession, '_ss_tuning.explVars.npy'));
    end
    if ~isempty(db(k).expGrayScreen)
        writeNPY(rhosPupilGray(valid), ...
            fullfile(folderSession, '_ss_corrsPupil.rhosGrayScreen.npy'));
        writeNPY(nullRhosPupilGray(valid,:), ...
            fullfile(folderSession, '_ss_corrsPupil.nullRhosGrayScreen.npy'));
        writeNPY(rhosRunningGray(valid), ...
            fullfile(folderSession, '_ss_corrsRunning.rhosGrayScreen.npy'));
        writeNPY(nullRhosRunningGray(valid,:), ...
            fullfile(folderSession, '_ss_corrsRunning.nullRhosGrayScreen.npy'));
    end
    if isfield(db, 'expDark') && ~isempty(db(k).expDark) && strcmp(units, 'boutons')
        writeNPY(rhosRunningDarkness(valid), ...
            fullfile(folderSession, '_ss_corrsRunning.rhosDark.npy'));
        writeNPY(nullRhosRunningDarkness(valid,:), ...
            fullfile(folderSession, '_ss_corrsRunning.nullRhosDark.npy'));
    end
    if ~isempty(db(k).expNoise)
        writeNPY(receptiveFields(valid,:,:,:,:), ...
            fullfile(folderSession, '_ss_rf.maps.npy'));
        writeNPY(rfExplVar(valid), ...
            fullfile(folderSession, '_ss_rf.explVars.npy'));
        writeNPY(rfEvRun(valid), ...
            fullfile(folderSession, '_ss_rf.explVarsRunning.npy'));
        writeNPY(rfEvStim(valid), ...
            fullfile(folderSession, '_ss_rf.explVarsStim.npy'));
        writeNPY(rfLambdaRun(valid), ...
            fullfile(folderSession, '_ss_rf.lambdasRunning.npy'));
        writeNPY(rfLambdaStim(valid), ...
            fullfile(folderSession, '_ss_rf.lambdasStim.npy'));
        writeNPY(rfPval(valid), ...
            fullfile(folderSession, '_ss_rf.pValues.npy'));
        writeNPY(rfTime, ...
            fullfile(folderSession, '_ss_rfDescr.timestamps.npy'));
        writeNPY(rfEdges, ...
            fullfile(folderSession, '_ss_rfDescr.edges.npy'));
        writeNPY(runningKernels(:,valid), ...
            fullfile(folderSession, '_ss_rfRunningKernels.dff.npy'));
        writeNPY(runKernelTime, ...
            fullfile(folderSession, '_ss_rfRunningKernels.timestamps.npy'));
    end
    
    writeNPY(runningSpeed, fullfile(folderSession, '_ss_running.speed.npy'));
    writeNPY(time_runningSpeed, fullfile(folderSession, ...
        '_ss_running.timestamps.npy'));
    
    writeNPY(pupilSize, fullfile(folderSession, 'eye.diameter.npy'));
    writeNPY(pupilXY, fullfile(folderSession, 'eye.xyPos.npy'));
    writeNPY(time_pupilSize, fullfile(folderSession, ...
        'eye.timestamps.npy'));
end
