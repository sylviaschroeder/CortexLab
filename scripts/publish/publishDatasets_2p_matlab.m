% units = 'boutons';
units = 'neurons';

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
    folderKernels = fullfile(folderBase, 'RESULTS\boutons\nonVisualEffects\kernelFit');
    d = load(fullfile(folderROIData, 'corrections_boutons.mat'));
    corrections = d.corrections;
    doCorrect = d.doCorrect;
    folderSave = fullfile(folderBase, 'DATA\DataToPublish\boutons');
else
    folderKernels = fullfile(folderBase, 'RESULTS\nonVisualEffects\modelGratingResp\kernelFit');
    corrections = [];
    folderSave = fullfile(folderBase, 'DATA\DataToPublish\sc neurons 2p');
end

%% Collect data
experiments = {'expGratings', 'expGrayScreen', 'expDark', 'expNoise'};
expNames = {'gratings','grayScreen','dark','noise'};

d = load(fullfile(folderKernels, 'results.mat'));
results = d.results;

for k = 1:length(db)
    fprintf('Dataset %d of %d\n', k, length(db))
    data = [];
    subject = db(k).subject;
    ind = strfind(subject, 'SS');
    subject = subject(ind:end);
    data.subject = subject;
    data.date = db(k).date;
    data.units.plane = [];
    data.units.xPos = [];
    data.units.yPos = [];
    data.units.zPos = [];
    data.units.isGad = [];
    
    valid = [];
    infoDone = false;
    for exp = 1:length(experiments)
        if ~isfield(db, experiments{exp}) || isempty(db(k).(experiments{exp}))
            continue
        end
        data.(expNames{exp}).F = [];
        runningSpeed = [];
        time_runningSpeed = [];
        pupilSize = [];
        time_pupilSize = [];
        
        folder = fullfile(folderROIData, db(k).subject, ...
            db(k).date, num2str(db(k).(experiments{exp})));
        file = [sprintf('%s_%d_%s', db(k).date, db(k).(experiments{exp}), ...
            db(k).subject) '_2P_plane%03d_ROI.mat'];
        valid_exp = [];
        for iPlane = 1:length(db(k).planes)
            % load data
            d = load(fullfile(folder, sprintf(file, db(k).planes(iPlane))));
            meta = d.meta;
            meta.folderTL = strrep(meta.folderTL, 'ioo.ucl.ac.uk', ...
                'cortexlab.net');
            numCells = size(meta.F_final,2);
            
            if ~infoDone
                xyz = cat(1, meta.ROI.CellXYZMicrons{:});
                data.units.plane = [data.units.plane; ones(numCells,1) .* iPlane];
                data.units.xPos = [data.units.xPos; xyz(:,1)];
                data.units.yPos = [data.units.yPos; xyz(:,2)];
                data.units.zPos = [data.units.zPos; xyz(:,3)];
                if strcmp(units, 'boutons')
                    g = NaN(size(xyz,1),1);
                else
                    g = meta.ROI.isGad;
                end
                data.units.isGad = [data.units.isGad; g];
            end
            
            if iPlane == 1
                ballData = nonVis.getRunningSpeed(meta);
                if ~isempty(ballData)
                    runningSpeed = ballData.total / median(diff(ballData.t)) / 53;
                    time_runningSpeed = ballData.t;
                end
                [pupilData, t] = nonVis.loadPupilData(meta);
                if ~isempty(pupilData) && ~strcmp(expNames{exp}, 'dark')
                    pupilSize = nonVis.getPupilDiam(pupilData);
                    time_pupilSize = t(1:length(pupilData.x));
                end
                frameTimes = ppbox.getFrameTimes(meta);
                frameDur = median(diff(frameTimes));
                timeShifts = (db(k).planes - db(k).planes(1)) .* ...
                    (frameDur / meta.nPlanes);
                
                data.(expNames{exp}).time_F = frameTimes;
                data.(expNames{exp}).runningSpeed = runningSpeed;
                data.(expNames{exp}).time_runningSpeed = time_runningSpeed;
                data.(expNames{exp}).pupilSize = pupilSize;
                data.(expNames{exp}).time_pupilSize = time_pupilSize;
                data.timeShiftPerPlane = timeShifts;
                
                switch expNames{exp}
                    case 'gratings'
                        [stimTimes, stimSeq, stimMatrix] = ...
                            ssLocal.getStimulusResponseInfo(meta);
                        [directions, blank] = gratings.getOrientations(stimSeq);
                        data.gratings.stimuli.onTimes = stimTimes.onset;
                        data.gratings.stimuli.offTimes = stimTimes.offset;
                        data.gratings.stimuli.sequence = stimSeq.seq;
                        dirs = directions(:,1);
                        dirs(blank) = NaN;
                        data.gratings.stimuli.directions = dirs;
                        data.gratings.responseAmplitudes = [];
                        data.gratings.kernels = [];
                    case 'noise'
                        stimTimes = ppbox.getStimTimes(meta);
                        [stimFrames, stimPosition] = whiteNoise.getStimulusFrames(meta);
                        stimFrameDur = mean(stimTimes.offset - stimTimes.onset) / size(stimFrames,3);
                        stimFrameTimes = ((1:size(stimFrames,3))-1) .* stimFrameDur;
                        allStimFrameTimes = reshape((stimTimes.onset + stimFrameTimes)', [], 1);
                        data.noise.stimuli.onTimes = allStimFrameTimes;
                        data.noise.stimuli.position = stimPosition;
                        data.noise.stimuli.sequence = repmat((1:size(stimFrames,3))', ...
                            length(stimTimes.onset), 1);
                        data.noise.stimuli.frames = stimFrames;
                end
            end
            
            F = meta.F_final;
            if ~isempty(corrections)
                a = corrections(k).plane(iPlane).a{db(k).(experiments{exp})};
                b = corrections(k).plane(iPlane).b{db(k).(experiments{exp})};
                F = doCorrect(a,b,F);
            end
            if size(data.(expNames{exp}).F,1) > size(F,1)
                F = padarray(F, size(data.(expNames{exp}).F,1)-size(F,1), ...
                    NaN, 'post');
            elseif size(data.(expNames{exp}).F,1)>0 && size(data.(expNames{exp}).F,1)<size(F,1)
                F(size(data.(expNames{exp}).F,1)+1:end,:) = [];
            end
            data.(expNames{exp}).F = [data.(expNames{exp}).F, F];
            
            if strcmp(expNames{exp}, 'gratings')
                amps = NaN(size(results(k).plane(iPlane).responses,3), ...
                    length(dirs), size(meta.F_final,2));
                kernels = NaN(length(results(k).plane(iPlane).kernelTime), ...
                    size(meta.F_final,2));
                for iCell = 1:length(results(k).plane(iPlane).cellIDs)
                    a = results(k).plane(iPlane).kernelFit(iCell).alphaEachTrial;
                    if isempty(a)
                        continue
                    end
                    amps(:,:,results(k).plane(iPlane).cellIDs(iCell)) = a;
                    kernels(:,results(k).plane(iPlane).cellIDs(iCell)) = ...
                        results(k).plane(iPlane).kernelFit(iCell).kernel;
                end
                if ~isempty(corrections)
                    a = corrections(k).plane(iPlane).a{db(k).(experiments{exp})};
                    a = permute(a, [1 3 2]);
                    b = corrections(k).plane(iPlane).b{db(k).(experiments{exp})};
                    b = permute(b, [1 3 2]);
                    amps = doCorrect(a,b,amps);
                end
                data.gratings.responseAmplitudes = ...
                    cat(3, data.gratings.responseAmplitudes, amps);
                data.gratings.kernels = [data.gratings.kernels, kernels];
                data.gratings.time_kernel = results(k).plane(iPlane).kernelTime;
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
    end
    sum(valid)
%     valid = all(valid,2);
    valid = any(valid,2);
    sum(valid)
    data.units.plane(~valid) = [];
    data.units.xPos(~valid) = [];
    data.units.yPos(~valid) = [];
    data.units.zPos(~valid) = [];
    data.units.isGad(~valid) = [];
    for exp = 1:length(experiments)
        if ~isfield(db, experiments{exp}) || isempty(db(k).(experiments{exp}))
            continue
        end
        data.(expNames{exp}).F(:,~valid) = [];
        if strcmp(expNames{exp},'gratings')
            data.gratings.responseAmplitudes(:,:,~valid) = [];
            data.gratings.kernels(:,~valid) = [];
        end
    end
    save(fullfile(folderSave, sprintf('%s_%s.mat', subject, ...
        db(k).date)), 'data')
end