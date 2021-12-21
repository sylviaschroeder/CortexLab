%% Parameters
timeGap = 600; %in s, gap between experiments

%% Folders
% UCL PC:
folderTools = 'C:\STORAGE\workspaces';
folderScript = 'C:\dev\workspace\CortexLab';
folderData = '\\zubjects.cortexlab.net\Subjects';
folderBase = 'C:\STORAGE\OneDrive - University of Sussex\';

folderROIData_Sylvia = fullfile(folderBase, 'Lab\DATA\InfoStructs');
folderROIData_Vane = fullfile(folderBase, 'Projects\2021_Vane_Columns\RawData');

%% Add paths
addpath(genpath(fullfile(folderScript)));
addpath('\\ZSERVER.cortexlab.net\Code\2photonPipeline')
addpath('\\ZSERVER.cortexlab.net\Code\Spikes')
addpath('\\ZSERVER.cortexlab.net\Code\Stimulus')
addpath('\\ZSERVER.cortexlab.net\Data\xfiles')
addpath(genpath(fullfile(folderTools, 'npy-matlab', 'npy-matlab')));

%% Load database
db = db_columns_vane;

%% Load corrections
d = load(fullfile(folderROIData_Sylvia, 'corrections_boutons.mat'));
corrections = d.corrections;
doCorrect = d.doCorrect;

%% Convert raw data (2P traces and stimulus info)
experiments = {'expGratings', 'expNoise', 'expStatic', 'expBars'};
expNames = {'gratingsDrifting','sparseNoise', 'gratingsStatic', 'bars'};

for k = 22:length(db)
    fprintf('Dataset %d of %d\n', k, length(db))
    subject = db(k).subject;
    ind = strfind(subject, 'SS');
    subject = subject(ind:end);
    if str2double(subject(end-1:end)) > 50 % boutons
        folderSave = fullfile(folderBase, 'Projects\2021_Vane_Columns\DataToPublish\boutons');
    else
        folderSave = fullfile(folderBase, 'Projects\2021_Vane_Columns\DataToPublish\neurons');
    end
    folderSession = fullfile(folderSave, subject, db(k).date);
    if ~isfolder(folderSession)
        mkdir(folderSession)
    end
    indCorrections = find(strcmp({corrections.subject}, db(k).subject) & ...
        strcmp({corrections.date}, db(k).date));
    
    valid = [];
    xyzPos = [];
    planePerUnit = [];
    cellIDs = [];
    isGad = [];
    infoDone = false;
    F_all = [];
    time = [];
    
    t0 = 0;
    for exp = 1:length(experiments)
        if ~isfield(db, experiments{exp}) || isempty(db(k).(experiments{exp}))
            continue
        end
        
        folder = fullfile(folderROIData_Sylvia, db(k).subject, ...
            db(k).date, num2str(db(k).(experiments{exp})));
        if ~isfolder(folder)
            folder = fullfile(folderROIData_Vane, db(k).subject, ...
                db(k).date, 'metaStructs', num2str(db(k).(experiments{exp})));
        end
        if ~isfolder(folder)
            continue
        end
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
                g = NaN(numCells,1);
                if isfield(meta.ROI, 'isGad')
                    g = meta.ROI.isGad;
                end
                isGad = [isGad; g];
            end
            
            if iPlane == 1
                frameTimes = ppbox.getFrameTimes(meta)';
                frameDur = median(diff(frameTimes));
                time = [time; frameTimes + t0];
                timeShifts = (db(k).planes - db(k).planes(1)) .* ...
                    (frameDur / meta.nPlanes);
                
                switch expNames{exp}
                    case 'gratingsDrifting'
                        [stimTimes, stimSeq] = ...
                            ssLocal.getStimulusResponseInfo(meta);
                        [directions, blank] = gratings.getOrientations(stimSeq);
                        dirs = directions(:,1);
                        dirs(blank) = NaN;
                        
                        writeNPY([stimTimes.onset, stimTimes.offset] + t0, ...
                            fullfile(folderSession, ...
                            sprintf('_ss_%s.intervals.npy', expNames{exp})));
                        writeNPY(stimSeq.seq, ...
                            fullfile(folderSession, ...
                            sprintf('_ss_%s._ss_%sID.npy', expNames{exp}, expNames{exp})));
                        writeNPY(dirs, ...
                            fullfile(folderSession, ...
                            sprintf('_ss_%sID.directions.npy', expNames{exp})));
                    case 'gratingsStatic'
                        [stimTimes, stimSeq] = ...
                            ssLocal.getStimulusResponseInfo(meta);
                        [orientations, spatPhases] = gratings.getStaticGratingPars(stimSeq);
                        
                        writeNPY([stimTimes.onset, stimTimes.offset] + t0, ...
                            fullfile(folderSession, ...
                            sprintf('_ss_%s.intervals.npy', expNames{exp})));
                        writeNPY(stimSeq.seq, ...
                            fullfile(folderSession, ...
                            sprintf('_ss_%s._ss_%sID.npy', expNames{exp}, expNames{exp})));
                        writeNPY(orientations, ...
                            fullfile(folderSession, ...
                            sprintf('_ss_%sID.orientations.npy', expNames{exp})));
                        writeNPY(spatPhases, ...
                            fullfile(folderSession, ...
                            sprintf('_ss_%sID.spatialPhases.npy', expNames{exp})));
                    case 'bars'
                        [stimTimes, stimSeq] = ...
                            ssLocal.getStimulusResponseInfo(meta);
                        [directions, blank] = bars.getOrientations(stimSeq);
                        dirs = directions(:,1);
                        dirs(blank) = NaN;
                        
                        writeNPY([stimTimes.onset, stimTimes.offset] + t0, ...
                            fullfile(folderSession, ...
                            sprintf('_ss_%s.intervals.npy', expNames{exp})));
                        writeNPY(stimSeq.seq, ...
                            fullfile(folderSession, ...
                            sprintf('_ss_%s._ss_%sID.npy', expNames{exp}, expNames{exp})));
                        writeNPY(dirs, ...
                            fullfile(folderSession, ...
                            sprintf('_ss_%sID.directions.npy', expNames{exp})));
                    case 'sparseNoise'
                        stimTimes = ppbox.getStimTimes(meta);
                        [stimFrames, stimPosition] = whiteNoise.getStimulusFrames(meta);
                        stimFrames = permute(stimFrames, [3 1 2]);
                        stimFrameDur = mean(stimTimes.offset - stimTimes.onset) / size(stimFrames,1);
                        stimFrameTimes = ((1:size(stimFrames,1))-1) .* stimFrameDur;
                        allStimFrameTimes = reshape((stimTimes.onset + stimFrameTimes)', [], 1);
                        
                        writeNPY(allStimFrameTimes + t0, ...
                            fullfile(folderSession, ...
                            sprintf('_ss_%s.times.npy', expNames{exp})));
                        writeNPY(repmat((1:size(stimFrames,1))', length(stimTimes.onset), 1), ...
                            fullfile(folderSession, ...
                            sprintf('_ss_%s._ss_%sID.npy', expNames{exp}, expNames{exp})));
                        writeNPY(stimPosition, ...
                            fullfile(folderSession, ...
                            sprintf('_ss_%sArea.edges.npy', expNames{exp})));
                        writeNPY(stimFrames, ...
                            fullfile(folderSession, ...
                            sprintf('_ss_%sID.map.npy', expNames{exp})));
                end
            end
            
            F = meta.F_final;
            
            if ~isempty(indCorrections)
                a = corrections(indCorrections).plane(iPlane).a{db(k).(experiments{exp})};
                b = corrections(indCorrections).plane(iPlane).b{db(k).(experiments{exp})};
                F = doCorrect(a,b,F);
            end
            if size(F_exp,1) > size(F,1)
                F = padarray(F, size(F_exp,1)-size(F,1), NaN, 'post');
            elseif size(F_exp,1)>0 && size(F_exp,1)<size(F,1)
                F(size(F_exp,1)+1:end,:) = [];
            end
            F_exp = [F_exp, F];
            
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
end
