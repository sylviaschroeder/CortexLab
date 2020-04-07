%% Description of useful variables
% caTime            [t x 1], time samples of Ca-traces, can be
%                   discontinuous if multiple experiments combined
% caTraces          [t x n], Ca-traces of all simultaneously recorded
%                   neurons
% cellIDs           [n x 2], [plane, cellID] for each recorded neuron
% tlTime            [t2 x 1], time samples of rotary encoder (wheel
%                   movements), higher resolution than caTime
% wheelVelocity     [t2 x 1], wheel velocity, if >0 movement to the
%                   right/chooses left, if <0 movement to the left/chosses
%                   right
% allMoves          [t3 x 2], [start stop] times of wheel movements
% afterStimMoves    [trial x 1], start of first wheel movements after stimulus
%                   onset for each trial, NaN if there was no movement
% stimOnTimes       [trial x 1], onset times of visual stimulus
% beepTimes         [trial x 1], times of go cue
% feedbackTimes     [trial x 1], times of feedback, water reward if correct,
%                   auditory noise if incorrect
% stimulus          [trial x 2], contrast of stimulus on left (1st column)
%                   and right (2nd column) side for each trial
% repeated          [trial x 1], true if trial condition is repeated
%                   (because previous trial was incorrect)
% choice            [trial x 1], 1: chose left, 2: chose right, 3: nogo
% outcome           [trial x 1], true if choice was correct

%% Folders and parameters
metaFolder = '\\zarchive.cortexlab.net\Data\Subjects';
timelineFolder = '\\ZSERVER.cortexlab.net\Data\expInfo';
configFolder = '\\ZSERVER.cortexlab.net\Code\Rigging\config';

addpath('\\ZSERVER.cortexlab.net\Code\2photonPipeline')

smoothing = .15; % s
velocityThreshold = 10; % mm/s
minBetweenMoves = .2; % s
reactionTime = .05; % s

%% Datasets
db_wheelTask;

for iSet = 1:length(db)
    fprintf('Dataset %d: %s %s\n', iSet, db(iSet).subject, db(iSet).date)
    
    %% Load data
    caTime = [];
    caTraces = [];
    tlTime = [];
    rotaryEncoder = [];
    stimOnTimes = [];
    stimOffTimes = [];
    beepTimes = [];
    feedbackTimes = [];
    stimulus = [];
    repeated = [];
    choice = [];
    outcome = [];
    t0 = 0;
    valid = [];
    for iExp = 1:length(db(iSet).exp)
        % load meta structures
        folder = fullfile(metaFolder, db(iSet).subject, db(iSet).date, ...
            'metaStructs', num2str(iExp));
        files = dir(fullfile(folder, '*_ROI.mat'));
        val = [];
        for i = 1:length(files) % different planes
            data = load(fullfile(folder, files(i).name));
            meta = data.meta;
            t = ppbox.getFrameTimes(meta);
            j = meta.ROI.isDuplicate==0 & meta.ROI.isSwitchOn==0;
            val = [val; j];
            if i == 1
                caTime_exp = t;
                caTraces_exp = meta.F_final;
                if iExp == 1
                    cellIDs = [ones(length(j),1) .* meta.iPlane, (1:length(j))'];
                end
            else
                indNaN = all(isnan(meta.F_final(:,j)),2);
                tr = NaN(length(caTime_exp), length(j));
                warning off
                tr(:,j) = interp1(t, meta.F_final(:,j), caTime_exp, 'pchip');
                warning on
                diffs = diff(indNaN);
                starts = find(diffs==1)+1;
                stops = find(diffs==-1)+1;
                if ~isempty(starts) || ~isempty(stops)
                    if length(starts)>length(stops)
                        stops(end+1) = length(indNaN);
                    end
                    if length(starts)<length(stops)
                        starts = [1 starts];
                    end
                    if starts(1)>stops(1)
                        starts = [1 starts];
                    end
                    if starts(end)>stops(end)
                        stops(end+1) = length(indNaN);
                    end
                    for k = 1:length(starts)
                        t1 = find(caTime_exp >= t(starts(k)), 1);
                        t2 = find(caTime_exp >= t(stops(k)), 1);
                        if t1 == 2, t1 = 1; end
                        if isempty(t2), t2 = length(caTime_exp); end
                        tr(t1:t2,:) = NaN;
                    end
                end
                %             ind2NaN = hist(t(indNaN), caTime_exp) > 0;
                %             tr(ind2NaN,:) = NaN;
                caTraces_exp = [caTraces_exp, tr];
                if iExp == 1
                    cellIDs = [cellIDs; ones(length(j),1) .* ...
                        meta.iPlane, (1:length(j))'];
                end
            end
        end
        valid = [valid, val];
        data = load(fullfile(folder, 'timeAlign.mat'));
        blockAlign = data.alignment;
        
        % load timeline (time stamps, rotary encoder)
        folder = fullfile(timelineFolder, db(iSet).subject, db(iSet).date, ...
            num2str(db(iSet).exp(iExp)));
        file = dir(fullfile(folder, '*_Timeline.mat'));
        if length(file) ~= 1
            disp('Timeline file missing or more than 1 exists!');
            return
        end
        data = load(fullfile(folder, file.name));
        timeline = data.Timeline;
        tlTime_exp = timeline.rawDAQTimestamps;
        rotaryEncoder_exp = timeline.rawDAQData(:, strcmp({timeline.hw.inputs.name}, ...
            'rotaryEncoder'));
        
        % block: stim times, conditions, beep times, choice times
        file = dir(fullfile(folder, '*_Block.mat'));
        if length(file) ~= 1
            disp('Block file missing or more than 1 exists!');
            return
        end
        data = load(fullfile(folder, file.name));
        block_exp = data.block;
        block_exp = beh.stripIncompleteTrials(block_exp);
        stimOnTimes_exp = blockAlign.stimOnTimes';
        beepTimes_exp = blockAlign.blockToPdTimeFrame([block_exp.trial.interactiveStartedTime])';
        feedbackTimes_exp = blockAlign.blockToPdTimeFrame([block_exp.trial.feedbackStartedTime])';
        stimulus_exp = zeros(length(block_exp.trial), 2);
        repeated_exp = false(length(block_exp.trial), 1);
        for i = 1:length(block_exp.trial)
            stimulus_exp(i,:) = block_exp.trial(i).condition.visCueContrast;
            repeated_exp(i) = block_exp.trial(i).condition.repeatNum > 1;
        end
        choice_exp = cat(1, block_exp.trial.responseMadeID);
        outcome_exp = cat(1, block_exp.trial.feedbackType) == 1;
        
        caTime = [caTime; caTime_exp'+t0];
        caTraces = [caTraces; caTraces_exp];
        tlTime = [tlTime; tlTime_exp'+t0];
        rotaryEncoder = [rotaryEncoder; rotaryEncoder_exp];
        stimOnTimes = [stimOnTimes; stimOnTimes_exp+t0];
        beepTimes = [beepTimes; beepTimes_exp+t0];
        feedbackTimes = [feedbackTimes; feedbackTimes_exp+t0];
        stimulus = [stimulus; stimulus_exp];
        repeated = [repeated; repeated_exp];
        choice = [choice; choice_exp];
        outcome = [outcome; outcome_exp];
        t0 = t0 + max(caTime_exp(end), tlTime_exp(end)) + .001;
    end
    valid = all(valid,2);
    caTraces = caTraces(:, valid);
    cellIDs = cellIDs(valid,:);
    outcome = logical(outcome);

    %% Get wheel movements
    switch db(iSet).microID
        case 'b'
            rigPC = 'ZMAZE';
            hardwarefile = 'hardware_20170519';
        case 'b2'
            rigPC = 'ZURPRISE';
            hardwarefile = 'hardware20170131.mat';
    end
    data = load(fullfile(configFolder, rigPC, hardwarefile), 'mouseInput');
    mmFactor = data.mouseInput.MillimetresFactor;
    
    rawPos = wheel.correctCounterDiscont(rotaryEncoder);
    sr = 1 / median(diff(tlTime));
    time = round(tlTime(1)*sr)/sr : 1/sr : round(tlTime(end)*sr)/sr;
    pos = interp1(tlTime, rawPos, time, 'pchip');
    pos = pos .* mmFactor;
    wheelVelocity = wheel.computeVelocity(pos, round(smoothing * sr), sr)';
    [allMoves, allAmps, peakVels] = wheel.findAllMoves(time, wheelVelocity, pos, velocityThreshold, ...
        minBetweenMoves); % velocity >0 if choosing left (wheel turned to the right),
    allMoves = allMoves';
    % velocity <0 if choosing right (wheel turned to the left)
    allOnsets = allMoves(:,1);
    startTimes = sort([stimOnTimes(~all(stimulus==0, 2)); ...
        beepTimes(all(stimulus==0, 2))] + reactionTime); % stim onset or beep (if contrast 0) + minimum reaction time
    moves = repmat(allOnsets, 1, length(startTimes));
    moves(bsxfun(@gt, moves, feedbackTimes')) = NaN;
    afterStim = bsxfun(@minus, moves, startTimes');
    afterStim(afterStim < 0) = NaN;
    [t,ind] = min(afterStim, [], 1);
    afterStimMoves = allOnsets(ind);
    afterStimMoves(isnan(t)) = NaN;
    
    variables(iSet).caTime = caTime;
    variables(iSet).caTraces = caTraces;
    variables(iSet).cellIDs = cellIDs;
    variables(iSet).tlTime = tlTime;
    variables(iSet).wheelVelocity = wheelVelocity;
    variables(iSet).allMoves = allMoves;
    variables(iSet).afterStimMoves = afterStimMoves;
    variables(iSet).stimOnTimes = stimOnTimes;
    variables(iSet).beepTimes = beepTimes;
    variables(iSet).feedbackTimes = feedbackTimes;
    variables(iSet).stimulus = stimulus;
    variables(iSet).repeated = repeated;
    variables(iSet).choice = choice;
    variables(iSet).outcome = outcome;
end