%% Folders and parameters
metaFolder = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\InfoStructs';
timelineFolder = '\\ZSERVER.cortexlab.net\Data\expInfo';
resultsFolder = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\wheelTask\kernelFit';
traceFolder = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\wheelTask\eventTraces';
wheelFolder = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\wheelTask\wheelMovements';
performanceFolder = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\wheelTask\performance';
configFolder = '\\ZSERVER.cortexlab.net\Code\Rigging\config';
RFfolder = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\wheelTask\receptiveFields';

velocityThreshold = 10; % mm/s
smoothing = .15; % s
minBetweenMoves = .2; % s
reactionTime = .05; % s

% kernel windows (in sec)
stimWin_krnl = [0 .8];
moveStimWin_krnl = [-.4 .4];
rewardWin_krnl = [0 3];
moveWin_krnl = [-.4 .8];
stimWin_eta = [-.5 3];
moveStimWin_eta = [-.5 2];
rewardWin_eta = [-2 3];
moveWin_eta = [-1 2];

% kernel regression
lambda = .01;
kernelNames = {'choice L corr', 'choice L incorr', 'choice R corr', ...
    'choice R incorr', 'choice NoGo corr', 'choice NoGo incorr',...
    'stim L home', 'stim R home', ...
    'stim L to front', 'stim L to back', 'stim R to front', 'stim R to back', ...
    'move wheel left', 'move wheel right'};

% cross-validation
partition = 50; % k-fold crossvalidation
modelNames = {{'-','stim home, pres','stim home, contrast','stim home, sqrt contrast'}, ...
    {'-','stim moves, pres','stim moves, contrast','stim moves, sqrt contrast'}, ...
    {'-','nonvisual'}};

%% Datasets
db_wheelTask;

% colormap
cm = ones(129, 3);
w = linspace(0, 1, 65)';
w(1) = [];
cm(1:64,:) = [0 0 1].*flip(w) + [1 1 1].*(1-flip(w));
cm(66:end,:) = [1 0 0].*w + [1 1 1].*(1-w);

for iSet = 1:length(db)
    fprintf('Dataset %d: %s %s\n', iSet, db(iSet).subject, db(iSet).date)
    
%     results(iSet).subject = db(iSet).subject;
%     results(iSet).date = db(iSet).date;
%     results(iSet).exp = db(iSet).exp;
    
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
        folder = fullfile(metaFolder, db(iSet).subject, db(iSet).date, num2str(iExp));
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
                    isGad = meta.ROI.isGad;
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
                    isGad = [isGad; meta.ROI.isGad];
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
        block_exp = psy.stripIncompleteTrials(block_exp);
        stimOnTimes_exp = blockAlign.stimOnTimes;
        stimOffTimes_exp = blockAlign.stimOffTimes;
        beepTimes_exp = blockAlign.blockToPdTimeFrame([block_exp.trial.interactiveStartedTime]);
        feedbackTimes_exp = blockAlign.blockToPdTimeFrame([block_exp.trial.feedbackStartedTime]);
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
        stimOnTimes = [stimOnTimes, stimOnTimes_exp+t0];
        stimOffTimes = [stimOffTimes, stimOffTimes_exp+t0];
        beepTimes = [beepTimes, beepTimes_exp+t0];
        feedbackTimes = [feedbackTimes, feedbackTimes_exp+t0];
        stimulus = [stimulus; stimulus_exp];
        repeated = [repeated; repeated_exp];
        choice = [choice; choice_exp];
        outcome = [outcome; outcome_exp];
        t0 = t0 + max(caTime_exp(end), tlTime_exp(end)) + .001;
        if iExp==1
            block = block_exp;
        else
            block = [block, block_exp];
        end
    end
    valid = all(valid,2);
    caTraces = caTraces(:, valid);
    isGad = isGad(valid);
    cellIDs = cellIDs(valid,:);
    outcome = logical(outcome);
    
    results(iSet).cellIDs = cellIDs;
    results(iSet).isGad = isGad;
    results(iSet).caTraces = caTraces;
    results(iSet).time = caTime;

    clear iExp folder files i data meta t ind caTime_exp caTraces_exp indNaN tr
    clear diffs starts stops k t1 t2 blockAlign timeline tlTime_exp
    clear rotaryEncoder_exp block_exp stimOnTimes_exp stimOffTimes_exp
    clear beepTimes_exp feedbackTimes_exp stimulus_exp repeated_exp choice_exp
    clear outcome_exp t0

    %% Plot psychometric function and reaction times
    folder = fullfile(performanceFolder, sprintf('%s_%s_%d', db(iSet).subject, ...
        db(iSet).date, db(iSet).exp(1)));
    if ~isdir(folder)
        mkdir(folder)
    end
    % find valid trials (calcium responses not NaN)
    indNaN = sum(isnan(caTraces), 2) > size(caTraces,2)/2;
    diffs = diff(indNaN);
    starts = find(diffs==1)+1;
    stops = find(diffs==-1)+1;
    if ~isempty([starts; stops])
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
    end
    invalidTrials = find(any((stimOnTimes>caTime(starts) & stimOnTimes<caTime(stops)) | ...
        (feedbackTimes>caTime(starts) & feedbackTimes<caTime(stops)), 1));
    bl = cell(1, length(block));
    k = 0;
    for i = 1:length(bl)
        bl{i} = block(i);
        ind = invalidTrials(invalidTrials>k & invalidTrials>k+length(bl{i}.trial)) - k;
        bl{i}.trial(ind) = [];
        bl{i}.numCompletedTrials = bl{i}.numCompletedTrials-length(ind);
        rew = [block(i).trial(ind).feedbackType] == 1;
        bl{i}.rewardDeliveredSizes(ind(rew)) = [];
        bl{i}.rewardDeliveryTimes(ind(rew)) = [];
        k = k + length(block(i).trial);
    end
    [~,f] = beh.generatePsychometric(bl, 'none');
    f(1).PaperPositionMode = 'auto';
    figure(f(1));
    print(fullfile(folder, 'psychometric.tiff'), '-dtiff', '-r0')
    f(2).PaperPositionMode = 'auto';
    figure(f(2));
    print(fullfile(folder, 'reactionTimes.tiff'), '-dtiff', '-r0')
    
    %%   close figures
%     close(f)
%     
%     clear folder indNaN diffs starts stops bl ind f

    %% Get wheel movements
    switch db(iSet).microID
        case 'b'
            rigPC = 'ZMAZE';
        case 'b2'
            rigPC = 'ZURPRISE';
    end
    data = load(fullfile(configFolder, rigPC, 'hardware.mat'), 'mouseInput');
    mmFactor = data.mouseInput.MillimetresFactor;
    
    rawPos = wheel.correctCounterDiscont(rotaryEncoder);
    sr = 1 / median(diff(tlTime));
    time = round(tlTime(1)*sr)/sr : 1/sr : round(tlTime(end)*sr)/sr;
    pos = interp1(tlTime, rawPos, time, 'pchip');
    pos = pos .* mmFactor;
    [vel, acc] = wheel.computeVelocity(pos, round(smoothing * sr), sr);
    [allMoves, allAmps, peakVels] = wheel.findAllMoves(time, vel, pos, velocityThreshold, ...
        minBetweenMoves); % velocity >0 if choosing left (wheel turned to the right),
    % velocity <0 if choosing right (wheel turned to the left)
    allOnsets = allMoves(1,:);
    startTimes = sort([stimOnTimes(~all(stimulus==0, 2)), ...
        beepTimes(all(stimulus==0, 2))] + reactionTime); % stim onset or beep (if contrast 0) + minimum reaction time
    moves = repmat(allOnsets', 1, length(startTimes));
    moves(bsxfun(@gt, moves, feedbackTimes)) = NaN;
    afterStim = bsxfun(@minus, moves, startTimes);
    afterStim(afterStim < 0) = NaN;
    [t,ind] = min(afterStim, [], 1);
    afterStimMoves = allOnsets(ind);
    afterStimMoves(isnan(t)) = NaN;
    afterStimAmps = allAmps(ind);
    afterStimAmps(isnan(t)) = NaN;
    
    clear rigPC data mmFactor rawPos acc allAmps startTimes moves afterStim t ind

    %% Plot wheel movements (absolute position and velocity)
    folder = fullfile(wheelFolder, sprintf('%s_%s_%d', db(iSet).subject, ...
        db(iSet).date, db(iSet).exp(1)));
    if ~isdir(folder)
        mkdir(folder)
    end
    clear f
    yLimPos = [-30 30];
    yLimVel = [-300 300];
    xLimPos = [-.5 1.5];
    xLimVel = [-.2 1];
    % (1) Movements between stim onset and beep
    ind = any(stimulus>0,2);
    durs = beepTimes(ind)-stimOnTimes(ind);
    win = [-1 max(durs)];
    t = stimOnTimes(ind)';
    nSamples = round(win .* sr);
    nSamples = nSamples(1) : nSamples(end);
    moves = NaN(length(t), length(nSamples));
    vels = NaN(length(t), length(nSamples));
    for k = 1:length(t)
        ind = round(t(k) * sr) + (nSamples(1) : round(durs(k) * sr));
        moves(k,1:length(ind)) = smooth(pos(ind), 51);
        vels(k,1:length(ind)) = smooth(vel(ind), 51);
    end
    moves = moves - moves(:,nSamples==0);
    figure
    subplot(1,2,1)
    plot((nSamples(1):nSamples(end))./sr, moves, 'k')
    hold on
    plot([0 0],yLimPos, 'k:')
    xlim(xLimPos)
    ylim(yLimPos)
    xlabel('Time from stim onset (s)')
    ylabel('Wheel position (mm)')
    title(sprintf('Early movements (between stim onset and beep, n = %d)', length(t)))
    set(gca, 'box', 'off')
    subplot(1,2,2)
    plot((nSamples(1):nSamples(end))./sr, vels, 'k')
    hold on
    plot([0 0],yLimVel, 'k:')
    xlim(xLimVel)
    ylim(yLimVel)
    xlabel('Time from stim onset (s)')
    ylabel('Wheel velocity (mm/s)')
    % title('Early movements (between stim onset and beep)')
    set(gca, 'box', 'off')
    set(gcf, 'Position', [5 400 940 690])
    f(1) = gcf;
    f(1).PaperPositionMode = 'auto';
    print(fullfile(folder, 'early.tiff'), '-dtiff', '-r0')
    
    % (2) Movements between beep and feedback (choice was NOT NoGo)
    ind = choice < 3;
    durs = feedbackTimes(ind)-beepTimes(ind);
    win = [-.5 max(durs)];
    t = beepTimes(ind)';
    nSamples = round(win .* sr);
    nSamples = nSamples(1) : nSamples(end);
    moves = NaN(length(t), length(nSamples));
    vels = NaN(length(t), length(nSamples));
    for k = 1:length(t)
        ind = round(t(k) * sr) + (nSamples(1) : round(durs(k) * sr));
        moves(k,1:length(ind)) = smooth(pos(ind), 51);
        vels(k,1:length(ind)) = smooth(vel(ind), 51);
    end
    moves = moves - moves(:,nSamples==0);
    figure
    subplot(1,2,1)
    hold on
    cols = lines(2);
    cols = cols([2 1],:);
    h = zeros(1,2);
    for k = 1:2
        h_ = plot((nSamples(1):nSamples(end))'./sr, ...
            moves(choice(choice<3)==k,:), 'Color', cols(k,:));
        h(k) = h_(1);
    end
    xlim(xLimPos)
    ylim(yLimPos)
    legend(h, 'left', 'right')
    xlabel('Time from beep (s)')
    ylabel('Wheel position (mm)')
    title(sprintf('Response movements (choice: left or right, n = %d)', length(t)))
    set(gca, 'box', 'off')
    subplot(1,2,2)
    hold on
    cols = lines(2);
    cols = cols([2 1],:);
    h = zeros(1,2);
    for k = 1:2
        h_ = plot((nSamples(1):nSamples(end))'./sr, ...
            vels(choice(choice<3)==k,:), 'Color', cols(k,:));
        h(k) = h_(1);
    end
    xlim(xLimVel)
    ylim(yLimVel)
    legend(h, 'left', 'right')
    xlabel('Time from beep (s)')
    ylabel('Wheel velocity (mm/s)')
    set(gca, 'box', 'off')
    set(gcf, 'Position', [950 400 940 690])
    f(2) = gcf;
    f(2).PaperPositionMode = 'auto';
    print(fullfile(folder, 'afterBeep.tiff'), '-dtiff', '-r0')
    
    % (3) Movements leading to decision (onset just before choice)
    indTrial = choice < 3;
    choiceTimes = feedbackTimes(indTrial);
    indMoves = sum(bsxfun(@lt, allOnsets', choiceTimes), 1);
    durs = choiceTimes - allOnsets(indMoves);
    win = [-.5 max(durs)];
    t = allOnsets(indMoves)';
    nSamples = round(win .* sr);
    nSamples = nSamples(1) : nSamples(end);
    moves = NaN(length(t), length(nSamples));
    vels = NaN(length(t), length(nSamples));
    for k = 1:length(t)
        ind = round(t(k) * sr) + (nSamples(1) : round(durs(k) * sr));
        moves(k,1:length(ind)) = smooth(pos(ind), 51);
        vels(k,1:length(ind)) = smooth(vel(ind), 51);
    end
    moves = moves - moves(:,nSamples==0);
    figure
    subplot(1,2,1)
    hold on
    cols = lines(2);
    cols = cols([2 1],:);
    h = zeros(1,2);
    for k = 1:2
        h_ = plot((nSamples(1):nSamples(end))'./sr, ...
            moves(choice(indTrial)==k,:), 'Color', cols(k,:));
        h(k) = h_(1);
    end
    xlim(xLimPos)
    ylim(yLimPos)
    legend(h, 'left', 'right')
    xlabel('Time from movement onset (s)')
    ylabel('Wheel position (mm)')
    title(sprintf('Choice movements (leading to choice, n = %d)', length(t)))
    set(gca, 'box', 'off')
    subplot(1,2,2)
    hold on
    cols = lines(2);
    cols = cols([2 1],:);
    h = zeros(1,2);
    for k = 1:2
        h_ = plot((nSamples(1):nSamples(end))'./sr, ...
            vels(choice(indTrial)==k,:), 'Color', cols(k,:));
        h(k) = h_(1);
    end
    xlim(xLimVel)
    ylim(yLimVel)
    legend(h, 'left', 'right')
    xlabel('Time from movement onset (s)')
    ylabel('Wheel velocity (mm/s)')
    set(gca, 'box', 'off')
    set(gcf, 'Position', [1930 400 940 690])
    f(3) = gcf;
    f(3).PaperPositionMode = 'auto';
    print(fullfile(folder, 'leadingToChoice.tiff'), '-dtiff', '-r0')
    
    % (4) Other movements (starting before stim onset or after feedback)
    indMoves = [find(allOnsets < stimOnTimes(1)), ...
        find(sum(bsxfun(@gt, allOnsets', feedbackTimes(1:end-1)) & ...
        bsxfun(@lt, allOnsets', stimOnTimes(2:end)),2)'>0), ...
        find(allOnsets > feedbackTimes(end))];
    durs = diff(allMoves(:,indMoves),1,1);
    win = [-.5 max(durs)];
    t = allOnsets(indMoves)';
    nSamples = round(win .* sr);
    nSamples = nSamples(1) : nSamples(end);
    moves = NaN(length(t), length(nSamples));
    vels = NaN(length(t), length(nSamples));
    for k = 1:length(t)
        ind = round(t(k) * sr) + (nSamples(1) : round(durs(k) * sr));
        j = ind>0;
        j(end+1:length(nSamples)) = false;
        moves(k,j) = smooth(pos(ind(j)), 51);
        vels(k,j) = smooth(vel(ind(j)), 51);
    end
    moves = moves - moves(:,nSamples==0);
    figure
    subplot(1,2,1)
    plot((nSamples(1):nSamples(end))./sr, moves, 'k')
    xlim(xLimPos)
    ylim(yLimPos)
    xlabel('Time from movement onset (s)')
    ylabel('Wheel position (mm)')
    title(sprintf('Out-of-trial movements (n = %d)', length(t)))
    set(gca, 'box', 'off')
    subplot(1,2,2)
    plot((nSamples(1):nSamples(end))./sr, vels, 'k')
    xlim(xLimVel)
    ylim(yLimVel)
    xlabel('Time from movement onset (s)')
    ylabel('Wheel velocity (mm/s)')
    set(gca, 'box', 'off')
    set(gcf, 'Position', [2875 400 940 690])
    f(4) = gcf;
    f(4).PaperPositionMode = 'auto';
    print(fullfile(folder, 'others.tiff'), '-dtiff', '-r0')

    %%   close figures
%     close(f)
%     
%     clear folder yLimPos yLimVel xLimPos xLimVel ind durs win t nSamples moves
%     clear vels k h h_ cols indTrial choiceTimes indMoves
    
    %% Plot event aligned traces for each cell
    % Define different movements
    moves = repmat(allOnsets', 1, length(stimOnTimes));
    indStimLpostBeep = sum(moves > beepTimes & moves < feedbackTimes & ...
        stimulus(:,1)' > 0, 2) > 0; % stim on left, bettween beep and choice
    indStimLpreBeep = sum(moves > stimOnTimes & moves < beepTimes & ...
        stimulus(:,1)' > 0, 2) > 0; % stim on left, bettween stim onset and beep
    indStimRpostBeep = sum(moves > beepTimes & moves < feedbackTimes & ...
        stimulus(:,2)' > 0, 2) > 0; % stim on right, bettween beep and choice
    indStimRpreBeep = sum(moves > stimOnTimes & moves < beepTimes & ...
        stimulus(:,2)' > 0, 2) > 0; % stim on right, bettween stim onset and beep
    indnoStimpostBeep = sum(moves > beepTimes & moves < feedbackTimes & ...
        all(stimulus==0,2)', 2) > 0; % no stim, between beep and choice
    indOthers = ~(indStimLpostBeep | indStimLpreBeep | indStimRpostBeep | ...
        indStimRpreBeep | indnoStimpostBeep); % all other movements (between choice and next stim (beep for nogo))
    % Define events windows
    stimWin = [-1 3];
    beepWin = [-1 2];
    choiceWin = [-2 2];
    moveWin = [-1 2];
    % Define event conditions
    contrastL = setdiff(unique(stimulus(:,1)), 0);
    contrastR = setdiff(unique(stimulus(:,2)), 0);
    contrastIndsL = cell(1, length(contrastL));
    contrastIndsR = cell(1, length(contrastR));
    for k = 1:length(contrastL)
        contrastIndsL{k} = stimulus(:,1) == contrastL(k);
    end
    for k = 1:length(contrastR)
        contrastIndsR{k} = stimulus(:,2) == contrastR(k);
    end
    % conditions: {{column plots {single row plot {division within plot} } }}
%     conditions = {{ ... % first column (stim)
%         contrastIndsL, ... % different contrasts on left
%         contrastIndsR}, ... % different contrasts on right
%         { ... % second column (beep)
%         {stimulus(:,1)>0, ...% stim on left
%         stimulus(:,2)>0, ...% stim on right
%         all(stimulus==0,2)}}, ...% no stim
%         { ... % third colum (choice)
%         {choice==1 & outcome, ...% time of choice to left, correct choice
%         choice==1 & ~outcome}, ...% time of choice to left, incorrect choice
%         {choice==2 & outcome, ...% time of choice to right, correct choice
%         choice==2 & ~outcome}, ...% time of choice to right, incorrect choice
%         {choice==3 & outcome, ...% time of choice for NoGo, correct choice
%         choice==3 & ~outcome}}, ...% time of choice for NoGo, incorrect choice
%         { ... % fourth column (left movements)
%         {peakVels'>0 & indStimLpostBeep, ...% left move, left stim, after beep
%         peakVels'>0 & indStimLpreBeep, ...% left move, left stim, before beep
%         peakVels'>0 & indStimRpostBeep, ...% left move, right stim, after beep
%         peakVels'>0 & indStimRpreBeep, ...% left move, right stim, before beep
%         peakVels'>0 & indnoStimpostBeep, ...% left move, no stim, after beep
%         peakVels'>0 & indOthers}}, ...% left move, outside trial
%         { ... % fifth column (right movements)
%         {peakVels'<0 & indStimLpostBeep, ...% right move, left stim, after beep
%         peakVels'<0 & indStimLpreBeep, ...% right move, left stim, before beep
%         peakVels'<0 & indStimRpostBeep, ...% right move, right stim, after beep
%         peakVels'<0 & indStimRpreBeep, ...% right move, right stim, before beep
%         peakVels'<0 & indnoStimpostBeep, ...% right move, no stim, after beep
%         peakVels'<0 & indOthers}} ...% right move, outside trial
%         };
    conditions = {{ ... % first column (stim)
        contrastIndsL, ... % different contrasts on left
        contrastIndsR}, ... % different contrasts on right
        { ... % second column (beep)
        {peakVels'>0 & indStimLpostBeep}, ...% left move, left stim, after beep
        {peakVels'<0 & indStimRpostBeep}}, ...% right move, right stim, after beep
        { ... % third colum (choice)
        {choice==1 & outcome, ...% time of choice to left, correct choice
        choice==1 & ~outcome}, ...% time of choice to left, incorrect choice
        {choice==2 & outcome, ...% time of choice to right, correct choice
        choice==2 & ~outcome}, ...% time of choice to right, incorrect choice
        {choice==3 & outcome, ...% time of choice for NoGo, correct choice
        choice==3 & ~outcome}}, ...% time of choice for NoGo, incorrect choice
        { ... % fourth column (movements)
        {peakVels'>0 & indStimLpreBeep, ...% left move, left stim, before beep
        peakVels'>0 & indStimRpostBeep, ...% left move, right stim, after beep
        peakVels'>0 & indStimRpreBeep, ...% left move, right stim, before beep
        peakVels'>0 & indnoStimpostBeep, ...% left move, no stim, after beep
        peakVels'>0 & indOthers}, ...% left move, outside trial
        {peakVels'<0 & indStimLpostBeep, ...% right move, left stim, after beep
        peakVels'<0 & indStimLpreBeep, ...% right move, left stim, before beep
        peakVels'<0 & indStimRpreBeep, ...% right move, right stim, before beep
        peakVels'<0 & indnoStimpostBeep, ...% right move, no stim, after beep
        peakVels'<0 & indOthers}} ...% right move, outside trial
        };
    % Define event classes
    % events: {{column plots {single row plot {division within plot} } }}
%     events = {{ ... % first column
%         stimOnTimes, stimOnTimes}, ...% stim onset on right, incorrect choice
%         { ... % second column
%         beepTimes}, ... % time of beep
%         { ... % third column
%         feedbackTimes, feedbackTimes, feedbackTimes}, ...% time of choice for NoGo, incorrect choice
%         { ... % fourth column
%         allOnsets}, ...% left move, outside trial
%         { ... % fifth column
%         allOnsets} ...% right move, outside trial
%         };
    events = {{ ... % first column
        stimOnTimes, stimOnTimes}, ...% stim onset on left, on right
        { ... % second column
        allOnsets, allOnsets}, ...% left move, right move
        { ... % third column
        feedbackTimes, feedbackTimes, feedbackTimes}, ...% time of feedback
        { ... % fourth column
        allOnsets, allOnsets}, ...% left move, right move
        };
%     windows = {{stimWin, stimWin}, {beepWin}, {choiceWin, choiceWin, choiceWin}, ...
%         {moveWin}, {moveWin}};
    windows = {{stimWin, stimWin}, {moveWin, moveWin}, ...
        {choiceWin, choiceWin, choiceWin}, {moveWin, moveWin}};
%     sortCrit = {{@(x)min(x(x>0)), @(x)min(x(x>0))}, ...
%         {@(x)min(x(x>0))}, ...
%         {@(x)max(x(x<0)), @(x)max(x(x<0)), @(x)max(x(x<0))}, ...
%         {@(x)min(x(x>0))}, {@(x)min(x(x>0))}};
    sortCrit = {{@(x)min(x(x>0)), @(x)min(x(x>0))}, ...
        {@(x)min(x(x>0)), @(x)min(x(x>0))}, ...
        {@(x)max(x(x<0)), @(x)max(x(x<0)), @(x)max(x(x<0))}, ...
        {@(x)min(x(x>0)), @(x)min(x(x>0))}};
%     titles = {{'Stimulus onset - left','Stimulus onset - right'}, ...
%         {'Beep'}, ...
%         {'Feedback - left choice','Feedback - right choice','Feedback - NoGo choice'}, ...
%         {'Movement - left'}, {'Movement - right'}};
    titles = {{'Stimulus onset - left','Stimulus onset - right'}, ...
        {'Left stimulus moves',' Right stimulus moves'}, ...
        {'Feedback - left choice','Feedback - right choice','Feedback - NoGo choice'}, ...
        {'Movement - left', 'Movement - right'}};
%     ylabels = {{cellfun(@num2str,num2cell(contrastL),'UniformOutput',false), ...
%         cellfun(@num2str,num2cell(contrastR),'UniformOutput',false)}, ...
%         {{'stim L','stim R','no stim'}}, ...
%         {{'correct','incorrect'}, {'correct','incorrect'}, {'correct','incorrect'}}, ...
%         {{'stim L\newlinepost beep','stim L\newlinepre beep','stim R\newlinepost beep', ...
%         'stim R\newlinepre beep', 'no stim\newlinepost beep','others'}}, ...
%         {{'stim L\newlinepost beep','stim L\newlinepre beep','stim R\newlinepost beep', ...
%         'stim R\newlinepre beep', 'no stim\newlinepost beep','others'}}};
    ylabels = {{cellfun(@num2str,num2cell(contrastL),'UniformOutput',false), ...
        cellfun(@num2str,num2cell(contrastR),'UniformOutput',false)}, ...
        {{''}, {''}}, ...
        {{'correct','incorrect'}, {'correct','incorrect'}, {'correct','incorrect'}}, ...
        {{'stim L\newlinepre beep','stim R\newlinepost beep', ...
        'stim R\newlinepre beep', 'no stim\newlinepost beep','others'}, ...
        {'stim L\newlinepost beep','stim L\newlinepre beep', ...
        'stim R\newlinepre beep', 'no stim\newlinepost beep','others'}}};

    % Get relative event times
    choiceTimes = cell(1, length(events));
    stimTimes = cell(1, length(events));
    moveOnTimes = cell(1, length(events));
    moveOffTimes = cell(1, length(events));
%     plotStims = {[0, 0], 1, [1, 1, 1], 1, 1};
%     plotChoices = {[1, 1], 1, [0, 0, 0], 1, 1};
%     plotOnMoves = {[1, 1], 1, [1, 1, 1], 1, 1};
%     plotOffMoves = {[1, 1], 1, [1, 1, 1], 1, 1};
    plotStims = {[0 0], [0 0], [0 0 0], [0 0]};
    plotChoices = {[0 0], [0 0], [0 0 0], [0 0]};
    plotOnMoves = {[0 0], [0 0], [0 0 0], [0 0]};
    plotOffMoves = {[0 0], [0 0], [0 0 0], [0 0]};
    for c = 1:length(events)
        choiceTimes{c} = cell(1, length(events{c}));
        stimTimes{c} = cell(1, length(events{c}));
        moveOnTimes{c} = cell(1, length(events{c}));
        moveOffTimes{c} = cell(1, length(events{c}));
        for r = 1:length(events{c})
            if plotChoices{c}(r) == 1
                choiceTimes{c}{r} = cell(1, length(conditions{c}{r}));
                for d = 1:length(conditions{c}{r})
                    t = events{c}{r}(conditions{c}{r}{d});
                    ev = repmat(feedbackTimes', 1, length(t));
                    ev = ev - t;
                    ev(ev < windows{c}{r}(1) | ev > windows{c}{r}(2)) = NaN;
                    n = sum(~isnan(ev), 1);
                    choiceTimes{c}{r}{d} = mat2cell(ev(~isnan(ev)), n);
                end
            end
            if plotStims{c}(r) == 1
                stimTimes{c}{r} = cell(1, length(conditions{c}{r}));
                for d = 1:length(conditions{c}{r})
                    t = events{c}{r}(conditions{c}{r}{d});
                    ev = repmat(stimOnTimes', 1, length(t));
                    ev = ev - t;
                    ev(ev < windows{c}{r}(1) | ev > windows{c}{r}(2)) = NaN;
                    n = sum(~isnan(ev), 1);
                    stimTimes{c}{r}{d} = mat2cell(ev(~isnan(ev)), n);
                end
            end
            moveOnTimes{c}{r} = cell(1, length(conditions{c}{r}));
            for d = 1:length(conditions{c}{r})
                t = events{c}{r}(conditions{c}{r}{d});
                ev = repmat(allOnsets', 1, length(t));
                ev = ev - t;
                ev(ev < windows{c}{r}(1) | ev > windows{c}{r}(2)) = NaN;
                n = sum(~isnan(ev), 1);
                moveOnTimes{c}{r}{d} = mat2cell(ev(~isnan(ev)), n);
            end
            moveOffTimes{c}{r} = cell(1, length(conditions{c}{r}));
            for d = 1:length(conditions{c}{r})
                t = events{c}{r}(conditions{c}{r}{d});
                ev = repmat(allMoves(2,:)', 1, length(t));
                ev = ev - t;
                ev(ev < windows{c}{r}(1) | ev > windows{c}{r}(2)) = NaN;
                n = sum(~isnan(ev), 1);
                moveOffTimes{c}{r}{d} = mat2cell(ev(~isnan(ev)), n);
            end
        end
    end
    
    % Sort trials
    trialOrders = cell(1, length(events));
%     sortVars = {moveOnTimes, moveOnTimes, moveOnTimes, moveOffTimes, moveOffTimes};
    sortVars = {moveOnTimes, moveOffTimes, moveOnTimes, moveOffTimes};
    for c = 1:length(events)
        trialOrders{c} = cell(1, length(events{c}));
        if isempty(sortVars{c})
            continue
        end
        for r = 1:length(events{c})
            trialOrders{c}{r} = cell(1, length(conditions{c}{r}));
            for d = 1:length(conditions{c}{r})
                relMoveTimes = NaN(length(sortVars{c}{c}{r}{d}), 1);
                for k = 1:length(sortVars{c}{c}{r}{d})
                    if isempty(sortCrit{c}{r}(sortVars{c}{c}{r}{d}{k}))
                        continue
                    end
                    relMoveTimes(k) = sortCrit{c}{r}(sortVars{c}{c}{r}{d}{k});
                end
                [~,trialOrders{c}{r}{d}] = sort(relMoveTimes);
            end
        end
    end
    
    % Get aligned traces
    eventTimes = {};
    wins = {};
    timeBin = median(diff(caTime));
    k = 1;
    for c = 1:length(events)
        for r = 1:length(events{c})
            for d = 1:length(conditions{c}{r})
                wins{k} = windows{c}{r};
                ev = events{c}{r}(conditions{c}{r}{d});
                ind = find(ev+wins{k}(1) < caTime(1) | ...
                    ev+wins{k}(2) > caTime(end) | [false, diff(ev) < 1/sr]);
                ev(ind) = [];
                trialOrders{c}{r}{d}(any(trialOrders{c}{r}{d}==ind,2)) = [];
                num = sum(ind<length(ev)/2);
                trialOrders{c}{r}{d} = trialOrders{c}{r}{d}-num;
                
                longTimeBins = find(diff(caTime) > 1.5 * median(diff(caTime))); % between experiments
                for b = 1:length(longTimeBins)
                    ind = find(ev > caTime(longTimeBins(b)) & ...
                        ev < caTime(longTimeBins(b)+1));
                    ev(ind) = [];
                    trialOrders{c}{r}{d}(any(trialOrders{c}{r}{d}==ind,2)) = [];
                    j = find(trialOrders{c}{r}{d} > max(ind));
                    trialOrders{c}{r}{d}(j) = trialOrders{c}{r}{d}(j) - length(ind);
                end
                eventTimes{k} = ev;
                k = k + 1;
            end
        end
    end
    [A, numSamples, windowTimes] = ...
        krnl.getToeplitz(caTime, eventTimes, wins);
    [~, traces] = krnl.getETA(caTraces, A, numSamples);
    
    % Plot for each cell
    folder = fullfile(traceFolder, sprintf('%s_%s_%d', db(iSet).subject, ...
        db(iSet).date, db(iSet).exp(1)));
    if ~isdir(folder)
        mkdir(folder)
    end
    
    l = lines(5);
%     colors = 'ybmrr';
%     colors = l([1:4,4,5],:);
    colors = zeros(6,3);
    cols = length(events);
    traceMeans = nanmean(caTraces, 1);
    traceSTDs = nanstd(caTraces, 0, 1);
    for iCell = 1:size(caTraces,2)
        figure('Position', [1 41 1920 1083])
        colormap(cm)
        allTraces = [];
        for ev = 1:length(events)
            allTraces = [allTraces; reshape(traces{ev}(:,iCell,:),[],1)];
        end
        allTraces = (allTraces - traceMeans(iCell)) ./ traceSTDs(iCell);
        mini = min(allTraces);
        maxi = max(allTraces);
        maxi = max(abs([mini maxi]));
        k = 0;
        for c = 1:cols
            rows = length(events{c});
            for r = 1:rows
                subplot(rows, cols, (r-1)*cols + c)
                img = [];
                divs = zeros(1, length(conditions{c}{r}));
                valid = cell(1, length(conditions{c}{r}));
                for d = 1:length(conditions{c}{r})
                    k = k+1;
                    if isempty(trialOrders{c}{r})
                        order = (1:size(traces{k},3))';
                    else
                        order = trialOrders{c}{r}{d};
                    end
                    tr = squeeze(traces{k}(:,iCell,order))';
                    valid{d} = ~any(isnan(tr),2);
                    tr(~valid{d},:) = [];
                    img = [img; tr];
                    trialOrders{c}{r}{d} = order(valid{d});
                    divs(d) = sum(valid{d});
                end
                if isempty(img)
                    continue
                end
                xLimPos = windowTimes{k}([1 end]);
                img = (img - traceMeans(iCell)) ./ traceSTDs(iCell);
                imagesc([xLimPos(1)+.5*timeBin, xLimPos(end)-.5*timeBin], ...
                    [1 size(img,1)], img, [-maxi maxi])
                hold on
                j = 0;
                for d = 1:length(divs)
                    if d < length(divs)
                        plot(windowTimes{k}([1 end]), [1 1] .* (j+divs(d)+.5), ...
                            ':', 'Color', 'k', 'LineWidth', 2)
                    end
                    if isempty(trialOrders{c}{r})
                        order = (1:size(traces{k},3))';
                    else
                        order = trialOrders{c}{r}{d};
                    end
                    for tr = 1:divs(d)
                        if plotOnMoves{c}(r) == 1
                            t = moveOnTimes{c}{r}{d}{order(tr)};
                            plot(t, ones(1,length(t)).*(tr+j), '.', ...
                                'Color', colors(4,:))
                        end
                        if plotOffMoves{c}(r) == 1
                            t = moveOffTimes{c}{r}{d}{order(tr)};
                            plot(t, ones(1,length(t)).*(tr+j), '.', ...
                                'Color', colors(6,:))
                        end
                        if plotStims{c}(r) == 1
                            t = stimTimes{c}{r}{d}{order(tr)};
                            if ~isempty(t)
                                plot(t, tr+j,'>','Color',colors(1,:),'MarkerSize',4)
                            end
                        end
                        if plotChoices{c}(r) == 1
                            t = choiceTimes{c}{r}{d}{order(tr)};
                            if ~isempty(t)
                                plot(t,tr+j,'<','Color',colors(3,:),'MarkerSize',4)
                            end
                        end
                    end
                    j = j + divs(d);
                end
                plot([0 0],[.5 size(img,1)+.5],':','Color',colors(c,:),'LineWidth',2)
                title(titles{c}{r})
                xlabel('Time (s)')
                ticks = divs ./ 2 + cumsum([1 divs(1:end-1)]);
                set(gca, 'YTick', ticks, 'YTickLabel', ylabels{c}{r}, ...
                    'YTickLabelRotation', 0)
                %             ylabel('# Events')
            end
        end
        c = colorbar;
        c.Label.String = 'Normalised response';
        c.Position = [.937 .11 .014 .342];
        
        f = gcf;
        f.PaperPositionMode = 'auto';
        print(fullfile(folder, sprintf('eventTraces_plane%d_cell%d.tiff', ...
            cellIDs(iCell,1), cellIDs(iCell,2))), '-dtiff','-r0')
        close gcf
    end
    
    %% Define events, windows and continuous vectors
    % kernels: (1) choice left + correct, (2) choice left + incorrect, (3)
    %          choice right + correct, (4) choice right + incorrect, (5) choice
    %          nogo + correct, (6) choice nogo + incorrect,
    %          (7) stim left home, (8) stim right home,
    %          (9) stim left moves to front, (10) stim left moves to back,
    %          (11) stim right moves to front, (12) stim right moves to back,
    %          (13) wheel moves
    % ETAs:    (1-C) onset of stim with certain contrast (left or right),
    %          (C+1) stim left moves to front, (C+2) stim left moves to back,
    %          (C+3) stim right moves to front, (C+4) stim right moves to back,
    %          (C+5) choice left + correct, (C+6) choice left + incorrect,
    %          (C+7) choice right + correct, (C+8) choice right + incorrect,
    %          (C+9) choice nogo + correct, (C+10) choice nogo + incorrect,
    %          (C+11) onset of wheel movement

    % vectors for stimulus at home position and for moving stimulus
    stimLHome = zeros(length(caTime),1);
    stimRHome = zeros(length(caTime),1);
    stimLMovesToFront = zeros(length(caTime),1);
    stimLMovesToBack = zeros(length(caTime),1);
    stimRMovesToFront = zeros(length(caTime),1);
    stimRMovesToBack = zeros(length(caTime),1);
    indL = find(stimulus(:,1)>0)';
    indR = find(stimulus(:,2)>0)';
    for k = indL
        contrast = stimulus(k,1);
        t1 = find(caTime >= stimOnTimes(k), 1); % stim onset
        onsets = find(allOnsets>beepTimes(k) & allOnsets<feedbackTimes(k)); % onset of movements between beep and time of choice
        if isempty(onsets)
            t2 = find(caTime <= feedbackTimes(k), 1, 'last');
        else
            t2 = find(caTime < allOnsets(onsets(1)), 1, 'last'); % first movement after beep (stimulus moves away from home position)
        end
        stimLHome(t1 : t2) = contrast;
        
        ind = peakVels(onsets) > 0; % movements to the front (choosing left)
        ons = onsets(ind);
        for j = ons
            t1 = find(caTime >= allOnsets(j), 1); % movement onset
            t2 = find(caTime <= min(allMoves(2,j),feedbackTimes(k)), 1, 'last'); % movement offset or time of choice
            stimLMovesToFront(t1 : t2) = contrast;
        end
        ons = onsets(~ind); % movements to the back (choosing right)
        for j = ons
            t1 = find(caTime >= allOnsets(j), 1); % movement onset
            t2 = find(caTime <= min(allMoves(2,j),feedbackTimes(k)), 1, 'last'); % movement offset or time of choice
            stimLMovesToBack(t1 : t2) = contrast;
        end
        
        offset = allMoves(2,:) > beepTimes(k) & allOnsets < beepTimes(k);
        if ~isempty(offset)
            t1 = find(caTime >= beepTimes(k), 1);
            t2 = find(caTime <= allMoves(2,offset), 1, 'last');
            if peakVels(offset) > 0 % movements to the front (choosing left)
                stimLMovesToFront(t1 : t2) = contrast;
            else % movements to the back (choosing right)
                stimLMovesToBack(t1 : t2) = contrast;
            end
        end
    end
    for k = indR
        contrast = stimulus(k,2);
        t1 = find(caTime >= stimOnTimes(k), 1); % stim onset
        onsets = find(allOnsets>beepTimes(k) & allOnsets<feedbackTimes(k)); % onset of movements between beep and time of choice
        if isempty(onsets)
            t2 = find(caTime <= feedbackTimes(k), 1, 'last');
        else
            t2 = find(caTime <= allOnsets(onsets(1)), 1, 'last'); % first movement after beep (stimulus moves away from home position)
        end
        stimRHome(t1 : t2) = contrast;
        
        ind = peakVels(onsets) < 0; % movements to the front (choosing right)
        ons = onsets(ind);
        for j = ons
            t1 = find(caTime >= allOnsets(j), 1);
            t2 = find(caTime <= min(allMoves(2,j),feedbackTimes(k)), 1, 'last');
            stimRMovesToFront(t1 : t2) = contrast;
        end
        ons = onsets(~ind); % movements to the back (choosing left)
        for j = ons
            t1 = find(caTime >= allOnsets(j), 1); % movement onset
            t2 = find(caTime <= min(allMoves(2,j),feedbackTimes(k)), 1, 'last'); % movement offset or time of choice
            stimRMovesToBack(t1 : t2) = contrast;
        end
        
        offset = allMoves(2,:) > beepTimes(k) & allOnsets < beepTimes(k);
        if ~isempty(offset)
            t1 = find(caTime >= beepTimes(k), 1);
            t2 = find(caTime <= allMoves(2,offset), 1, 'last');
            if peakVels(offset) < 0 % movements to the front (choosing right)
                stimRMovesToFront(t1 : t2) = contrast;
            else % movements to the back (choosing right)
                stimRMovesToBack(t1 : t2) = contrast;
            end
        end
    end
    % vector for wheel movements
    moveL = zeros(length(caTime),1);
    moveR = zeros(length(caTime),1);
    for k = 1:length(allOnsets)
        t1 = find(caTime >= allOnsets(k), 1);
        t2 = find(caTime <= allMoves(2,k), 1, 'last');
        if peakVels(k) > 0 % choosing left
            moveL(t1 : t2) = 1;
        else
            moveR(t1 : t2) = 1;
        end
    end
    
    % vectors for successive steps of fitting
    % only stim in home position
    v = {stimLHome, stimRHome};
    vectors1 = {cellfun(@double, cellfun(@gt, v, num2cell(zeros(size(v))), ...
        'UniformOutput', false), 'UniformOutput', false), ... % stim present
        v, ...                                                % stim contrast
        cellfun(@sqrt, v, 'UniformOutput', false)};    % sqrt contrast
    vector1Windows_krnl = {stimWin_krnl, stimWin_krnl};
    % only moving stim
    v = {stimLMovesToFront, stimLMovesToBack, stimRMovesToFront, stimRMovesToBack};
    vectors2 = {cellfun(@double, cellfun(@gt, v, num2cell(zeros(size(v))), ...
        'UniformOutput', false), 'UniformOutput', false), ... % stim present
        v, ...                                                % stim contrast
        cellfun(@sqrt, v, 'UniformOutput', false)};    % sqrt contrast
    vector2Windows_krnl = {moveStimWin_krnl, moveStimWin_krnl, ...
        moveStimWin_krnl, moveStimWin_krnl};
    vectors3 = {moveL, moveR}; % only non-visual events
    vector3Windows_krnl = {moveWin_krnl, moveWin_krnl};
    
    % events for successive steps of fitting
    events1 = {};
    event1Windows_krnl = {};
    events2 = {};
    event2Windows_krnl = {};
    events3 = {feedbackTimes(choice==1 & outcome), ...                      % time of choice to left
        feedbackTimes(choice==1 & ~outcome), ...
        feedbackTimes(choice==2 & outcome), ...                            % time of choice to right
        feedbackTimes(choice==2 & ~outcome), ...
        feedbackTimes(choice==3 & outcome), ...                            % time of choice to hold still
        feedbackTimes(choice==3 & ~outcome)};
    event3Windows_krnl = {rewardWin_krnl, rewardWin_krnl, rewardWin_krnl, ...
        rewardWin_krnl, rewardWin_krnl, rewardWin_krnl};
    
    % events for stimulus onset and onset of stim moving fro ETAs
    contrasts = setdiff(unique(diff(stimulus,1,2)), 0);
    evC = cell(1, length(contrasts));
    for c = 1:length(contrasts)
        if contrasts(c) < 0
            evC{c} = stimOnTimes(stimulus(:,1)==-contrasts(c));
        else
            evC{c} = stimOnTimes(stimulus(:,2)==contrasts(c));
        end
    end
    evLToFront = [];
    evLToBack = [];
    evRToFront = [];
    evRToBack = [];
    indL = find(stimulus(:,1)>0)';
    indR = find(stimulus(:,2)>0)';
    for k = indL
        ind = allOnsets>beepTimes(k) & allOnsets<feedbackTimes(k);
        evLToFront = [evLToFront; allOnsets(ind & peakVels>0)'];
        evLToBack = [evLToBack; allOnsets(ind & peakVels<0)'];
    end
    for k = indR
        ind = allOnsets>beepTimes(k) & allOnsets<feedbackTimes(k);
        evRToFront = [evRToFront; allOnsets(ind & peakVels<0)'];
        evRToBack = [evRToBack; allOnsets(ind & peakVels>0)'];
    end
    eventsETA = [evC, {evLToFront, evLToBack, evRToFront, evRToBack}, ...
        events3, {allOnsets(peakVels>0), allOnsets(peakVels<0)}];
    eventWindowsETA = [repmat({stimWin_eta},1,length(contrasts)), ...
        repmat({moveStimWin_eta}, 1, 4), ...
        repmat({rewardWin_eta}, 1, 6), ...
        repmat({moveWin_eta}, 1, 2)];


    %% Select model with crossvalidation
    partition = 20; % k-fold crossvalidation
    ev = cell(1, 3);
    ev{1} = NaN(size(caTraces,2), 4); % stim home: (1) none, (2) pres, (3) contr, (4) sqrt contr
    ev{2} = NaN(size(caTraces,2), 4); % stim moves: (1) none, (2) pres, (3) contr, (4) sqrt contr
    ev{3} = NaN(size(caTraces,2), 2); % non-vis: (1) none, (2) pres
    A = krnl.getToeplitz(caTime, [events1, events2, events3], ...
        [event1Windows_krnl, event2Windows_krnl, event3Windows_krnl], ...
        [vectors1{1}, vectors2{1}, vectors3{1}], ...
        [vector1Windows_krnl, vector2Windows_krnl, vector3Windows_krnl]);
    tValid = sum(A,2) > 0;
    fprintf('  Crossvalidation, cell (of %d):', size(caTraces,2))
    for iCell = 1:size(caTraces,2)
        fprintf(' %d', iCell)
        
        residuals1 = cell(1, 4);
        residuals1{1} = caTraces(:,iCell);
        % (0) No kernels, just mean response
        residuals = krnl.crossvalidate(caTraces(:,iCell), caTime, partition, ...
            {}, {}, {}, {}, lambda);
        ind = tValid & ~isnan(residuals) & ~isnan(caTraces(:,iCell));
        ev{1}(iCell,1) = 1 - sum(residuals(ind).^2) / ...
            sum((caTraces(ind,iCell)-mean(caTraces(ind,iCell))).^2);
        
        % (1) Stimulus in home position
        for k = 1:3
            residuals1{k+1} = krnl.crossvalidate(caTraces(:,iCell), caTime, partition, ...
                events1, event1Windows_krnl, vectors1{k}, vector1Windows_krnl, lambda);
            ind = tValid & ~isnan(residuals1{k+1}) & ~isnan(caTraces(:,iCell));
            ev{1}(iCell,k+1) = 1 - sum(residuals1{k+1}(ind).^2) / ...
                sum((caTraces(ind,iCell)-mean(caTraces(ind,iCell))).^2);
        end
        % select best stim model
        [~,bestStim] = max(ev{1}(iCell,:));
        
        % (2) Stimulus moving
        residuals2 = cell(1, 4);
        sig = residuals1{bestStim};
        % (0) no kernels
        residuals2{1} = sig;
        ev{2}(iCell,1) = ev{1}(iCell, bestStim);
        for k = 1:3
            residuals2{1+k} = krnl.crossvalidate(sig, caTime, partition, ...
                events2, event2Windows_krnl, vectors2{k}, vector2Windows_krnl, lambda);
            ind = tValid & ~isnan(residuals2{1+k}) & ~isnan(sig);
            ev{2}(iCell,1+k) = 1 - sum(residuals2{1+k}(ind).^2) / ...
                sum((caTraces(ind,iCell)-mean(caTraces(ind,iCell))).^2);
        end
        % select best stim model
        [~,bestStim] = max(ev{2}(iCell,:));
        
        % (3) Non-visually related kernels
        sig = residuals2{bestStim};
        % (0) no kernels
        ev{3}(iCell,1) = ev{2}(iCell, bestStim);
        % (a) kernels present
        residuals = krnl.crossvalidate(sig, caTime, partition, ...
            events3, event3Windows_krnl, vectors3, vector3Windows_krnl, lambda);
        ind = tValid & ~isnan(residuals) & ~isnan(sig);
        ev{3}(iCell,2) = 1 - sum(residuals(ind).^2) / ...
            sum((caTraces(ind,iCell)-mean(caTraces(ind,iCell))).^2);
    end
    fprintf('\n')
    
    results(iSet).expVar = ev;


    %% Perform linear regression and get event triggered averages
    % (1) get event-triggered averages
    [A, numSamples, windowTimes_eta] = ...
        krnl.getToeplitz(caTime, eventsETA, eventWindowsETA);
    etas = krnl.getETA(caTraces, A, numSamples);
    
    % (2) get kernels
    [maxEV, bestModel] = cellfun(@max, ev, {[],[],[]}, {2,2,2}, 'UniformOutput', false);
    maxEV = cat(2, maxEV{:});
    bestModel = cat(2, bestModel{:});
    events = { { cell(size(events1)), events1, events1, events1 }, ...
        { cell(size(events2)), events2, events2, events2 }, ...
        { cell(size(events3)), events3, events3, events3 } };
    evWins = [event1Windows_krnl, event2Windows_krnl, event3Windows_krnl];
    vectors = { [ {cell(size(vectors1{1}))} , vectors1 ], ...
        [ {cell(size(vectors2{1}))}, vectors2 ], ...
        [ {cell(size(vectors3))}, repmat({vectors3},1,3) ] };
    vecWins = [vector1Windows_krnl, vector2Windows_krnl, vector3Windows_krnl];
    [~, numSamples] = krnl.getToeplitz(caTime, ...
        [events{1}{end}, events{2}{end}, events{3}{end}], evWins, ...
        [vectors{1}{end}, vectors{2}{end}, vectors{3}{end}], vecWins);
    
    fitKernels = cell(1, length(numSamples));
    predictions = NaN(size(caTraces));
    for k = 1:length(fitKernels)
        fitKernels{k} = NaN(numSamples(k), size(caTraces,2));
    end
    fprintf('  Fitting models (of %d):', size(ev{1},2)*size(ev{2},2)*size(ev{3},2))
    j = 0;
    for b1 = 1:size(ev{1},2)
        for b2 = 1:size(ev{2},2)
            for b3 = 1:size(ev{3},2)
                j = j + 1;
                fprintf(' %d', j)
                if b1==1 && b2==1 && b3==1 % no kernels is best
                    continue
                end
                ind = all(bestModel == [b1 b2 b3], 2);
                [krnls, ~, predSig, ~, windowTimes_krnl] = ...
                    krnl.kernelRegression(caTraces(:,ind), caTime, ...
                    [events{1}{b1}, events{2}{b2}, events{3}{b3}], evWins, ...
                    [vectors{1}{b1}, vectors{2}{b2}, vectors{3}{b3}], vecWins, ...
                    lambda, true);
                for k = 1:length(fitKernels)
                    fitKernels{k}(:,ind) = krnls{k};
                end
                predictions(:,ind) = predSig;
            end
        end
    end
    fprintf('\n')

    % (3) get event-triggered averages from predicted signal
    [A, numSamples] = krnl.getToeplitz(caTime, eventsETA, eventWindowsETA);
    predETAs = krnl.getETA(predictions, A, numSamples);
    
    results(iSet).fitKernels = fitKernels;
    results(iSet).windowTimes_krnl = windowTimes_krnl;
    results(iSet).bestModel = bestModel;
    results(iSet).predictions = predictions;
    results(iSet).predETAs = predETAs;
    results(iSet).etas = etas;
    results(iSet).windowTimes_eta = windowTimes_eta;
    results(iSet).contrasts = contrasts;
    
    %% Save results
    save(fullfile(resultsFolder, 'kernels.mat'), 'results', 'kernelNames', ...
        'modelNames');

end

%% Plot kernels and event triggered averages for each neuron
% load kernel fit results
data = load(fullfile(resultsFolder, 'kernels.mat'));
results = data.results;
modelNames = data.modelNames;

% define plotting
subplotsKrnl = {[7 8], 9:12, [1 3 5], [2 4 6], 13:14}; % what kernels to draw in each plot
% subplotsKrnl = {[7 8], [9 11], [1 3 5 6], 13:14};
subplotETAsDiff = {-12, -11:-8, [-7 -5 -3], [-6 -4 -2], -1:0};
swapHemispheres = {[], [3 4 1 2], [2 1 3], [2 1 3], [2 1]};
titles = {'Stimulus onset', 'Stimulus moves', 'Reward', 'Penalty', 'Wheel moves'};
labels = {{'contra front','contra back','ipsi front','ipsi back'}, ...
        {'contra','ipsi','nogo'}, {'contra','ipsi','nogo'}, {'contra','ipsi'}};
% labels = {{'left front','right front'}, ...
%     {'L 1', 'R 1', 'nogo 1', 'nogo 0'}, {'left','right'}};
cols = lines(5);
colors = {[1 0 0;1 0 0;0 0 1;0 0 1], ...
    [1 0 0;0 0 1;cols(3,:)], [1 0 0;0 0 1;cols(3,:)], cols(4:5,:)};
% colors = {[1 0 0;0 0 1], [1 0 0;0 0 1;cols(3,:);cols(3,:)], cols(4:5,:)};
lins = {{'-',':','-',':'}, {'-','-','-'}, {'-','-','-'}, {'-','-'}};
% lins = {{'-','-'}, {'-','-','-',':'}, {'-','-'}};


for iSet = 1:length(results)
    cellIDs = results(iSet).cellIDs;
    isGad = results(iSet).isGad;
    expVar = results(iSet).expVar;
    fitKernels = results(iSet).fitKernels;
    windowTimes_krnl = results(iSet).windowTimes_krnl;
    bestModel = results(iSet).bestModel;
    predictions = results(iSet).predictions;
    predETAs = results(iSet).predETAs;
    etas = results(iSet).etas;
    windowTimes_eta = results(iSet).windowTimes_eta;
    contrasts = results(iSet).contrasts;
    % define plotting
    k = length(etas);
%     subplotsETA = {1:k-12, k-11:k-8, k-7:k-2, k-1:k};
    subplotsETA = [{1:k+subplotETAsDiff{1}}, ...
        cellfun(@plus, num2cell(ones(1, length(subplotETAsDiff)-1).*k), ...
        subplotETAsDiff(2:end), 'UniformOutput', false)];
    spK = subplotsKrnl;
    if strcmp(results(iSet).hemisphere, 'left')
        subplotsETA{1} = flip(subplotsETA{1});
        spK{1} = flip(spK{1});
        for k = 2:length(subplotsETA)
            subplotsETA{k} = subplotsETA{k}(swapHemispheres{k});
            spK{k} = spK{k}(swapHemispheres{k});
        end
    end
    labs = [{num2str(contrasts(:))}, labels];
    
    n = sum(contrasts<0);
    w = linspace(1, .2, n-1)';
    blues = [0 0 .8; [0 0 1] .* w + (1-w)];
    n = sum(contrasts>0);
    w = linspace(1, .2, n-1)';
    reds = [.8 0 0; [1 0 0] .* w + (1-w)];
    colorsKrnl = [{[1 0 0; 0 0 1]}, colors];
    colorsETA = [{[reds;flip(blues,1)]}, colors];
    linesETA = [{repmat({'-'},length(contrasts),1)}, lins];

    % Plot results for each cell
    folder = fullfile(resultsFolder, sprintf('%s_%s_%d', results(iSet).subject, ...
        results(iSet).date, results(iSet).exp(1)));
    if ~isdir(folder)
        mkdir(folder)
    end
    
    limits = cell(size(subplotsETA));
    for k = 1:length(limits)
        limits{k} = [min(windowTimes_eta{subplotsETA{k}(1)}(1), ...
            windowTimes_krnl{spK{k}(1)}(1)), ...
            max(windowTimes_eta{subplotsETA{k}(1)}(end), ...
            windowTimes_krnl{spK{k}(1)}(end))];
    end
    allKernels = cat(1, fitKernels{:});
    allETAs = cat(1, etas{:});
    allPredETAs = cat(1, predETAs{:});

    for iCell = 1:size(results(iSet).cellIDs,1)
        figure('Position', [1 41 1920 1083])
        
        annotation('textbox', [0 0.96 1 0.03], 'String', ...
            sprintf('Best model: %s | %s | %s, expl. var.: %.4f', ...
            modelNames{1}{bestModel(iCell,1)}, modelNames{2}{bestModel(iCell,2)}, ...
            modelNames{3}{bestModel(iCell,3)}, expVar{3}(iCell, bestModel(iCell,3))), ...
            'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center')
        
        maxi = max([allPredETAs(:,iCell); allETAs(:,iCell)]);
        mini = min([allPredETAs(:,iCell); allETAs(:,iCell)]);
        rng = maxi - mini;
        maxi = maxi + .05*rng;
        mini = mini - .05*rng;
        if rng == 0
            maxi = 1;
            mini = -1;
        end
        for sp = 1:length(subplotsETA)
            subplot(3,length(subplotsETA),sp)
            hold on
            plot([0 0], [mini maxi], 'k:')
            h = zeros(1,length(subplotsETA{sp}));
            for l = 1:length(subplotsETA{sp})
                h(l) = plot(windowTimes_eta{subplotsETA{sp}(l)}, ...
                    etas{subplotsETA{sp}(l)}(:,iCell), linesETA{sp}{l}, ...
                    'Color', colorsETA{sp}(l,:), 'Linewidth', 2);
            end
            if ~isempty(labs{sp})
                leg = legend(h, labs{sp});
                if sp == 1
                    leg.Position = [.058 .709 .035 .216];
                else
                    leg.Color = 'none';
                    leg.Location = 'NorthWest';
                end
            end
            if sp == 1
                ylabel('Event triggered average')
            end
            set(gca, 'XTick', ceil(limits{sp}(1)) : floor(limits{sp}(end)))
            xlim(limits{sp})
            ylim([mini maxi])
            title(titles{sp})
        end
        
        if all(bestModel(iCell,:) == 1)
            f = gcf;
            f.PaperPositionMode = 'auto';
            print(fullfile(folder, sprintf('kernelFit_plane%d_cell%d', ...
                cellIDs(iCell,1), cellIDs(iCell,2))), '-djpeg', '-r0')
            close gcf
            continue
        end
        
        maxi = max([allPredETAs(:,iCell); allETAs(:,iCell)]);
        mini = min([allPredETAs(:,iCell); allETAs(:,iCell)]);
        rng = maxi - mini;
        maxi = maxi + .05*rng;
        mini = mini - .05*rng;
        if rng == 0
            maxi = 1;
            mini = -1;
        end
        for sp = 1:length(subplotsETA)
            subplot(3, length(subplotsETA), length(subplotsETA)+sp)
            hold on
            plot([0 0], [mini maxi], 'k:')
            for l = 1:length(subplotsETA{sp})
                plot(windowTimes_eta{subplotsETA{sp}(l)}, predETAs{subplotsETA{sp}(l)}(:,iCell), ...
                    linesETA{sp}{l}, 'Color', colorsETA{sp}(l,:), 'LineWidth', 2);
                %             plot(windowTimes_eta{subplots{sp}(l)}, etas{subplots{sp}(l)}(:,iCell), ...
                %                 ':', 'Color', colors{sp}(l,:), 'LineWidth', 2);
            end
            if sp == 1
                ylabel('Predicted event triggered average')
            end
            set(gca, 'XTick', ceil(limits{sp}(1)) : floor(limits{sp}(end)))
            xlim(limits{sp})
            ylim([mini maxi])
            xlabel('Time from event onset')
        end
        
        maxi = max(allKernels(:,iCell));
        mini = min(allKernels(:,iCell));
        rng = maxi - mini;
        maxi = maxi + .05*rng;
        mini = mini - .05*rng;
        if rng == 0
            maxi = 1;
            mini = -1;
        end
        for sp = 1:length(spK)
            subplot(3, length(spK), 2*length(subplotsETA) + sp)
            hold on
            plot([0 0], [mini maxi], 'k:')
            for l = 1:length(spK{sp})
                plot(windowTimes_krnl{spK{sp}(l)}, ...
                    fitKernels{spK{sp}(l)}(:,iCell), linesETA{sp}{l}, ...
                    'Color', colorsKrnl{sp}(l,:), 'LineWidth', 2)
            end
            if sp == 1
                ylabel('Fitted kernel')
            end
            set(gca, 'XTick', ceil(limits{sp}(1)) : floor(limits{sp}(end)))
            xlim(limits{sp})
            ylim([mini maxi])
        end
        
        f = gcf;
        f.PaperPositionMode = 'auto';
        print(fullfile(folder, sprintf('kernelFit_plane%d_cell%d', ...
            cellIDs(iCell,1), cellIDs(iCell,2))), '-djpeg', '-r0')
        close gcf
    end
end

%% Plot kernel sizes for each dataset

% load kernel fit results
data = load(fullfile(resultsFolder, 'kernels.mat'));
results = data.results;
kernelNames = data.kernelNames;

for iSet = 1:length(results)
    expVar = results(iSet).expVar;
    fitKernels = results(iSet).fitKernels;
    bestModel = results(iSet).bestModel;

    valid = ~all(bestModel == 1, 2);

    % summarise each kernel as one number (integral of absolute values)
    integrals = NaN(length(fitKernels), length(valid));
    for k = 1:length(fitKernels)
        integrals(k,valid) = sum(abs(fitKernels{k}(:,valid))) ./ ...
            size(fitKernels{k},1);
    end
    % z-score integrals
    integrals = (integrals - mean(integrals,1)) ./ std(integrals,0,1);

    % sort neurons according to the max. kernel
    sortedNeurons = cell(1, length(fitKernels));
    [m,maxKernel] = max(integrals, [], 1);
    for k = 1:length(fitKernels)
        n = find(maxKernel==k & ~isnan(m));
        [~,s] = sort(m(n), 'descend');
        sortedNeurons{k} = n(s);
    end
    sn = cat(2,sortedNeurons{:});
    
    figure('Position', [585 395 1290 705])
    subplot(5,1,1:4)
    imagesc(integrals(:,sn))
    set(gca,'YTick',1:length(fitKernels), 'YTickLabel', kernelNames)
    c = colorbar;
    c.Label.String = 'Kernel integral (norm.)';
    c.Position = [.92 .28 .018 .642];
    ax(1) = gca;
    expStr = sprintf('%d ', results(iSet).exp);
    expStr(end) = [];
    title(sprintf('%s %s exp. %s, %s hemisphere', results(iSet).subject, ...
        results(iSet).date, expStr, results(iSet).hemisphere))
    
    subplot(5,1,5)
    ind = sub2ind(size(expVar{3}), (1:size(expVar{3},1))', bestModel(:,3));
    expV = expVar{3}(ind);
    expV = expV(sn);
    stem(expV, 'filled', 'MarkerSize', 2)
    ax(2) = gca;
    linkaxes(ax,'x')
    xlim([.5 length(sn)+.5])
    ylim([min(expV) max(expV)])
    set(gca, 'box', 'off')
    xlabel('Neurons')
    ylabel('Expl. var.')
    
    folder = fullfile(resultsFolder, sprintf('%s_%s_%d', results(iSet).subject, ...
        results(iSet).date, results(iSet).exp(1)));
    if ~isdir(folder)
        mkdir(folder)
    end
    f = gcf;
    f.PaperPositionMode = 'auto';
    print(fullfile(folder, 'allKernels'), '-djpeg', '-r0')
    close gcf
end

%% Plot kernel sizes and ETAs for all dataset together
minEV = .05;
timeBin = .1;

% load kernel fit results
data = load(fullfile(resultsFolder, 'kernels.mat'));
results = data.results;

kernelNames = {'choice contra corr', 'choice contra incorr', 'choice ipsi corr', ...
    'choice ipsi incorr', 'choice NoGo corr', 'choice NoGo incorr',...
    'stim contra home', 'stim ipsi home', ...
    'stim contra to front', 'stim contra to back', 'stim ipsi to front', 'stim ipsi to back', ...
    'move wheel contra', 'move wheel ipsi'};
% kernelOrder1 = [7:12, 1:6, 13:14];
% kernelOrder1 = [7, 13:14, 9, 11, 1, 3, 5:6];
kernelOrder1 = [7, 13:14, 9, 1, 3, 5];

etaNames = {'stim contra on', 'stim ipsi on', ...
    'stim contra moves', 'stim ipsi moves', ...
    'reward contra', 'reward ipsi', 'reward nogo', 'punish nogo', ...
    'wheel turn contra', 'wheel turn ipsi'};
% etaNames = {'stim contra on', ...
%     'wheel turn contra', 'wheel turn ipsi', ...
%     'stim contra moves', 'stim ipsi moves', ...
%     'reward contra', 'reward ipsi', 'reward nogo', 'punish nogo'};
etaOrder = [1 9 10 3 5 6 7];

integrals = [];
explainedVariance = [];
etas = cell(length(results), length(etaNames));
predETAs = cell(length(results), length(etaNames));
times_eta = cell(length(results), length(etaNames));
RFdists = zeros(0, 2); % if NaN, RFs in dataset were mapped but neuron 
                       % does not have an RF; if Inf, RFs of neurons were
                       % not mapped
RFdists_global = zeros(0,2);

leftHem = find(strcmp('left', {results.hemisphere}));
for iSet = 1:length(results)
    expVar = results(iSet).expVar;
    fitKernels = results(iSet).fitKernels;
    bestModel = results(iSet).bestModel;
    et = results(iSet).etas;
    pe = results(iSet).predETAs;
    contrasts = results(iSet).contrasts;
    cellIDs = results(iSet).cellIDs;
    
    % get RF distances to stimulus home position
    folder = dir(fullfile(RFfolder, sprintf('%s_%s_*', ...
        results(iSet).subject, results(iSet).date)));
    if isempty(folder)
        RF = [Inf Inf];
    else
        data = load(fullfile(RFfolder, folder.name, 'RFdist.mat'));
        dists = data.RFdistance;
        if iscell(dists)
            planes = unique(cellIDs(:,1));
            RF = [];
            for pl = 1:length(planes)
                ind = cellIDs(:,1) == planes(pl);
                RF = [RF; dists{pl}(cellIDs(ind,2),:)];
            end
        else
            RF = dists;
        end
    end
    
    ind = sub2ind(size(expVar{3}), (1:size(expVar{3},1))', bestModel(:,3));
    ev = expVar{3}(ind);
    valid = ~all(bestModel == 1, 2) & ev >= minEV;

    % summarise each kernel as one number (integral of absolute values)
    intgrls = NaN(length(fitKernels), sum(valid));
    for k = 1:length(fitKernels)
        intgrls(k,:) = sum(abs(fitKernels{k}(:,valid))) ./ ...
            size(fitKernels{k},1);
    end
    % z-score integrals
    intgrls = (intgrls - mean(intgrls,1)) ./ std(intgrls,0,1);
    
    % get ETAs and predicted ETAs
    c = length(contrasts);
    indETAs = [1 c c+[1 3 5 7 9 10 11 12]];
    traceMean = nanmean(results(iSet).caTraces(:,valid), 1);
    traceSTD = nanstd(results(iSet).caTraces(:,valid), 0, 1);
    for k = 1:length(indETAs)
        etas{iSet,k} = (et{indETAs(k)}(:,valid) - traceMean) ./ traceSTD;
        times_eta{iSet,k} = results(iSet).windowTimes_eta{indETAs(k)};
        predETAs{iSet,k} = (pe{indETAs(k)}(:,valid) - traceMean) ./ traceSTD;
    end
    
    % sort kernels and ETAs according to contra- and ipsilateral, instead 
    % of left and right
    if ismember(iSet, leftHem)
        intgrls = intgrls([3 4 1 2 5 6 8 7 11 12 9 10 14 13],:);
        etas(iSet,:) = etas(iSet,[2 1 4 3 6 5 7 8 9 10]);
        predETAs(iSet,:) = predETAs(iSet,[2 1 4 3 6 5 7 8 9 10]);
    end
    
    integrals = [integrals, intgrls];
    explainedVariance = [explainedVariance; ev(valid)];
    if size(RF,1) > 1
        RFdists = [RFdists; RF(valid,:)];
        RFglobal = nanmean(RF, 1);
    else
        RFdists = [RFdists; Inf(sum(valid),2)];
        RFglobal = RF;
    end
    RFdists_global = [RFdists_global; repmat(RFglobal, sum(valid), 1)];
end

% make timing of ETAs consistent across datasets
for k = 1:size(etas,2)
    mini = round(max(cellfun(@min, times_eta(:,k))) / (timeBin/2)) * (timeBin/2);
    maxi = round(min(cellfun(@max, times_eta(:,k))) / (timeBin/2)) * (timeBin/2);
    t = mini : timeBin : maxi;
    for j = 1:size(etas,1)
        if all(all(isnan(etas{j,k}),1),2)
            etas{j,k} = NaN(length(t), size(etas{j,k},2));
            predETAs{j,k} = NaN(length(t), size(predETAs{j,k},2));
        else
            etas{j,k} = interp1(times_eta{j,k}, etas{j,k}, t, 'pchip');
            predETAs{j,k} = interp1(times_eta{j,k}, predETAs{j,k}, t, 'pchip');
        end
    end
    etas{1,k} = cat(2, etas{:,k});
    predETAs{1,k} = cat(2, predETAs{:,k});
    times_eta{1,k} = t;
end
etas = etas(1,etaOrder);
predETAs = predETAs(1,etaOrder);
times_eta = times_eta(1,etaOrder);

% sort neurons according to the max. kernel
sortedNeurons = cell(1, size(integrals,1));
[m,maxKernel] = max(integrals, [], 1);
intgr = [];
expVar = [];
for k = 1:size(integrals,1)
    n = find(maxKernel==k & ~isnan(m));
    [~,s] = sort(m(n), 'descend');
    sortedNeurons{k} = n(s);
    intgr = [intgr, NaN(size(integrals,1),3), integrals(:,n(s))];
    expVar = [expVar; NaN(3,1); explainedVariance(n(s))];
end
intgr(:,1:3) = [];
expVar(1:3) = [];
rejected = setdiff(1:length(sortedNeurons), kernelOrder1);
sortedNeurons = [sortedNeurons(kernelOrder1), {cat(2,sortedNeurons{rejected})}];
sn = cat(2,sortedNeurons{:});

integrals = integrals(kernelOrder1,:);

% plot kernel sizes
figure('Position', [1 41 1920 1083])
subplot(1,6,1:5)
imagesc(integrals(:,sn)')
% imagesc(intgr')
set(gca,'XTick',1:length(fitKernels), 'XTickLabel', ...
    kernelNames(kernelOrder1), 'XTickLabelRotation', 90, 'box', 'off')
c = colorbar('Ticks', -1:3);
c.Label.String = 'Kernel integral (norm.)';
c.Position = [.073 .11 .02 .815];
ylabel('Neurons')

subplot(1,6,6)
stem(expVar, 'k', 'filled', 'MarkerSize', 2)
view(90,90)
ylim([0 max(explainedVariance)])
xlim([.5 length(expVar)+.5])
set(gca, 'box', 'off', 'XTick', [])
ylabel('Explained variance')

% Plot kernel sizes only of neurons whose RFs are within stimulus home
% position, only of neurons whose RFs are outside, and only neurons with a
% population RF within the stimulus home position but that don't a RF
% themselves (non-visual?)
maxDist = 2; % 2 sigma from center of stimulus position
minDist = 2;
dists = sqrt(sum(RFdists.^2, 2));
dists = dists(sn);
globalDists = sqrt(sum(RFdists_global.^2, 2));
globalDists = globalDists(sn);
ints = integrals(:,sn);
ev = explainedVariance(sn);
groupSizes = cellfun(@length, sortedNeurons);
allInts = ints(:,~isinf(globalDists));
mini = min(allInts(:));
maxi = max(allInts(:));
cm = flip(gray,1);

% Neuron RFs are within stimulus position
figure('Position', [1 41 1920 1083])
colormap(cm)
subplot(1,6,1:5)
valid = dists < maxDist;
imagesc(ints(:, valid)',[mini maxi])
gs = cellfun(@sum,mat2cell(valid, groupSizes));
gs = unique(cumsum(gs));
gs(end) = [];
hold on
plot([.5 size(ints,1)+.5]', repmat(gs+.5,1,2), 'r:', 'LineWidth', 2)
set(gca,'XTick',1:length(fitKernels), 'XTickLabel', ...
    kernelNames(kernelOrder1), 'XTickLabelRotation', 90, 'box', 'off')
c = colorbar('Ticks', -1:3);
c.Label.String = 'Kernel integral (norm.)';
c.Position = [.073 .11 .02 .815];
ylabel('Neurons')
title('Neurons'' RFs are within stimulus home position')
subplot(1,6,6)
stem(ev(valid), 'k', 'filled', 'MarkerSize', 2)
view(90,90)
ylim([0 max(ev)])
xlim([.5 sum(valid)+.5])
set(gca, 'box', 'off', 'XTick', [])
ylabel('Explained variance')

% Neuron RFs are outside stimulus position
figure('Position', [1 41 1920 1083])
colormap(cm)
subplot(1,6,1:5)
valid = dists>minDist & ~isinf(dists) & (RFdists(:,1)>0 | RFdists(:,2)>minDist);
imagesc(ints(:, valid)',[mini maxi])
gs = cellfun(@sum,mat2cell(valid, groupSizes));
gs = unique(cumsum(gs));
gs(end) = [];
hold on
plot([.5 size(ints,1)+.5]', repmat(gs+.5,1,2), 'r:', 'LineWidth', 2)
set(gca,'XTick',1:length(fitKernels), 'XTickLabel', ...
    kernelNames(kernelOrder1), 'XTickLabelRotation', 90, 'box', 'off')
c = colorbar('Ticks', -1:3);
c.Label.String = 'Kernel integral (norm.)';
c.Position = [.073 .11 .02 .815];
ylabel('Neurons')
title('Neurons'' RFs are outside stimulus home position')
subplot(1,6,6)
stem(ev(valid), 'k', 'filled', 'MarkerSize', 2)
view(90,90)
ylim([0 max(ev)])
xlim([.5 sum(valid)+.5])
set(gca, 'box', 'off', 'XTick', [])
ylabel('Explained variance')

% Neurons have no RF, but population RF is within stimulus position
figure('Position', [1 41 1920 1083])
colormap(cm)
subplot(1,6,1:5)
valid = globalDists < maxDist & isnan(dists);
imagesc(ints(:, valid)',[mini maxi])
gs = cellfun(@sum,mat2cell(valid, groupSizes));
gs = unique(cumsum(gs));
gs(end) = [];
hold on
plot([.5 size(ints,1)+.5]', repmat(gs+.5,1,2), 'r:', 'LineWidth', 2)
set(gca,'XTick',1:length(fitKernels), 'XTickLabel', ...
    kernelNames(kernelOrder1), 'XTickLabelRotation', 90, 'box', 'off')
c = colorbar('Ticks', -1:3);
c.Label.String = 'Kernel integral (norm.)';
c.Position = [.073 .11 .02 .815];
ylabel('Neurons')
title('Nonvisual neurons with population RF within stimulus home position')
subplot(1,6,6)
stem(ev(valid), 'k', 'filled', 'MarkerSize', 2)
view(90,90)
ylim([0 max(ev)])
xlim([.5 sum(valid)+.5])
set(gca, 'box', 'off', 'XTick', [])
ylabel('Explained variance')


maxDist = 2; % 2 sigma from center of stimulus position
minDist = 3;

% Population RF is within stimulus position
figure('Position', [1 41 1920 1083])
colormap(cm)
subplot(1,6,1:5)
valid = globalDists < maxDist;
imagesc(ints(:, valid)',[mini maxi])
gs = cellfun(@sum,mat2cell(valid, groupSizes));
gs = unique(cumsum(gs));
gs(end) = [];
hold on
plot([.5 size(ints,1)+.5]', repmat(gs+.5,1,2), 'r:', 'LineWidth', 2)
set(gca,'XTick',1:length(fitKernels), 'XTickLabel', ...
    kernelNames(kernelOrder1), 'XTickLabelRotation', 90, 'box', 'off')
c = colorbar('Ticks', -1:3);
c.Label.String = 'Kernel integral (norm.)';
c.Position = [.073 .11 .02 .815];
ylabel('Neurons')
title('Population RF is within stimulus home position')
subplot(1,6,6)
stem(ev(valid), 'k', 'filled', 'MarkerSize', 2)
view(90,90)
ylim([0 max(ev)])
xlim([.5 sum(valid)+.5])
set(gca, 'box', 'off', 'XTick', [])
ylabel('Explained variance')

% Population RF is outside stimulus position
figure('Position', [1 41 1920 1083])
colormap(cm)
subplot(1,6,1:5)
valid = globalDists>minDist & ~isinf(globalDists) & ...
    (RFdists_global(:,1)>0 | RFdists_global(:,2)>minDist);
imagesc(ints(:, valid)',[mini maxi])
gs = cellfun(@sum,mat2cell(valid, groupSizes));
gs = unique(cumsum(gs));
gs(end) = [];
hold on
plot([.5 size(ints,1)+.5]', repmat(gs+.5,1,2), 'r:', 'LineWidth', 2)
set(gca,'XTick',1:length(fitKernels), 'XTickLabel', ...
    kernelNames(kernelOrder1), 'XTickLabelRotation', 90, 'box', 'off')
c = colorbar('Ticks', -1:3);
c.Label.String = 'Kernel integral (norm.)';
c.Position = [.073 .11 .02 .815];
ylabel('Neurons')
title('Population RF is outside stimulus home position')
subplot(1,6,6)
stem(ev(valid), 'k', 'filled', 'MarkerSize', 2)
view(90,90)
ylim([0 max(ev)])
xlim([.5 sum(valid)+.5])
set(gca, 'box', 'off', 'XTick', [])
ylabel('Explained variance')


% plot ETAs

allETAs = cat(1, etas{:});
allPredETAs = cat(1, predETAs{:});
allETAs = allETAs(:,sn);
allPredETAs = allPredETAs(:,sn);
miniETA = min(allETAs(:));
maxiETA = max(allETAs(:));
maxiETA = max(abs([maxiETA, miniETA]));
miniPredETA = min(allPredETAs(:));
maxiPredETA = max(allPredETAs(:));
maxiPredETA = max(abs([maxiPredETA, miniPredETA]));
cm = ones(129, 3);
w = linspace(0, 1, 65)';
w(1) = [];
cm(1:64,:) = [0 0 1].*flip(w) + [1 1 1].*(1-flip(w));
cm(66:end,:) = [1 0 0].*w + [1 1 1].*(1-w);
gs = unique(cumsum(groupSizes));
gs(end) = [];
% cm = gray;
% cm = flip(cm, 1);

data = {etas, predETAs, cellfun(@minus, predETAs, etas, 'UniformOutput', false)};

% (1) All neurons
strings = {'Event triggered averages of data', ...
    'Event triggered averages of predictions', ...
    'Event triggered averages of residuals'};
for m = 1:3
    figure('Position', [1 41 1920 1083])
    annotation('textbox', [0 0.96 1 0.03], 'String', strings{m}, ...
        'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center')
    colormap(cm);
    for k = 1:length(etas)
        subplot(1,length(etas),k)
        et = data{m}{k}(:,sn)';
        et(isnan(et)) = 0;
        imagesc(times_eta{k}([1 end]), [1 size(et,1)], ...
            et, [-maxiETA maxiETA])
        hold on
        plot(times_eta{k}([1 end])', repmat(gs'+.5,1,2), 'r:', 'LineWidth', 2)
        plot([0 0], [.5 size(et,1)+.5], 'k:', 'LineWidth', 2)
        xlabel('Time (s)')
        if k==1
            ylabel('# Neurons')
        end
        set(gca, 'YTick', [])
        set(gca, 'box', 'off')
        title(etaNames{etaOrder(k)})
    end
    c = colorbar('Position', [.925 .11 .014 .815]);
    c.Label.String = 'normalised response';
end

% (2) Neurons with population RFs within stimulus position
valid = globalDists < maxDist;
strings = {'Event triggered averages of data (population RF is within stimulus position)', ...
    'Event triggered averages of predictions (population RF is within stimulus position)', ...
    'Event triggered averages of residuals (population RF is within stimulus position)'};
for m = 1:3
    figure('Position', [1 41 1920 1083])
    annotation('textbox', [0 0.96 1 0.03], 'String', strings{m}, ...
        'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center')
    colormap(cm);
    for k = 1:length(etas)
        subplot(1,length(etas),k)
        et = data{m}{k}(:,sn)';
        et = et(valid,:);
        et(isnan(et)) = 0;
        imagesc(times_eta{k}([1 end]), [1 size(et,1)], ...
            et, [-maxiETA maxiETA])
        gs = cellfun(@sum,mat2cell(valid, groupSizes));
        gs = unique(cumsum(gs));
        gs(end) = [];
        hold on
        plot(times_eta{k}([1 end])', repmat(gs+.5,1,2), 'r:', 'LineWidth', 2)
        plot([0 0], [.5 size(et,1)+.5], 'k:', 'LineWidth', 2)
        xlabel('Time (s)')
        if k==1
            ylabel('# Neurons')
        end
        set(gca, 'YTick', [])
        set(gca, 'box', 'off')
        title(etaNames{etaOrder(k)})
    end
    c = colorbar('Position', [.925 .11 .014 .815]);
    c.Label.String = 'normalised response';
end

% (3) Neurons with population RFs outside stimulus position
valid = globalDists>minDist & ~isinf(globalDists) & ...
    (RFdists_global(:,1)>0 | RFdists_global(:,2)>minDist);
strings = {'Event triggered averages of data (population RF is outside stimulus position)', ...
    'Event triggered averages of predictions (population RF is outside stimulus position)', ...
    'Event triggered averages of residuals (population RF is outside stimulus position)'};
for m = 1:3
    figure('Position', [1 41 1920 1083])
    annotation('textbox', [0 0.96 1 0.03], 'String', strings{m}, ...
        'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center')
    colormap(cm);
    for k = 1:length(etas)
        subplot(1,length(etas),k)
        et = data{m}{k}(:,sn)';
        et = et(valid,:);
        et(isnan(et)) = 0;
        imagesc(times_eta{k}([1 end]), [1 size(et,1)], ...
            et, [-maxiETA maxiETA])
        gs = cellfun(@sum,mat2cell(valid, groupSizes));
        gs = unique(cumsum(gs));
        gs(end) = [];
        hold on
        plot(times_eta{k}([1 end])', repmat(gs+.5,1,2), 'r:', 'LineWidth', 2)
        plot([0 0], [.5 size(et,1)+.5], 'k:', 'LineWidth', 2)
        xlabel('Time (s)')
        if k==1
            ylabel('# Neurons')
        end
        set(gca, 'YTick', [])
        set(gca, 'box', 'off')
        title(etaNames{etaOrder(k)})
    end
    c = colorbar('Position', [.925 .11 .014 .815]);
    c.Label.String = 'normalised response';
end