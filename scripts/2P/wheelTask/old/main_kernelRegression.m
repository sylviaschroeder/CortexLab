%% Folders and parameters
dataFolder = 'C:\RESULTS\Matteo A\Q';
resultsFolder = 'C:\RESULTS\wheelTask\kernelFit';
traceFolder = 'C:\RESULTS\wheelTask\eventTraces';
wheelFolder = 'C:\RESULTS\wheelTask\wheelVelocity';
configFolder = '\\ZSERVER.cortexlab.net\Code\Rigging\config';

rigPC = 'ZURPRISE';
velocityThreshold = 10; % mm/s
smoothing = .15; % s
minBetweenMoves = .2; % s
reactionTime = .05; % s


%% Load data
db_wheelTask;

dataset = 3;
data = load(fullfile(dataFolder, db(dataset).subject, db(dataset).date, num2str(db(dataset).exp), ...
    sprintf('%s_%d_%s_integrated.mat', db(dataset).date, db(dataset).exp, db(dataset).subject)));
q = data.q;

%% Get wheel movements
data = load(fullfile(configFolder, rigPC, 'hardware.mat'), 'mouseInput');
mmFactor = data.mouseInput.MillimetresFactor;

rawPos = wheel.correctCounterDiscont(q.timeline.rawDAQData(:, ...
    strcmp({q.timeline.hw.inputs.name}, 'rotaryEncoder')));
sr = 1 / median(diff(q.timeline.rawDAQTimestamps));
time = round(q.timeline.rawDAQTimestamps(1)*sr)/sr : 1/sr : ...
    round(q.timeline.rawDAQTimestamps(end)*sr)/sr;
pos = interp1(q.timeline.rawDAQTimestamps, rawPos, time, 'pchip');
pos = pos .* mmFactor;
[vel, acc] = wheel.computeVelocity(pos, round(smoothing * sr), sr);
[allMoves, allAmps] = wheel.findAllMoves(time, vel, pos, velocityThreshold, ...
    minBetweenMoves);
allOnsets = allMoves(1,:);
stimTimes = sort([q.stimTimes(1,q.trialConditions(:,3)~=0), ...
    q.beepTimes(q.trialConditions(:,3)==0)] + reactionTime); % stim onset or beep (if contrast 0) + minimum reaction time
moves = repmat(allOnsets', 1, length(stimTimes));
moves(bsxfun(@gt, moves, q.choiceTimes)) = NaN;
afterStim = bsxfun(@minus, moves, stimTimes);
afterStim(afterStim < 0) = NaN;
[t,ind] = min(afterStim, [], 1);
afterStimMoves = allOnsets(ind);
afterStimMoves(isnan(t)) = NaN;
afterStimAmps = allAmps(ind);
afterStimAmps(isnan(t)) = NaN;

%% Plot wheel velocity and calcium trace for each neuron
vel = wheel.computeVelocity(pos, round(1 * sr), sr);
vel = interp1(time, vel, q.frameTimes, 'pchip');

contrasts = unique(q.trialConditions(:,3));
contrasts = setdiff(contrasts, 0);
stimTimes = cell(1, length(contrasts));
for c = 1:length(contrasts)
    stimTimes{c} = q.stimTimes(:,q.trialConditions(:,3)==contrasts(c));
end

folder = fullfile(wheelFolder, sprintf('%s_%s_%d', db(dataset).subject, ...
    db(dataset).date, db(dataset).exp));
if ~isdir(folder)
    mkdir(folder)
end

colors = jet(length(contrasts));
subplots = {1:3,5:7};
mini = [0 0];
maxi = [0 0];
mini(1) = min(vel);
maxi(1) = max(vel);
rng = maxi(1)-mini(1);
mini(1) = mini(1)-.05*rng;
maxi(1) = maxi(1)+.05*rng;
for iCell = 1:size(q.cellMat,2)
    F = smooth(q.cellMat(:,iCell), 5);
    mini(2) = min(F);
    maxi(2) = max(F);
    rng = maxi(2)-mini(2);
    mini(2) = mini(2)-.05*rng;
    maxi(2) = maxi(2)+.05*rng;
    figure('Position', [2 380 1917 735])
    subplot(2,4,subplots{1})
    hold on
    subplot(2,4,subplots{2})
    hold on
    h = zeros(1, length(contrasts));
    for c = 1:length(contrasts)
        for sb = 1:2
            subplot(2,4,subplots{sb})
            p = fill(stimTimes{c}([1 1 2 2],:), [mini(sb) maxi(sb) maxi(sb) mini(sb)]', 'k', ...
                'FaceColor', colors(c,:), 'FaceAlpha', .1, 'EdgeColor', 'none');
            h(c) = p(1);
        end
    end
    lgn = legend(h, {num2str(contrasts)}, 'Position', [.07 .285 .03 .16]);
    title(lgn, 'Contrast')
    subplot(2,4,subplots{1})
    plot(q.frameTimes, vel, 'k')
    ylim([mini(1) maxi(1)])
    ax1 = gca;
    title('Wheel velocity')
    ylabel('Velocity (mm/s)')
    subplot(2,4,subplots{2})
    plot(q.frameTimes, F, 'k')
    ylim([mini(2) maxi(2)])
    ax2 = gca;
    title('Calcium trace')
    xlabel('Time (s)')
    ylabel('\DeltaF/F')
    linkaxes([ax1 ax2], 'x')
    xlim(time([1 end]))
    
    subplot(2,4,4)
    plot(vel, F, 'k.')
    xlabel('Wheel velocity')
    ylabel('Calcium trace')
    
    maxLag = round(.1*sr);
    [r, lags] = xcorr(vel, F, maxLag, 'coeff');
    subplot(2,4,8)
    plot(lags ./ sr, r, 'k')
    xlabel('Time lag (s)')
    ylabel('Cross-correlation')
    
    f = gcf;
    f.PaperPositionMode = 'auto';
    print(fullfile(folder, sprintf('wheelVel_%03d.tiff', iCell)), ...
        '-dtiff','-r0')
    close gcf
end

%% Define events and windows
stimWin = [0 3];
beepWin = [0 2];
moveWin = [-1 2];
choiceWin = [-1 1];
fbWin = [0 2];

contrasts = unique(q.trialConditions(:,3));
eventTimes0 = cell(1, 3*length(contrasts));
for c = 1:length(contrasts)
    eventTimes0{c} = q.stimTimes(1,q.trialConditions(:,3)==contrasts(c));
    eventTimes0{c+length(contrasts)} = afterStimMoves( ...
        q.trialConditions(:,3)==contrasts(c) & ~isnan(afterStimMoves'));
    eventTimes0{c+2*length(contrasts)} = ...
        q.choiceTimes(q.trialConditions(:,3)==contrasts(c));
end
windows0 = [mat2cell(repmat(stimWin, length(contrasts), 1), ...
    ones(1, length(contrasts)), 2); ...
    mat2cell(repmat(moveWin, length(contrasts), 1), ...
    ones(1, length(contrasts)), 2); ...
    mat2cell(repmat(choiceWin, length(contrasts), 1), ...
    ones(1, length(contrasts)), 2)];

% eventTimes1 = {q.stimTimes(1,q.trialConditions(:,3)<0), ...% stim onset on left
%     q.stimTimes(1,q.trialConditions(:,3)>0), ...           % stim onset on right
%     q.beepTimes, ...                                       % beep
%     afterStimMoves(afterStimAmps>0), ...                   % movement onsets to left
%     afterStimMoves(afterStimAmps<0), ...                   % movement onsets to right
%     q.choiceTimes(q.trialConditions(:,4)==1), ...          % time of choice to left
%     q.choiceTimes(q.trialConditions(:,4)==2), ...          % time of choice to right
%     q.choiceTimes(q.trialConditions(:,4)==3), ...          % time of choice to hold still
%     q.choiceTimes(q.trialConditions(:,5)==1), ...          % time of reward
%     q.choiceTimes(q.trialConditions(:,5)==-1)};            % time of neg. feedback
eventTimes1 = {q.stimTimes(1,q.trialConditions(:,3)<0), ...% stim onset on left
    q.stimTimes(1,q.trialConditions(:,3)>0), ...           % stim onset on right
    q.beepTimes, ...                                       % beep
    afterStimMoves(~isnan(afterStimMoves')&q.trialConditions(:,3)<0), ... % movement onsets to left
    afterStimMoves(~isnan(afterStimMoves')&q.trialConditions(:,3)>0), ... % movement onsets to right
    q.choiceTimes(q.trialConditions(:,4)==1), ...          % time of choice to left
    q.choiceTimes(q.trialConditions(:,4)==2), ...          % time of choice to right
    q.choiceTimes(q.trialConditions(:,4)==3), ...          % time of choice to hold still
    q.choiceTimes(q.trialConditions(:,5)==1), ...          % time of reward
    q.choiceTimes(q.trialConditions(:,5)==-1)};            % time of neg. feedback
eventTimes2 = {q.stimTimes(1,q.trialConditions(:,3)~=0 & q.trialConditions(:,5)==1), ...% stim onset + correct
    q.stimTimes(1,q.trialConditions(:,3)~=0 & q.trialConditions(:,5)==-1), ...          % stim onset + incorrect
    q.beepTimes(q.trialConditions(:,5)==1), ...                                         % beep + correct
    q.beepTimes(q.trialConditions(:,5)==-1), ...                                        % beep + incorrect
    afterStimMoves(q.trialConditions(:,5)==1 & ~isnan(afterStimMoves')), ...             % move onset + correct
    afterStimMoves(q.trialConditions(:,5)==-1 & ~isnan(afterStimMoves')), ...            % move onset + incorrect
    q.choiceTimes(q.trialConditions(:,4)<3 & q.trialConditions(:,5)==1), ...            % time of choice (moved) + correct
    q.choiceTimes(q.trialConditions(:,4)<3 & q.trialConditions(:,5)==-1), ...           % time of choice (moved) + incorrect
    q.choiceTimes(q.trialConditions(:,4)==3 & q.trialConditions(:,5)==1), ...           % time of choice (still) + correct
    q.choiceTimes(q.trialConditions(:,4)==3 & q.trialConditions(:,5)==-1)};             % time of choice (still) + incorrect

windows1 = {stimWin, stimWin, beepWin, moveWin, moveWin, choiceWin, ...
    choiceWin, choiceWin, fbWin, fbWin};
windows2 = {stimWin, stimWin, beepWin, beepWin, moveWin, moveWin, ...
    choiceWin, choiceWin, choiceWin, choiceWin};
lambda = .01;

%% Perform linear regression and get event triggered averages
[A, numSamples, windowTimes0] = ...
    krnl.getToeplitz(q.frameTimes, eventTimes0, windows0);
eta0 = krnl.getETA(q.cellMat, A, numSamples);

[fitKernels, eta1, predSig, predETA1, windowTimes1, snippets1] = ...
    krnl.kernelRegression(q.cellMat, q.frameTimes, eventTimes1, ...
    windows1, lambda);

[A, numSamples, windowTimes2] = ...
    krnl.getToeplitz(q.frameTimes, eventTimes2, windows2);
eta2 = krnl.getETA(q.cellMat, A, numSamples);
predETA2 = krnl.getETA(predSig, A, numSamples);

%% Define plotting
subplots0 = {1:length(contrasts), [], length(contrasts)+(1:length(contrasts)), ...
    2*length(contrasts)+(1:length(contrasts)), []};
subplots1 = {[1 2], 3, [4 5], [6 7 8], [9 10]};
subplots2 = {[1 2], [3 4], [5 6], [7 8], [9 10]};
titles0 = {'Stimulus onset', [], 'Movement onset', 'Time of choice'};
titles1 = {'Stimulus onset', 'Beep', 'Movement onset', 'Time of choice', ...
    'Time of feedback'};
titles2 = {'Stimulus onset', 'Beep', 'Movement onset', 'Time of choice when moved', ...
    'Time of choice when still'};
labels0 = {num2cell(num2str(contrasts(:)),2), [], num2cell(num2str(contrasts(:)),2), ...
    num2cell(num2str(contrasts(:)),2)};
labels1 = {{'left','right'},[],{'left','right'},{'left','right','still'},...
    {'positive','negative'}};
labels2 = {{'correct','incorrect'},{'correct','incorrect'},{'correct','incorrect'}, ...
    {'correct','incorrect'},{'correct','incorrect'}};

allKernels = cat(1, fitKernels{:});
allETAs0 = cat(1, eta0{:});
allETAs1 = cat(1, eta1{:}, predETA1{:});
allETAs2 = cat(1, eta2{:}, predETA2{:});

%% Plot results for each cell
folder = fullfile(resultsFolder, sprintf('%s_%s_%d', db(dataset).subject, ...
    db(dataset).date, db(dataset).exp));
if ~isdir(folder)
    mkdir(folder)
end

cols = lines(max(cellfun(@length, subplots1)));
for iCell = 1:size(q.cellMat,2)
    figure('Position', [1 41 1920 1083])
    
    maxi = max(allETAs0(:,iCell));
    mini = min(allETAs0(:,iCell));
    rng = maxi - mini;
    maxi = maxi + .05*rng;
    mini = mini - .05*rng;
    c0 = jet(length(contrasts));
    k = 0;
    for sp = 1:length(subplots0)
        if isempty(subplots0{sp})
            continue
        end
        subplot(4,length(subplots1),sp)
        hold on
        for c = 1:length(contrasts)
            plot(windowTimes0{c+k*length(contrasts)}, ...
                eta0{c+k*length(contrasts)}(:,iCell), 'Color', c0(c,:))
        end
        legend(labels0{sp})
        xlim(windowTimes0{1+k*length(contrasts)}([1 end]))
        ylim([mini maxi])
        title(titles0{sp})
        ylabel('Event triggered average')
        k = k + 1;
    end
    
    maxi = max(allKernels(:,iCell));
    mini = min(allKernels(:,iCell));
    if maxi-mini==0
        f = gcf;
        f.PaperPositionMode = 'auto';
        print(fullfile(folder, sprintf('kernelFit_%03d.jpg', iCell)), ...
            '-djpeg','-r0')
        close gcf
        continue
    end
    rng = maxi - mini;
    maxiK = maxi + .05*rng;
    miniK = mini - .05*rng;
    maxi = max(allETAs1(:,iCell));
    mini = min(allETAs1(:,iCell));
    rng = maxi - mini;
    maxi = maxi + .05*rng;
    mini = mini - .05*rng;
    for sp = 1:length(subplots1)
        subplot(4, length(subplots1), length(subplots1) + sp)
        hold on
        for l = 1:length(subplots1{sp})
            m = subplots1{sp}(l);
            plot(windowTimes1{m}, fitKernels{m}(:,iCell))
        end
        if ~isempty(labels1{sp})
            legend(labels1{sp})
        end
        xlim(windowTimes1{m}([1 end]))
        ylim([miniK maxiK])
        title(titles1{sp})
        ylabel('Fitted kernel')
        
        subplot(4, length(subplots1), 2*length(subplots1)+sp)
        h = zeros(1,length(subplots1{sp}));
        hold on
        for l = 1:length(subplots1{sp})
            m = subplots1{sp}(l);
            h(l) = plot(windowTimes1{m}, eta1{m}(:,iCell), 'Color', cols(l,:));
            plot(windowTimes1{m}, predETA1{m}(:,iCell), ':', 'Color', cols(l,:))
        end
        if ~isempty(labels1{sp})
            legend(h,labels1{sp})
        end
        xlim(windowTimes1{m}([1 end]))
        ylim([mini maxi])
        xlabel('Time from event onset')
        ylabel('Event triggered average')
    end
    
    maxi = max(allETAs2(:,iCell));
    mini = min(allETAs2(:,iCell));
    rng = maxi - mini;
    maxi = maxi + .05*rng;
    mini = mini - .05*rng;
    for sp = 1:length(subplots2)
        subplot(4, length(subplots2), 3*length(subplots2)+sp)
        h = zeros(1,length(subplots2{sp}));
        hold on
        for l = 1:length(subplots2{sp})
            m = subplots2{sp}(l);
            h(l) = plot(windowTimes2{m}, eta2{m}(:,iCell), 'Color', cols(l,:));
            plot(windowTimes2{m}, predETA2{m}(:,iCell), ':', 'Color', cols(l,:))
        end
        if ~isempty(labels2{sp})
            legend(h,labels2{sp})
        end
        xlim(windowTimes2{m}([1 end]))
        ylim([mini maxi])
        xlabel('Time from event onset')
        ylabel('Event triggered average')
        title(titles2{sp})
    end
    f = gcf;
    f.PaperPositionMode = 'auto';
    print(fullfile(folder, sprintf('kernelFit_%03d.jpg', iCell)), ...
        '-djpeg','-r0')
    close gcf
    
%     figure('position',[3 678 1916 420])
%     plot(q.frameTimes, q.cellMat(:,iCell))
%     hold on
%     plot(q.frameTimes, predSig(:,iCell))
%     legend('data','prediction')
%     xlim(q.frameTimes([1 end]))
%     title('Calcium trace')
%     xlabel('Time (s)')
%     ylabel('\DeltaF/F')
%     f = gcf;
%     f.PaperPositionMode = 'auto';
%     print(fullfile(folder, sprintf('prediction_%03d.jpg', iCell)), ...
%         '-djpeg','-r0')
%     close gcf
end

%% Plot event aligned traces for each cell
% Define events and windows
stimWin = [-.5 3];
choiceWin = [-3 1];
eventTimes = {q.stimTimes(1,q.trialConditions(:,3)<0), ...% stim onset on left
    q.stimTimes(1,q.trialConditions(:,3)>0), ...          % stim onset on right
    q.choiceTimes(q.trialConditions(:,4)==1), ...         % time of choice to left
    q.choiceTimes(q.trialConditions(:,4)==2)};            % time of choice to right
windows = {stimWin, stimWin, choiceWin, choiceWin};
sortCrit = {@(x)min(x(x>0)); @(x)min(x(x>0)); @(x)min(abs(x(x<0))); @(x)min(abs(x(x<0)))};
titles = {'Stimulus onset - left','Stimulus onset - right','Choice - left','Choice - right'};
plotRows = 2;
plotCols = 2;
subplots = 1:4;

% Get choice times
choiceTimes = cell(1, length(eventTimes));
choiceTimes{1} = q.choiceTimes(q.trialConditions(:,3)<0) - ...
    q.stimTimes(1,q.trialConditions(:,3)<0);
choiceTimes{1}(choiceTimes{1}>windows{1}(2)) = NaN;
choiceTimes{2} = q.choiceTimes(q.trialConditions(:,3)>0) - ...
    q.stimTimes(1,q.trialConditions(:,3)>0);
choiceTimes{2}(choiceTimes{2}>windows{2}(2)) = NaN;

% Get stim onset times
stimTimes = cell(1, length(eventTimes));
stimTimes{3} = q.stimTimes(1,q.trialConditions(:,4)==1) - ...
    q.choiceTimes(q.trialConditions(:,4)==1);
stimTimes{3}(stimTimes{3}<windows{3}(1)) = NaN;
stimTimes{4} = q.stimTimes(1,q.trialConditions(:,4)==2) - ...
    q.choiceTimes(q.trialConditions(:,4)==2);
stimTimes{4}(stimTimes{4}<windows{4}(1)) = NaN;

% Get movement onsets for each trial of each event category
moveOnTimes = cell(1, length(eventTimes));
moveOffTimes = cell(1, length(eventTimes));
for ev = 1:length(eventTimes)
    moves = repmat(allOnsets', 1, length(eventTimes{ev}));
    moves = bsxfun(@minus, moves, eventTimes{ev});
    moves(moves < windows{ev}(1) | moves > windows{ev}(2)) = NaN;
    n = sum(~isnan(moves),1);
    moveOnTimes{ev} = mat2cell(moves(~isnan(moves)), n);
    moves = repmat(allMoves(2,:)', 1, length(eventTimes{ev}));
    moves = bsxfun(@minus, moves, eventTimes{ev});
    moves(moves < windows{ev}(1) | moves > windows{ev}(2)) = NaN;
    n = sum(~isnan(moves),1);
    moveOffTimes{ev} = mat2cell(moves(~isnan(moves)), n);
end

% Sort trials
trialOrders = cell(1, length(eventTimes));
for ev = 1:length(eventTimes)
    relMoveTimes = NaN(length(moveOnTimes{ev}), 1);
    for k = 1:length(moveOnTimes{ev})
        if isempty(sortCrit{ev}(moveOnTimes{ev}{k}))
            continue
        end
        relMoveTimes(k) = sortCrit{ev}(moveOnTimes{ev}{k});
    end
    [~,trialOrders{ev}] = sort(relMoveTimes);
end

% Get aligned traces
[A, numSamples, windowTimes] = ...
    krnl.getToeplitz(q.frameTimes, eventTimes, windows);
[~, traces] = krnl.getETA(q.cellMat, A, numSamples);

% Plot for each cell
folder = fullfile(traceFolder, sprintf('%s_%s_%d', db(dataset).subject, ...
    db(dataset).date, db(dataset).exp));
if ~isdir(folder)
    mkdir(folder)
end

cols = 'yymm';
for iCell = 1:size(q.cellMat,2)
    figure('Position', [1 41 1920 1083])
    allTraces = [];
    for ev = 1:length(eventTimes)
        allTraces = [allTraces; reshape(traces{ev}(:,iCell,:),[],1)];
    end
    mini = min(allTraces);
    maxi = max(allTraces);
    for ev = 1:length(eventTimes)
        subplot(plotRows, plotCols, subplots(ev))
        imagesc(windowTimes{ev}([1 end]), [1 size(traces{ev},3)], ...
            squeeze(traces{ev}(:,iCell,trialOrders{ev}))', [mini maxi])
        colormap gray
        hold on
        for tr = 1:length(trialOrders{ev})
            t = moveOnTimes{ev}{trialOrders{ev}(tr)};
            plot(t, ones(1,length(t)).*tr, '.r')
            t = moveOffTimes{ev}{trialOrders{ev}(tr)};
            plot(t, ones(1,length(t)).*tr, '.c')
            
            if ~isempty(stimTimes{ev})
                plot(stimTimes{ev}(trialOrders{ev}(tr)),tr,'>y','MarkerSize',4)
            end
            if~isempty(choiceTimes{ev})
                plot(choiceTimes{ev}(trialOrders{ev}(tr)),tr,'<m','MarkerSize',4)
            end
        end
        plot([0 0],[.5 size(traces{ev},3)+.5],[cols(ev) ':'],'LineWidth',2)
        title(titles{ev})
        xlabel('Time (s)')
        ylabel('# Events')
    end
    
    f = gcf;
    f.PaperPositionMode = 'auto';
    print(fullfile(folder, sprintf('eventTraces_%03d.tiff', iCell)), ...
        '-dtiff','-r0')
    close gcf
end

%% Plot population results
