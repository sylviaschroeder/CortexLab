%% Folders
folderData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\NPY\task_2p';
folderTools = 'C:\STORAGE\workspaces';
folderCode = 'C:\dev\workspace\CortexLab';
folderPlots = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\wheelTask\pupil';

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(fullfile(folderCode))

%% Align pupil to task events
pre = 2;
post = 5;

subjects = cell(0,1);
dates = cell(0,1);
algnmnts = cell(0, 3);
timeAlgn = cell(0,1);
binSz = NaN(0,1);
stimTimes = cell(0,1);
contrasts = cell(0,1);
choices = cell(0,1);
feedbacks = cell(0,1);

subjDirs = dir(fullfile(folderData, 'SS*'));
for subj = 1:length(subjDirs)
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folderData, name, '2*'));
    for dt = 1:length(dateDirs)
        date = dateDirs(dt).name;
        if ~isfile(fullfile(folderData, name, date, 'eye.diameter.npy'))
            continue
        end
        subjects{end+1,1} = name;
        dates{end+1,1} = date;
        
        pupil = readNPY(fullfile(folderData, name, date, 'eye.diameter.npy'));
        pupilTime = readNPY(fullfile(folderData, name, date, 'eye.timestamps.npy'));
        stimT = readNPY(fullfile(folderData, name, date, '_ibl_trials.stim_intervals.npy'));
        contrastL = readNPY(fullfile(folderData, name, date, '_ibl_trials.contrastLeft.npy'));
        contrastR = readNPY(fullfile(folderData, name, date, '_ibl_trials.contrastRight.npy'));
        beeps = readNPY(fullfile(folderData, name, date, '_ibl_trials.goCue_times.npy'));
        choice = readNPY(fullfile(folderData, name, date, '_ibl_trials.choice.npy'));
        fb = readNPY(fullfile(folderData, name, date, '_ibl_trials.feedbackType.npy'));
        fbTimes = readNPY(fullfile(folderData, name, date, '_ibl_trials.feedback_times.npy'));
        
        stimTimes{end+1,1} = stimT;
        contrasts{end+1,1} = [contrastL, contrastR];
        choices{end+1,1} = choice;
        feedbacks{end+1,1} = fb;
        feedbackTimes{end+1,1} = fbTimes;
        
        sampleDur = diff(pupilTime(1:2));
        preSamples = round(pre / sampleDur);
        postSamples = round(post / sampleDur);
        
        binSz(end+1,1) = sampleDur;
        timeAlgn{end+1,1} = (-preSamples:postSamples) .* sampleDur;
        
        events = [stimT(:,1), beeps, fbTimes];
        n = length(subjects);
        for ev = 1:size(events,2)
            ind = NaN(size(events,1),1);
            for t = 1:size(events,1)
                ind(t) = find(pupilTime >= events(t,ev), 1, 'first');
            end
            ind = ind + (-preSamples : postSamples);
            ind2 = ind;
            invalid = ind < 1 | ind > length(pupilTime);
            ind2(invalid) = 1;
            algnmnts{n,ev} = pupil(ind2);
            algnmnts{n,ev}(invalid) = NaN;
        end
    end
end

%% Plot pupil aligned to task events
plotSingleTrials = false;
pos = [[500 550 560 420]; [1100 550 560 420]; ...
    [500 100 560 420]; [1100 100 560 420]];
evLabels = {'stimulus onset', 'beep', 'feedback'};
evColors = lines(length(evLabels));
trialTypeLabels = {{'Stimulus contrast > 0', 'Stimulus contrast == 0'}, ...
    {'Go trials', 'NoGo trials'}, ...
    {'Rewarded trials', 'Non-rewarded trials'}, ...
    {'Go & rewarded trials', 'Go & non-rewarded trials', ...
    'NoGo & rewarded trials', 'NoGo & non-rewarded trials'}};

for iExp = 1:length(subjects)
    trialTypes = {[contrasts{iExp}(:,1)>0 | contrasts{iExp}(:,2)>0, ...
        contrasts{iExp}(:,1)==0 & contrasts{iExp}(:,2)==0], ... % contrast>0, contrast=0 (both sides)
        [choices{iExp} < 3, choices{iExp} == 3], ... % Go trials, NoGo trials
        [feedbacks{iExp} == 1, feedbacks{iExp} == 0], ... % correct, incorrect
        [choices{iExp} < 3 & feedbacks{iExp} == 1, choices{iExp} < 3 & feedbacks{iExp} == 0, ...
        choices{iExp} == 3 & feedbacks{iExp} == 1, choices{iExp} == 3 & feedbacks{iExp} == 0]}; % Go/NoGo & correct/incorrect
    mini = Inf;
    maxi = -Inf;
    
    ax = [];
    figs = [];
    for ev = 1:size(algnmnts,2)
        for type = 1:length(trialTypes)
            for cond = 1:size(trialTypes{type},2)
                ppl = algnmnts{iExp,ev}(trialTypes{type}(:,cond),:);
                figs(end+1) = figure('Position', pos(cond,:));
                hold on
                m = nanmean(ppl,1);
                if plotSingleTrials
                    plot(timeAlgn{iExp}, ppl', 'k')
                else
                    sem = nanstd(ppl) ./ sqrt(sum(~isnan(ppl)));
                    fill([timeAlgn{iExp} flip(timeAlgn{iExp})], [m + sem, flip(m - sem)], ...
                        'k', 'EdgeColor', 'none', 'FaceColor', evColors(ev,:), ...
                        'FaceAlpha', 0.3)
                    mini = min([mini, m-sem]);
                    maxi = max([maxi, m+sem]);
                end
                plot(timeAlgn{iExp}, m, 'Color', evColors(ev,:), 'LineWidth', 2)
                xlim(timeAlgn{iExp}([1 end]))
                ylim([mini maxi])
                xlabel(sprintf('Time from %s (s)', evLabels{ev}))
                ylabel('Pupil size')
                title(trialTypeLabels{type}{cond})
                ax(end+1) = gca;
            end
        end
    end
    if plotSingleTrials
        mini = min(min(cat(1, algnmnts{iExp,:})));
        maxi = max(max(cat(1, algnmnts{iExp,:})));
    end
    set(ax, 'YLim', [mini maxi])
    folder = fullfile(folderPlots, 'pupilAlignedToEvents', ...
        subjects{iExp}, dates{iExp});
    if ~isfolder(folder)
        mkdir(folder)
    end
    for k = 1:length(ax)
        saveas(figs(k), fullfile(folder, sprintf('figure%02d.png', k)))
        close(figs(k))
    end
end

%% Plot pupil baseline, Go, NoGo, correct Go, incorrect NoGo rate across time
baseTime = [-1 0];
win = normpdf(-7:7, 0, 2);
for iExp = 1:length(subjects)
    bt = timeAlgn{iExp} >= baseTime(1) & timeAlgn{iExp} <= baseTime(2);
    baseline = nanmean(algnmnts{iExp,1}(:, bt), 2);
    go = choices{iExp} < 3;
    nogo = choices{iExp} == 3;
    corrGo = choices{iExp} < 3 & feedbacks{iExp} == 1;
    incorrNogo = choices{iExp} == 3 & feedbacks{iExp} == 0;
    corr = feedbacks{iExp} == 1;
    
    ax = [0 0];
    figure
    subplot(2,1,1)
    plot(baseline, 'k', 'LineWidth', 2)
%     plot(conv(baseline, win, 'same'), 'k', 'LineWidth', 2)
    ax(1) = gca;
    subplot(2,1,2)
    hold on
    h = [0 0 0 0];
    h(1) = plot(conv(go, win, 'same'), 'LineWidth', 2);
    h(2) = plot(conv(nogo, win, 'same'), 'LineWidth', 2);
    h(3) = plot(conv(corrGo, win, 'same'), 'LineWidth', 2);
    h(4) = plot(conv(incorrNogo, win, 'same'), 'LineWidth', 2);
    legend(h, {'Go','Nogo','corr Go','incorr Nogo'})
    ax(2) = gca;
    linkaxes(ax, 'x')
    xlim([1 length(baseline)])
end

%% Plot scatters: pupil baseline vs performance
baseTime = [-1 0];
cols = lines(4);
outcomes = {'hit','wrong dir.','false alarm','corr. rej.','miss'};
for iExp = 1:length(subjects)
    bt = timeAlgn{iExp} >= baseTime(1) & timeAlgn{iExp} <= baseTime(2);
    baseline = nanmean(algnmnts{iExp,1}(:, bt), 2);
    hits = choices{iExp} < 3 & feedbacks{iExp} == 1;
    incorrectDirs = choices{iExp} < 3 & feedbacks{iExp} == 0 & max(contrasts{iExp}, [], 2) > 0;
    falseAlarms = choices{iExp} < 3 & feedbacks{iExp} == 0 & max(contrasts{iExp}, [], 2) == 0;
    corrRejects = choices{iExp} == 3 & feedbacks{iExp} == 1;
    misses = choices{iExp} == 3 & feedbacks{iExp} == 0;
    performance = sum([hits, incorrectDirs, falseAlarms, corrRejects, misses] * (1:5), 2);
    
    % previous trial performance vs pupil baseline
    figure %('Position', [4 543 1913 420])
    scatter(performance(1:end-1)+randn(length(corr)-1,1).*0.1, baseline(2:end))
    set(gca, 'XTick', 1:size(performance,2), 'XTickLabel', outcomes)
    xlabel('Previous trial performance')
    ylabel('Pupil diameter')
end

% 1. align pupil to task events (stimulus onset, beep, reward/penalty, wheel
% movement(?), ...)
% 2. Like 1. but differentiate between conditions (contrast of stimulus,
% side of stimulus, correct/incorrect move, correct/incorrect nogo, 
% left/right, ...)
% 3. quantify how pupil baseline/dilation changes across different trials
% (correct/incorrect, move/nogo, outcome of previous trial)