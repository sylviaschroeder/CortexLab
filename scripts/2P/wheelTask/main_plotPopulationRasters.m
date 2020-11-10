%% Folders
folderData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\NPY\task_2p';
folderTools = 'C:\STORAGE\workspaces';
folderCode = 'C:\dev\workspace\CortexLab';
folderPlots = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\wheelTask\Plots';

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(fullfile(folderCode))

%% Parameters
sigma = 1;

%% Plot population rasters
subjDirs = dir(fullfile(folderData, 'SS*'));
for subj = 1:length(subjDirs)
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folderData, name, '2*'));
    for dt = 1:length(dateDirs)
        date = dateDirs(dt).name;
        if ~isfile(fullfile(folderData, name, date, 'eye.diameter.npy'))
            continue
        end
        
        % load data
        caTraces = readNPY(fullfile(folderData, name, date, '_ss_2pCalcium.dff.npy'));
        caTime = readNPY(fullfile(folderData, name, date, '_ss_2pCalcium.timestamps.npy'));
        planes = readNPY(fullfile(folderData, name, date, '_ss_2pRois._ss_2pPlanes.npy'));
        planeDelays = readNPY(fullfile(folderData, name, date, '_ss_2pPlanes.delay.npy'));
        pupil = readNPY(fullfile(folderData, name, date, 'eye.diameter.npy'));
        pupilTime = readNPY(fullfile(folderData, name, date, 'eye.timestamps.npy'));
        stimT = readNPY(fullfile(folderData, name, date, '_ibl_trials.stim_intervals.npy'));
        contrastL = readNPY(fullfile(folderData, name, date, '_ibl_trials.contrastLeft.npy'));
        contrastR = readNPY(fullfile(folderData, name, date, '_ibl_trials.contrastRight.npy'));
        wheelVel = readNPY(fullfile(folderData, name, date, '_ibl_wheel.velocity.npy'));
        wheelTime = readNPY(fullfile(folderData, name, date, '_ibl_wheel.timestamps.npy'));
        cueInteractiveDelays = readNPY(fullfile(folderData, name, date, '_ibl_trials.cueInteractiveDelay.npy'));
        beeps = readNPY(fullfile(folderData, name, date, '_ibl_trials.goCue_times.npy'));
        choice = readNPY(fullfile(folderData, name, date, '_ibl_trials.choice.npy'));
        fb = readNPY(fullfile(folderData, name, date, '_ibl_trials.feedbackType.npy'));
        fbTimes = readNPY(fullfile(folderData, name, date, '_ibl_trials.feedback_times.npy'));
        
        % probability of go choice over time
        probGo = choice < 3;
        probGo = smooth(probGo, 9);
        % probability of correct choice
        probCorrect = smooth(fb, 9);
        % interpolate and convolve pupil trace
        ppl = traces.downSample(pupil, pupilTime, sigma, caTime);
        % convolve neurons
        dtime = median(diff(caTime));
        sigSamples = round(sigma / dtime);
        win = normpdf(-4*sigSamples : 4*sigSamples, 0, sigSamples);
        caNorm = normalize(caTraces);
        for n = 1:size(caTraces,2)
            caNorm(:,n) = conv(caNorm(:,n), win, 'same');
        end
%         matr = convmtx(win', size(caTraces,1));
%         caNorm = matr * caNorm;
%         caNorm(1:ceil(length(win)/2)-1, :) = [];
%         caNorm(end-floor(length(win)/2)+1:end, :) = [];
        % sort neurons according to correlation with pupil
        valid = ~any(isnan([caNorm ppl]), 2);
        rho = corr(caNorm(valid,:), ppl(valid));
        [~,order] = sort(rho,'descend');
        ax = [0 0 0];
        figure('Position',[1 1 1920 1000])
        subplot(6,1,1)
        hold on
        plot(beeps, probGo, 'Color', [0.5 0 0], 'LineWidth', 2)
        plot(beeps, probCorrect, 'Color', [0 0.5 0], 'LineWidth', 2)
        ylim([0 1])
        ylabel('Probability')
        title(sprintf('%s %s', name, date))
        legend('go choice','correct choice')
        ax(1) = gca;
        subplot(6,1,2)
        plot(caTime, normalize(ppl), 'k', 'LineWidth', 1)
        set(gca, 'box', 'off')
        ylabel('Pupil diameter')
        ax(2) = gca;
        subplot(6,1,3:6)
        imagesc(caTime([1 end]), [1 size(caTraces,2)], caNorm(:,order)', [-3 7])
        colormap(gca, flip(gray(500)))
        set(gca, 'box', 'off')
        xlabel('Time (s)')
        ylabel('#Neurons')
        ax(3) = gca;
        linkaxes(ax, 'x')
        axis tight
        
        savefig(fullfile(folderPlots, 'populationRasters', ...
            sprintf('%s_%s', name, date)))
        saveas(gcf, fullfile(folderPlots, 'populationRasters', ...
            sprintf('%s_%s.png', name, date)))
        close gcf
    end
end

%% Collect data
timePeriods = {'pre-stim', 'stim', 'feedback'};

% collect predictors and responses per dataset
contrast = {}; % [left right]
choice = {}; % -1 left, 0 nogo, 1 right
outcome = {}; % 0 incorrect, 1 correct
pupilSize = {}; % [trials x 3]; columns: pre-stim, stim, feedback
responses = {}; % [trials x 3 x neurons]; columns: pre-stim, stim, feedback
subjects = {};
dates = {};

earlyMove = {}; % 1 if wheel moved before beep
wheelVelThresh = 30;

subjDirs = dir(fullfile(folderData, 'SS*'));
k = 1;
for subj = 1:length(subjDirs)
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folderData, name, '2*'));
    for dt = 1:length(dateDirs)
        date = dateDirs(dt).name;
        if ~isfile(fullfile(folderData, name, date, 'eye.diameter.npy'))
            continue
        end
        subjects{k} = name;
        dates{k} = date;
        
        % load data
        caTraces = readNPY(fullfile(folderData, name, date, '_ss_2pCalcium.dff.npy'));
        caTime = readNPY(fullfile(folderData, name, date, '_ss_2pCalcium.timestamps.npy'));
        planes = readNPY(fullfile(folderData, name, date, '_ss_2pRois._ss_2pPlanes.npy'));
        planeDelays = readNPY(fullfile(folderData, name, date, '_ss_2pPlanes.delay.npy'));
        pupil = readNPY(fullfile(folderData, name, date, 'eye.diameter.npy'));
        pupilTime = readNPY(fullfile(folderData, name, date, 'eye.timestamps.npy'));
        stimT = readNPY(fullfile(folderData, name, date, '_ibl_trials.stim_intervals.npy'));
        contrastL = readNPY(fullfile(folderData, name, date, '_ibl_trials.contrastLeft.npy'));
        contrastR = readNPY(fullfile(folderData, name, date, '_ibl_trials.contrastRight.npy'));
        wheelVel = readNPY(fullfile(folderData, name, date, '_ibl_wheel.velocity.npy'));
        wheelTime = readNPY(fullfile(folderData, name, date, '_ibl_wheel.timestamps.npy'));
        cueInteractiveDelays = readNPY(fullfile(folderData, name, date, '_ibl_trials.cueInteractiveDelay.npy'));
        ch = readNPY(fullfile(folderData, name, date, '_ibl_trials.choice.npy'));
        fb = readNPY(fullfile(folderData, name, date, '_ibl_trials.feedbackType.npy'));
        fbTimes = readNPY(fullfile(folderData, name, date, '_ibl_trials.feedback_times.npy'));
        
        numNeurons = size(caTraces,2);
        numTrials = length(stimT);
        
        % collect task events
        contrast{k} = [contrastL contrastR];
        ch(ch == 1) = -1;
        ch(ch == 2) = 1;
        ch(ch == 3) = 0;
        choice{k} = ch;
        outcome{k} = fb;
        
        % determine early moves (before earliest occurrence of beep)
        minDelay = min(cueInteractiveDelays(:,1));
        wheelAlgn = traces.getAlignedTraces(wheelVel, wheelTime, ...
            stimT(:,1), [0 minDelay]);
        earlyMove{k} = max(abs(wheelAlgn), [], 1)' >= wheelVelThresh;
        
        % collect pupil and neural responses for different time periods
        events = [stimT(:,1) stimT(:,1) fbTimes];
        windows = [-0.6 0; 0 minDelay; 0 1];
        pupilSize{k} = [];
        responses{k} = [];
        for ev = 1:size(events,2)
            ppl = traces.getAlignedTraces(pupil, pupilTime, ...
                events(:,ev), windows(ev,:));
            pupilSize{k} = [pupilSize{k}, nanmean(ppl,1)'];
            
            ca = [];
            planesUni = unique(planes)';
            for p = 1:length(planesUni)
                ind = planes == planesUni(p);
                caP = traces.getAlignedTraces(caTraces(:,ind), ...
                    caTime+planeDelays(p), events(:,ev), windows(ev,:));
                if isempty(ca) && sum(ind)>0
                    ca = NaN(size(caP,1), size(caP,2), numNeurons);
                end
                ca(:,:,ind) = caP;
            end
            responses{k} = [responses{k}, permute(nanmean(ca,1), [2 1 3])];
        end
        
        k = k + 1;
    end
end

%% Plot binned responses, pupil, choice and outcome across time for each data

for k = 1:length(responses)
    resp = responses{k};
    resp = resp ./ nanstd(resp(:,1,:));
    ppl = pupilSize{k};
    for period = 1:size(resp,2)
        % probability of go choice over time
        probGo = abs(choice{k}) > 0;
        probGo = smooth(probGo, 9);
        % probability of correct choice
        probCorrect = smooth(outcome{k}, 9);
        % find correlation between pupil and responses
        numTrials = size(resp,1);
        r = reshape(resp, numTrials, []);
        p = reshape(ppl, numTrials, []);
        valid = ~any(isnan([r p]),2);
        rho = corr(squeeze(resp(valid,1,:)), ppl(valid,1));
        [~,order] = sort(rho, 'descend');
        % plot
        ax = [0 0 0];
        figure('Position',[1 1 1920 1000])
        subplot(6,1,1)
        hold on
        plot(probGo, 'Color', [0.5 0 0])
        plot(probCorrect, 'Color', [0 0.5 0])
        ylim([0 1])
        ylabel('Probability')
        title(sprintf('%s %s: %s', subjects{k}, dates{k}, timePeriods{period}))
        legend('Go choices','Correct choices')
        ax(1) = gca;
        subplot(6,1,2)
        hold on
        plot((1:numTrials)'+[-0.3 0 0.3], ppl, 'LineWidth', 1)
        ylabel('Pupil diameter')
        ax(2) = gca;
        subplot(6,1,3:6)
        r = resp(:,period,order);
        mini = min(r(:));
        maxi = max(r(:));
        rng = maxi - mini;
        imagesc([1 numTrials], [1 size(resp,3)], ...
            squeeze(resp(:,period,order))', [mini+.1*rng maxi-.3*rng])
        colormap(gca, flip(gray(500)))
        set(gca, 'box', 'off')
        xlabel('Time (s)')
        ylabel('#Neurons')
        ax(3) = gca;
        linkaxes(ax, 'x')
        axis tight
        
        savefig(fullfile(folderPlots, 'populationRasters', ...
            sprintf('%s_%s_%s', subjects{k}, dates{k}, timePeriods{period})))
        saveas(gcf, fullfile(folderPlots, 'populationRasters', ...
            sprintf('%s_%s_%s.png', subjects{k}, dates{k}, timePeriods{period})))
        close gcf
    end
end