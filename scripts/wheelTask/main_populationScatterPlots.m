%% Folders
folderData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\NPY\task_2p';
folderTools = 'C:\STORAGE\workspaces';
folderCode = 'C:\dev\workspace\CortexLab';
folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\wheelTask';
folderPlots = fullfile(folderResults, 'Plots');

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(fullfile(folderCode))

%% Parameters
contraLateralPerSubject = [-1 -1 -1 1];

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

%% Plot scatters: small vs large pupil, go vs nogo choice, correct vs incorrect choice
optionNames = {{'small pupil', 'large pupil'}, ...
    {'nogo choice', 'go choice'}, ...
    {'incorrect choice', 'correct choice'}, ...
    {'prev. nogo choice', 'prev. go choice'}, ...
    {'prev. incorrect choice', 'prev. correct choice'}};
criteriaNames = {'pupil','choice','outcome','prevChoice','prevOutcome'};

for period = 1:size(responses{1},2)
    allResp = [];
    for k = 1:length(responses)
        pupilThresh = nanmean(pupilSize{k}(:,period));
        criteria = [pupilSize{k}(:,period) > pupilThresh, ...
            abs(choice{k})>0, outcome{k}==1, ...
            [0;abs(choice{k}(1:end-1))>0], [0;outcome{k}(1:end-1)==1]];
        
        r = NaN(size(responses{k},3), 2, size(criteria,2));
        for crit = 1:size(criteria,2)
            if period == 2 % stim
                if contraLateralPerSubject(strcmp(subjects{k}, unique(subjects))) < 0
                    stims = contrast{k}(:,1);
                else
                    stims = contrast{k}(:,2);
                end
                st = setdiff(unique(stims), 0);
                rst = NaN(size(responses{k},3), 2, length(st));
                for s = 1:length(st)
                    rst(:,1,s) = nanmean(responses{k}(...
                        stims==st(s) & criteria(:,crit)==0,period,:),1);
                    rst(:,2,s) = nanmean(responses{k}(...
                        stims==st(s) & criteria(:,crit)==1,period,:),1);
                end
                r(:,:,crit) = nanmean(rst,3);
            else
                r(:,1,crit) = nanmean(responses{k}(criteria(:,crit)==0,period,:),1);
                r(:,2,crit) = nanmean(responses{k}(criteria(:,crit)==1,period,:),1);
            end
        end
        r = r ./ squeeze(nanstd(responses{k}(:,1,:),0,1));
        allResp = [allResp; r];
    end
    
    for crit = 1:size(allResp,3)
        mini = min(reshape(allResp(:,:,crit),[],1));
        maxi = max(reshape(allResp(:,:,crit),[],1));
        figure
        hold on
        scatter(allResp(:,1,crit), allResp(:,2,crit), [], 'k', 'filled', ...
            'MarkerFaceAlpha', 0.2)
        plot([mini maxi], [mini maxi], 'r')
        axis([mini maxi mini maxi])
        axis square
        xlabel(optionNames{crit}{1})
        ylabel(optionNames{crit}{2})
        title(timePeriods{period})
        
        savefig(fullfile(folderPlots, 'populationScatters', ...
            sprintf('%s_%s', timePeriods{period}, criteriaNames{crit})))
        saveas(gcf, fullfile(folderPlots, 'populationScatters', ...
            sprintf('%s_%s.png', timePeriods{period}, criteriaNames{crit})))
        close gcf
    end
end
