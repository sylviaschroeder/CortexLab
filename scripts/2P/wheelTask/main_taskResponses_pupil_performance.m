%% Folders
folderData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\NPY\task_2p';
folderTools = 'C:\STORAGE\workspaces';
folderCode = 'C:\dev\workspace\CortexLab';
folderPlots = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\wheelTask\Plots';

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(fullfile(folderCode))

%% Plot responses of populations, separated by pupil, choice and performance
% critNames = {'Go + corr','Go + incorr','NoGo + corr','NoGo + incorr'};
critNames = {'Hit','False al.','Corr. rej.','Miss'};
titles = {'beep to feedback', 'after feedback'};

folder = fullfile(folderPlots, 'populationEventResponses');
if ~isfolder(folder)
    mkdir(folder)
end

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
        
        numTrials = length(beeps);
        
        % trial periods of interest
        events = {beeps+.1, fbTimes};
        windows = {[zeros(numTrials,1) fbTimes-beeps-.1], [0 1]};
        
        % define separation criteria
        criteria = {choice<3 & fb>0, choice<3 & fb==0 & contrastL==0 & contrastR==0, ...
            choice==3 & fb>0, choice==3 & fb==0};
        
        for ev = 1:length(events)
            % pupil
            pupilAlgn = traces.getAlignedTraces(pupil, pupilTime, ...
                events{ev}, windows{ev});
            % neural traces
            caAlgnBsln = []; % [t x trials x neurons]
            caAlgnEv = []; % [t x trials x neurons]
            planesUni = unique(planes)';
            for p = 1:length(planesUni)
                ind = planes == planesUni(p);
                ca = traces.getAlignedTraces(caTraces(:,ind), ...
                    caTime+planeDelays(p), stimT, [-0.6 0]);
                if isempty(caAlgnBsln) && sum(ind)>0
                    caAlgnBsln = NaN(size(ca,1), size(ca,2), size(caTraces,2));
                end
                caAlgnBsln(:,:,ind) = ca;
                ca = traces.getAlignedTraces(caTraces(:,ind), ...
                    caTime+planeDelays(p), events{ev}, windows{ev});
                if isempty(caAlgnEv) && sum(ind)>0
                    caAlgnEv = NaN(size(ca,1), size(ca,2), size(caTraces,2));
                end
                caAlgnEv(:,:,ind) = ca;
            end
            
            % pupil threshold
            pupilEv = nanmean(pupilAlgn,1)';
            pupilThresh = nanmedian(pupilEv);
            pupilCond = pupilEv > pupilThresh;
            
            % neural response during baseline (pre stim) averaged across time within trial
            caBsln = squeeze(nanmean(caAlgnBsln,1)); % [trials x neurons];
            % neural response to event averaged across time within trial
            caEv = squeeze(nanmean(caAlgnEv,1)); % [trials x neurons]
            
            % exclude invalid trials
            validTrials = ~isnan(pupilEv) & ~all(isnan(caEv),2);
%             validTrials = ~all(isnan(caEv),2);
            
            means = NaN(size(caTraces,2), length(criteria), 2);
            for crit = 1:length(criteria)
                means(:,crit,1) = nanmean(caEv( ...
                    validTrials & criteria{crit} & ~pupilCond,:),1);
                means(:,crit,2) = nanmean(caEv( ...
                    validTrials & criteria{crit} & pupilCond,:),1);
            end
                
            % normalize neural responses by baseline variation
            stdPre = nanstd(caBsln, 1, 1)';
            meansNorm = means ./ stdPre;
            resp = reshape(permute(meansNorm, [1 3 2]), size(meansNorm,1), []);
            
            % sort neurons by responses averaged across criteria
            [~, order] = sort(nanmean(resp,2), 'descend');
            
            % plot
            figure('Position', [100 42 560 954])
            hold on
            imagesc(resp(order,:))
            plot([2 2;4 4;6 6]'+.5, [0 size(resp,1)]+0.5, 'k')
            colormap(flip(gray))
            set(gca, 'XTick', (1:length(criteria))*2-.5, 'XTickLabel', critNames, 'box', 'off')
            axis ij tight
            ylabel('Neuron')
            title(sprintf('%s %s: %s', name, date, titles{ev}))
            
            savefig(fullfile(folder, sprintf('%s_%s_%s', ...
                name, date, titles{ev})))
            saveas(gcf, fullfile(folder, sprintf('%s_%s_%s.png', ...
                name, date, titles{ev})))
            close gcf
        end
    end
end
