%% Folders
folderData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\NPY\task_2p';
folderTools = 'C:\STORAGE\workspaces';
folderCode = 'C:\dev\workspace\CortexLab';
folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\wheelTask';
folderPlots = fullfile(folderResults, 'Plots');

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(fullfile(folderCode))

%% Collect data
timePeriods = {'pre-stim', 'stim', 'feedback'};

% collect predictors and responses per dataset
contrast = {}; % [left right]
choice = {}; % -1 left, 0 nogo, 1 right
outcome = {}; % 0 incorrect, 1 correct
pupilSize = {}; % [trials x 3]; columns: pre-stim, stim, feedback
responses = {}; % [trials x 3 x neurons]; columns: pre-stim, stim, feedback

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

%% Fit linear models
% pre-stim: abs(choice), outcome, pupil[, abs(choice(t-1)), outcome(t-1)]
% stim: stimulus, choice, abs(choice), outcome, pupil
% feedback: stimulus choice, abs(choice), outcome, pupil
modelNames = {{'full','no action','no outcome','no pupil','no prev action','no prev outcome'}, ...
    {'full','no stim','no action','no outcome','no pupil','no prev action','no prev outcome'}, ...
    {'full','no stim','no action','no outcome','no pupil','no prev action','no prev outcome'}};
lambda = 10.^(-5:0.5:1);

MSEs = cell(length(responses),size(responses{1},2)); % [neuron x model x lambda]
for k = 1:length(responses)
    fprintf('dataset %d\n', k)
    cUni = unique(contrast{k}, 'rows');
    stimulus = zeros(size(contrast{k},1), size(cUni,1));
    for c = 1:size(cUni)
        stimulus(all(contrast{k} == cUni(c,:),2), c) = 1;
    end
    pplSz = normalize(pupilSize{k});
    for period = 2:size(responses{k},2)
        switch period
            case 1 % pre-stim
                X{1} = [abs(choice{k}) outcome{k} abs(outcome{k}-1) ...
                    pplSz(:,period) abs([0; choice{k}(1:end-1)]) ...
                    [0; outcome{k}(1:end-1)] [0; abs(outcome{k}(1:end-1)-1)]];
                X{2} = [outcome{k} abs(outcome{k}-1) ...
                    pplSz(:,period) abs([0; choice{k}(1:end-1)]) ...
                    [0; outcome{k}(1:end-1)] [0; abs(outcome{k}(1:end-1)-1)]];
                X{3} = [abs(choice{k}) ...
                    pplSz(:,period) abs([0; choice{k}(1:end-1)]) ...
                    [0; outcome{k}(1:end-1)] [0; abs(outcome{k}(1:end-1)-1)]];
                X{4} = [abs(choice{k}) outcome{k} abs(outcome{k}-1) ...
                    abs([0; choice{k}(1:end-1)]) ...
                    [0; outcome{k}(1:end-1)] [0; abs(outcome{k}(1:end-1)-1)]];
                X{5} = [abs(choice{k}) outcome{k} abs(outcome{k}-1) ...
                    pplSz(:,period) ...
                    [0; outcome{k}(1:end-1)] [0; abs(outcome{k}(1:end-1)-1)]];
                X{6} = [abs(choice{k}) outcome{k} abs(outcome{k}-1) ...
                    pplSz(:,period) abs([0; choice{k}(1:end-1)])];
            case {2, 3} % stim, feedback
%                 X{1} = [stimulus choice{k} abs(choice{k}) ...
%                     outcome{k} abs(outcome{k}-1) pplSz(:,period)];
%                 X{2} = [choice{k} abs(choice{k}) ...
%                     outcome{k} abs(outcome{k}-1) pplSz(:,period)];
%                 X{3} = [stimulus abs(choice{k}) ...
%                     outcome{k} abs(outcome{k}-1) pplSz(:,period)];
%                 X{4} = [stimulus choice{k} ...
%                     outcome{k} abs(outcome{k}-1) pplSz(:,period)];
%                 X{5} = [stimulus choice{k} abs(choice{k}) ...
%                     pupilSize{k}(:,period)];
%                 X{6} = [stimulus choice{k} abs(choice{k}) ...
%                     outcome{k} abs(outcome{k}-1)];
                X{1} = [stimulus abs(choice{k}) ...
                    outcome{k} abs(outcome{k}-1) pplSz(:,period) ...
                    abs([0; choice{k}(1:end-1)]) ...
                    [0; outcome{k}(1:end-1)] [0; abs(outcome{k}(1:end-1)-1)]];
                X{2} = [abs(choice{k}) ...
                    outcome{k} abs(outcome{k}-1) pplSz(:,period) ...
                    abs([0; choice{k}(1:end-1)]) ...
                    [0; outcome{k}(1:end-1)] [0; abs(outcome{k}(1:end-1)-1)]];
                X{3} = [stimulus ...
                    outcome{k} abs(outcome{k}-1) pplSz(:,period) ...
                    abs([0; choice{k}(1:end-1)]) ...
                    [0; outcome{k}(1:end-1)] [0; abs(outcome{k}(1:end-1)-1)]];
                X{4} = [stimulus abs(choice{k}) ...
                    pplSz(:,period) ...
                    abs([0; choice{k}(1:end-1)]) ...
                    [0; outcome{k}(1:end-1)] [0; abs(outcome{k}(1:end-1)-1)]];
                X{5} = [stimulus abs(choice{k}) ...
                    outcome{k} abs(outcome{k}-1) ...
                    abs([0; choice{k}(1:end-1)]) ...
                    [0; outcome{k}(1:end-1)] [0; abs(outcome{k}(1:end-1)-1)]];
                X{6} = [stimulus abs(choice{k}) ...
                    outcome{k} abs(outcome{k}-1) pplSz(:,period) ...
                    [0; outcome{k}(1:end-1)] [0; abs(outcome{k}(1:end-1)-1)]];
                X{7} = [stimulus abs(choice{k}) ...
                    outcome{k} abs(outcome{k}-1) pplSz(:,period) ...
                    abs([0; choice{k}(1:end-1)])];
        end
        MSEs{k,period} = NaN(size(responses{k},3), length(X), length(lambda));
        valid = ~isnan(pplSz(:,period));
        if period == 2 % during stim
            valid = valid & ~earlyMove{k};
        end
        for m = 1:length(X)
            for n = 1:size(responses{k},3)
                respN = normalize(responses{k}(:,period,n));
                v = valid & ~isnan(respN);
                CVMdl = fitrlinear(X{m}(v,:), respN(v), ...
                    'Learner', 'leastsquares', 'Solver', 'bfgs', ...
                    'Regularization', 'ridge', 'Lambda', lambda, ...
                    'KFold', 10);
%                 CVMdl = fitrlinear(X{m}(v,:), respN(v), ...
%                     'Learner', 'leastsquares', 'Solver', 'asgd', ...
%                     'Regularization', 'ridge', 'Lambda', lambda, ...
%                     'KFold', 10);
                MSEs{k,period}(n,m,:) = kfoldLoss(CVMdl);
            end
        end
        save(fullfile(folderResults, 'predictResponses', 'mses.mat'), 'MSEs')
    end
end

%% Plot results
folder = fullfile(folderPlots, 'predictNeuralResponses');
for period = 2:size(MSEs,2)
    mse = cat(1, MSEs{:,period}); % [neuron x model x lambda]
    [~,bestLam] = min(mean(mse,2),[],3);
    mseLam = NaN(size(mse,1), size(mse,2)); % [neuron x model]
    for n = 1:size(mse,1)
        mseLam(n,:) = mse(n,:,bestLam(n));
    end
    explVar = 1 - mseLam;
    for m = 2:size(mse,2)
        difference = explVar(:,1) - explVar(:,m);
        p = signrank(difference);
        
        mini = min([explVar(:,1); explVar(:,m)]);
        maxi = max([explVar(:,1); explVar(:,m)]);
        figure
        hold on
        plot([mini maxi],[mini maxi], 'r')
        scatter(explVar(:,1), explVar(:,m), [], 'k', 'filled', ...
            'MarkerFaceAlpha', .2)
        axis([mini maxi mini maxi])
        axis square
        xlabel(sprintf('Expl. var.: %s', modelNames{period}{1}))
        ylabel(sprintf('Expl. var.: %s', modelNames{period}{m}))
        title(sprintf('%s: full - incomplete = %.4f, p=%.4f', timePeriods{period}, ...
            median(difference), p))
        
        savefig(fullfile(folder, sprintf('explVar_%s_full_vs_%s', ...
            timePeriods{period}, modelNames{period}{m})))
        saveas(gcf, fullfile(folder, sprintf('explVar_%s_full_vs_%s.png', ...
            timePeriods{period}, modelNames{period}{m})))
        close gcf
    end
end

%% Plot residuals against predictor
% see whether predictor has a nonlinear relation to 