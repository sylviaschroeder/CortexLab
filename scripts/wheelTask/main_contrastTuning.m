%% Folders
folderData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\NPY\task_2p';
folderTools = 'C:\STORAGE\workspaces';
folderCode = 'C:\dev\workspace\CortexLab';
folderPlots = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\wheelTask';

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(fullfile(folderCode))

%% Parameters
wheelVelThresh = 30;

%% Plot tuning curves
subjDirs = dir(fullfile(folderData, 'SS*'));
dataset = [];
p_stim = [];
p_pupil = [];
k = 0;
for subj = 1:length(subjDirs)
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folderData, name, '2*'));
    for dt = 1:length(dateDirs)
        date = dateDirs(dt).name;
        if ~isfile(fullfile(folderData, name, date, 'eye.diameter.npy'))
            continue
        end
        k = k + 1;
        
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
        
        % get data aligned to stimulus onset
        periodEnd = min(cueInteractiveDelays(:,1));
        window = [0 periodEnd];
        % pupil
        pupilAlgn = traces.getAlignedTraces(pupil, pupilTime, ...
            stimT, window);
        % wheel
        wheelAlgn = traces.getAlignedTraces(wheelVel, ...
            wheelTime, stimT, window);
        % neural traces
        caAlgnPost = [];
        caAlgnPre = [];
        planesUni = unique(planes)';
        for p = 1:length(planesUni)
            ind = planes == planesUni(p);
            ca = traces.getAlignedTraces(caTraces(:,ind), ...
                caTime+planeDelays(p), stimT, window);
            if isempty(caAlgnPost) && sum(ind)>0
                caAlgnPost = NaN(size(ca,1), size(ca,2), size(caTraces,2));
            end
            caAlgnPost(:,:,ind) = ca;
            ca = traces.getAlignedTraces(caTraces(:,ind), ...
                caTime+planeDelays(p), stimT, [-0.6 0]);
            if isempty(caAlgnPre) && sum(ind)>0
                caAlgnPre = NaN(size(ca,1), size(ca,2), size(caTraces,2));
            end
            caAlgnPre(:,:,ind) = ca;
        end
        
        % pupil threshold
        pupilStim = nanmean(pupilAlgn,1)';
        pupilThresh = nanmedian(pupilStim);
        pupilCond = pupilStim > pupilThresh;
        
        % exclude trials where wheel was moved early or where pupil could
        % not be measured
        validTrials = max(abs(wheelAlgn), [], 1)' < wheelVelThresh & ...
            ~isnan(pupilStim);
        
        % mean neural response to stim
        caStim = squeeze(nanmean(caAlgnPost,1) - nanmean(caAlgnPre,1)); % [trials x neurons]
        
        % generate stim IDs: different IDs for each possible combination of
        % contrasts
        contrasts = [contrastL, contrastR];
        cUni = unique(contrasts, 'rows');
        contrIDs = NaN(length(contrastL),1);
        for c = 1:size(cUni,1)
            contrIDs(all(contrasts == cUni(c,:), 2)) = c;
        end
        
        % plot tuning curves
        folder = fullfile(folderPlots, 'contrastTuning', name, date);
        if ~isfolder(folder)
            mkdir(folder)
        end
        pVals = NaN(size(caTraces,2),2);
        for n = 1:size(caTraces,2)
            tbl = table(caStim(validTrials,n), nominal(contrIDs(validTrials)), ...
                pupilStim(validTrials), 'VariableNames', {'resp','stim','pupil'});
            lme = fitlme(tbl, 'resp ~ stim + pupil');
            stats = anova(lme);
            pVals(n,:) = stats.pValue(2:3);
            
            %get mean and sem for each stimulus and pupil size
            [mResp, semResp, cL, cR] = task.getContrastResponses(caStim(validTrials,n), ...
                contrastL(validTrials), contrastR(validTrials), ...
                pupilCond(validTrials));
            stim = [-flip(cL); cR(2:end)];
            mStim = [flip(squeeze(mResp(:,1,:))); squeeze(mResp(1,2:end,:))];
            semStim = [flip(squeeze(semResp(:,1,:))); squeeze(semResp(1,2:end,:))];
            
            % plot
            figure
            hold on
            errorbar(stim, mStim(:,1), semStim(:,1), 'k', 'CapSize', 0, ...
                'LineWidth', 1)
            errorbar(stim, mStim(:,2), semStim(:,2), 'r', 'CapSize', 0, ...
                'LineWidth', 1)
            xlabel('Stimulus contrast')
            ylabel('\DeltaF/F')
            title(sprintf('P %d, N %d (p_{stim}=%.4f, p_{pupil}=%.4f)', ...
                planes(n), n, stats.pValue(2), stats.pValue(3)))
            legend({'small pupil','large pupil'})
            
            saveas(gcf, fullfile(folder, sprintf('plane%02d_roi%03d.png', ...
                planes(n), n)))
            close(gcf)
        end
        dataset = [dataset; ones(size(pVals,1),1).*k];
        p_stim = [p_stim; pVals(:,1)];
        p_pupil = [p_pupil; pVals(:,2)];
    end
end

figure
pie([sum(p_stim<.05&p_pupil<.05), sum(p_stim<.05&p_pupil>=.05), ...
    sum(p_stim>=.05&p_pupil<.05), sum(p_stim>=.05&p_pupil>=.05)], ...
    {'stim & pupil','stim','pupil','-'})

%% Plot correlations between neural response and pupil before stimulus
window = [-1 0];

subjDirs = dir(fullfile(folderData, 'SS*'));
dataset = [];
rho = [];
p_pupil = [];
k = 0;
for subj = 1:length(subjDirs)
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folderData, name, '2*'));
    for dt = 1:length(dateDirs)
        date = dateDirs(dt).name;
        if ~isfile(fullfile(folderData, name, date, 'eye.diameter.npy'))
            continue
        end
        k = k + 1;
        
        % load data
        caTraces = readNPY(fullfile(folderData, name, date, '_ss_2pCalcium.dff.npy'));
        caTime = readNPY(fullfile(folderData, name, date, '_ss_2pCalcium.timestamps.npy'));
        planes = readNPY(fullfile(folderData, name, date, '_ss_2pRois._ss_2pPlanes.npy'));
        planeDelays = readNPY(fullfile(folderData, name, date, '_ss_2pPlanes.delay.npy'));
        pupil = readNPY(fullfile(folderData, name, date, 'eye.diameter.npy'));
        pupilTime = readNPY(fullfile(folderData, name, date, 'eye.timestamps.npy'));
        stimT = readNPY(fullfile(folderData, name, date, '_ibl_trials.stim_intervals.npy'));
        
        % pupil
        pupilAlgn = traces.getAlignedTraces(pupil, pupilTime, ...
            stimT, window);
        pupilPre = nanmean(pupilAlgn,1)';
        % neural traces
        caAlgnPre = [];
        planesUni = unique(planes)';
        for p = 1:length(planesUni)
            ind = planes == planesUni(p);
            ca = traces.getAlignedTraces(caTraces(:,ind), ...
                caTime+planeDelays(p), stimT, window);
            if isempty(caAlgnPre) && sum(ind)>0
                caAlgnPre = NaN(size(ca,1), size(ca,2), size(caTraces,2));
            end
            caAlgnPre(:,:,ind) = ca;
        end
        caPre = squeeze(nanmean(caAlgnPre,1)); % [trials x neurons]
        
        % plot tuning curves
%         folder = fullfile(folderPlots, 'preStimCorrWPupil', name, date);
%         if ~isfolder(folder)
%             mkdir(folder)
%         end
        r = NaN(size(caTraces,2),1);
        pVals = NaN(size(caTraces,2),1);
        for n = 1:size(caTraces,2)
            ind = ~isnan(caPre(:,n)) & ~isnan(pupilPre);
            [r(n), pVals(n)] = corr(caPre(ind,n), pupilPre(ind));
            
%             % plot
%             figure
%             scatter(pupilPre, caPre(:,n), [], 'k')
%             xlabel('Pupil size (1s pre stim)')
%             ylabel('\DeltaF/F')
%             title(sprintf('P %d, N %d (rho=%.3f, p=%.4f)', ...
%                 planes(n), n, r(n), pVals(n)))
%             
%             saveas(gcf, fullfile(folder, sprintf('plane%02d_roi%03d.png', ...
%                 planes(n), n)))
%             close(gcf)
        end
        dataset = [dataset; ones(size(pVals,1),1).*k];
        rho = [rho; r];
        p_pupil = [p_pupil; pVals];
    end
end

edges = floor(min(rho)*10)/10 : 0.1 : ceil(max(rho)*10)/10;
bins = edges(1:end-1)+0.05;
figure
N1 = histcounts(rho(p_pupil<0.05), edges);
N2 = histcounts(rho(p_pupil>=0.05), edges);
b = bar(bins, [N1; N2], 'stacked');
b(1).FaceColor = 'k';
b(2).FaceColor = 'w';
xlabel('Pearson''s correlation')
ylabel('#Neurons')
legend({'p < 0.05', 'p \geq 0.05'})
title(sprintf('n = %d', sum(~isnan(rho))))

%% Predict pupil baseline from neural population
window = [-1 0];
lambda = 10.^(-4:.5:3);
folder = fullfile(folderPlots, 'predictPupil', 'preStim');
if ~isfolder(folder)
    mkdir(folder)
end

subjDirs = dir(fullfile(folderData, 'SS*'));
dataset = [];
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
        
        % pupil
        pupilAlgn = traces.getAlignedTraces(pupil, pupilTime, ...
            stimT, window);
        pupilPre = nanmean(pupilAlgn,1)';
        % neural traces
        caAlgnPre = [];
        planesUni = unique(planes)';
        for p = 1:length(planesUni)
            ind = planes == planesUni(p);
            ca = traces.getAlignedTraces(caTraces(:,ind), ...
                caTime+planeDelays(p), stimT, window);
            if isempty(caAlgnPre) && sum(ind)>0
                caAlgnPre = NaN(size(ca,1), size(ca,2), size(caTraces,2));
            end
            caAlgnPre(:,:,ind) = ca;
        end
        caPre = squeeze(nanmean(caAlgnPre,1)); % [trials x neurons]
        
        % remove NaNs, and z-score neural responses and pupil
        ind = ~isnan(pupilPre) & ~any(isnan(caPre),2);
        pupilPre = normalize(pupilPre(ind));
        caPre = normalize(caPre(ind,:));
        
        CVMdl = fitrlinear(caPre, pupilPre, 'Lambda', lambda, 'Solver', 'bfgs', ...
            'KFold', 10, 'Learner', 'leastsquares', 'Regularization', 'ridge');
        CVmse = kfoldLoss(CVMdl);
        CVpred = kfoldPredict(CVMdl);
        [~,bestLam] = min(CVmse);
        CVpred = CVpred(:,bestLam);
        residuals = pupilPre - CVpred;
        Mdl = fitrlinear(caPre, pupilPre, 'Lambda', lambda(bestLam), ...
            'Solver', 'bfgs', 'Learner', 'leastsquares', 'Regularization', 'ridge');
        beta = Mdl.Beta;
        
        ax = [0 0];
        figure('Position', [10 170 1230 810])
        % lambda vs MSE
        subplot(3,3,1)
        plot(log10(lambda), CVmse)
        xlabel('Lambda')
        ylabel('MSE (cross-val)')
        % data vs prediction
        subplot(3,3,2)
        mini = min([pupilPre; CVpred]);
        maxi = max([pupilPre; CVpred]);
        hold on
        plot([mini maxi], [mini maxi], 'r')
        scatter(pupilPre, CVpred, [], 'k')
        axis([mini maxi mini maxi])
        axis square
        xlabel('Pupil size (pre-stim)')
        ylabel('Prediction (cross-val)')
        % histogram of coefficients
        subplot(3,3,3)
        histogram(beta, 'FaceColor', 'w');
        xlabel('Fitted coefficients')
        ylabel('# Neurons')
        % residuals across time
        subplot(3,1,2)
        plot(residuals, 'ko')
        xlabel('# Trial')
        ylabel('Residuals')
        set(gca, 'box', 'off')
        ax(1) = gca;
        % pupil size across time
        subplot(3,1,3)
        plot(pupilPre, 'o')
        xlabel('# Trial')
        ylabel('Pupil size (pre-stim)')
        set(gca, 'box', 'off')
        ax(2) = gca;
        
        linkaxes(ax, 'xy')
        xlim([1 length(pupilPre)])
        annotation('textbox', [0 0.97 1 0.03], 'String', ...
            sprintf('%s %s', name, date), 'HorizontalAlignment', 'center', ...
            'LineStyle', 'none', 'FontWeight', 'bold', 'FontSize', 12)
        
        saveas(gcf, fullfile(folder, sprintf('%s_%s.png', name, date)))
        close gcf
    end
end