%% Folders
folderData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\NPY\task_2p';
folderTools = 'C:\STORAGE\workspaces';
folderCode = 'C:\dev\workspace\CortexLab';
folderPlots = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\wheelTask\Plots';

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(fullfile(folderCode))

%% Parameters
wheelVelThresh = 30;

%% Plot tuning curves for single neurons
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
            lm = fitlme(tbl, 'resp ~ stim + pupil');
            stats = anova(lm);
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

%% Plot tuning curves for populations, small vs large pupil, each dataset separately
showPreStimCol = true;
subtractPreStim = false;
critNames = {'pupil','choice','outcome'};
titles = {{'Small pupil (during stim)', 'Large pupil (during stim)'}, ...
    {'Go choice', 'NoGo choice'}, ...
    {'Correct', 'Incorrect'}};

folder = fullfile(folderPlots, 'populationTuningCurves');
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
        choice = readNPY(fullfile(folderData, name, date, '_ibl_trials.choice.npy'));
        fb = readNPY(fullfile(folderData, name, date, '_ibl_trials.feedbackType.npy'));
        
        % get data aligned to stimulus onset
        periodEnd = min(cueInteractiveDelays(:,1));
        window = [0 periodEnd];
        % pupil
        pupilAlgn = traces.getAlignedTraces(pupil, pupilTime, ...
            stimT(:,1), window);
        % wheel
        wheelAlgn = traces.getAlignedTraces(wheelVel, ...
            wheelTime, stimT(:,1), window);
        % neural traces
        caAlgnPost = [];
        caAlgnPre = [];
        planesUni = unique(planes)';
        for p = 1:length(planesUni)
            ind = planes == planesUni(p);
            ca = traces.getAlignedTraces(caTraces(:,ind), ...
                caTime+planeDelays(p), stimT(:,1), window);
            if isempty(caAlgnPost) && sum(ind)>0
                caAlgnPost = NaN(size(ca,1), size(ca,2), size(caTraces,2));
            end
            caAlgnPost(:,:,ind) = ca;
            ca = traces.getAlignedTraces(caTraces(:,ind), ...
                caTime+planeDelays(p), stimT(:,1), [-0.6 0]);
            if isempty(caAlgnPre) && sum(ind)>0
                caAlgnPre = NaN(size(ca,1), size(ca,2), size(caTraces,2));
            end
            caAlgnPre(:,:,ind) = ca;
        end
        
        % pupil threshold
        pupilStim = nanmean(pupilAlgn,1)';
        pupilThresh = nanmedian(pupilStim);
        pupilCond = pupilStim > pupilThresh;
        
        % neural response to stim averaged across time within trial
        caPre = squeeze(nanmean(caAlgnPre,1)); % [trials x neurons]
        if subtractPreStim
            caStim = squeeze(nanmean(caAlgnPost,1) - nanmean(caAlgnPre,1)); % [trials x neurons]
        else
            caStim = squeeze(nanmean(caAlgnPost,1)); % [trials x neurons]
        end
        
        % exclude trials where wheel was moved early or where pupil could
        % not be measured
        validTrials = max(abs(wheelAlgn), [], 1)' < wheelVelThresh & ...
            ~isnan(pupilStim) & ~all(isnan(caStim),2);
        
        % generate stim IDs: different IDs for each possible combination of
        % contrasts
        contrasts = [contrastL, contrastR];
        cUni = unique(contrasts, 'rows');
        [stim, cOrder] = sort(sum(cUni.*[-1 1],2),'ascend');
        cUni = cUni(cOrder,:);
        contrIDs = NaN(length(contrastL),1);
        for c = 1:size(cUni,1)
            contrIDs(all(contrasts == cUni(c,:), 2)) = c;
        end
        
        % define separation criteria
        criteria = {pupilCond, choice<3, fb>0};
        
        for crit = 1:length(criteria)
            % get stimulus responses (mean, sem) and p-values for stim and
            % pupil
            means = NaN(size(caTraces,2), size(cUni,1), 2);
            sems = NaN(size(caTraces,2), size(cUni,1), 2);
            
            for n = 1:size(caTraces,2)
                %get mean and sem for each stimulus and criterium
                [mResp, semResp, cL, cR] = task.getContrastResponses(caStim(validTrials,n), ...
                    contrastL(validTrials), contrastR(validTrials), ...
                    criteria{crit}(validTrials));
                for st = 1:size(cUni,1)
                    means(n,st,:) = mResp(cL==cUni(st,1), cR==cUni(st,2),:);
                    sems(n,st,:) = semResp(cL==cUni(st,1), cR==cUni(st,2),:);
                end
            end
            
            % normalize neural responses by baseline variation
            stdPre = nanstd(caPre, 1, 1)';
            meansNorm = means ./ stdPre;
            preNorm = [nanmean(caPre(~pupilCond,:))' nanmean(caPre(pupilCond,:))'] ./ stdPre;
            
            % sort neurons by maximal response first condition of criterium
            [~, order] = sort(max(meansNorm(:,:,1),[],2), 'descend');
            
            % plot
            if showPreStimCol
                data = [permute(preNorm, [1 3 2]) meansNorm];
                xTicks = [1 (1:2:length(stim))+1];
                xTL = [{'pre'}; num2str(stim(1:2:end))];
            else
                data = meansNorm;
                xTicks = (1:2:length(stim));
                xTL = stim(1:2:end);
            end
            
            mini = min(data(:));
            maxi = max(data(:));
            figure('Position', [100 42 560 954])
            
            subplot(1,2,1)
            hold on
            imagesc(data(order,:,1), [mini maxi])
            if showPreStimCol
                plot([1.5 1.5], [0.5 size(caTraces,2)+0.5], 'k')
            end
            colormap(flip(gray))
            set(gca, 'XTick', xTicks, 'XTickLabel', xTL, 'box', 'off')
            xlabel('Contrast')
            ylabel('Neuron')
            title(titles{crit}{1})
            axis tight ij
            
            subplot(1,2,2)
            hold on
            imagesc(data(order,:,2), [mini maxi])
            if showPreStimCol
                plot([1.5 1.5], [0.5 size(caTraces,2)+0.5], 'k')
            end
            colormap(flip(gray))
            set(gca, 'XTick', xTicks, 'XTickLabel', xTL, 'box', 'off')
            xlabel('Contrast')
            title(titles{crit}{2})
            axis tight ij
            
            annotation('textbox', [0 0.97 1 0.03], 'String', ...
                sprintf('%s %s', name, date), 'HorizontalAlignment', 'center', ...
                'LineStyle', 'none', 'FontWeight', 'bold', 'FontSize', 12)
            
            savefig(fullfile(folder, sprintf('%s_%s_stim_%s', ...
                name, date, critNames{crit})))
            saveas(gcf, fullfile(folder, sprintf('%s_%s_stim_%s.png', ...
                name, date, critNames{crit})))
            close gcf
        end
    end
end

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
