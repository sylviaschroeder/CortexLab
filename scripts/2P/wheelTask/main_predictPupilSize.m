% Model pupil size (before stimulus onset, i.e. baseline) based on
% responses of all recorded neurons.

%% Folders
folderData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\NPY\task_2p';
folderTools = 'C:\STORAGE\workspaces';
folderCode = 'C:\dev\workspace\CortexLab';
folderPlots = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\wheelTask';

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(fullfile(folderCode))

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