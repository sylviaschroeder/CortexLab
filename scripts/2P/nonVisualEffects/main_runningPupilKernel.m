%% Parameters
neuralUnits = 'boutons';
% Smoothing running and neural data
sigmaXCorr = 0.1;
sigma = 1; % in sec
% lags of cross-correlation
lags = [-10 10];
% Length of kernel
krnlLims = [-0.5 0.5];

%% Folders
folderBase = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\DataToPublish\arousal_NYP_matlab';
folderTools = 'C:\STORAGE\workspaces';
folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS';

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath('C:\dev\workspace\CortexLab')

%% Plot cross-correlograms between neural activity and running
subjDirs = dir(fullfile(folderBase, neuralUnits, 'SS*'));
fRes = fullfile(folderResults, 'boutons', 'nonVisualEffects', ...
    'plots_running', 'crossCorrelations', sprintf('sigma %.1f', sigmaXCorr));
if ~isfolder(fRes)
    mkdir(fRes)
end
allXCorrs = {};
allLags = {};
count = 1;
for subj = 1:length(subjDirs)
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folderBase, neuralUnits, name, '2*'));
    for dt = 1:length(dateDirs)
        date = dateDirs(dt).name;
        folder = fullfile(folderBase, neuralUnits, name, date, '001');
        
        darkData = io.getDarknessData(folder);
        if isempty(darkData.interval)
            continue
        end
        caData = io.getCalciumData(folder);
        runData = io.getRunningData(folder);
        rhos = readNPY(fullfile(folder, '_ss_corrsRunning.rhosDark.npy'));
        nullRhos = readNPY(fullfile(folder, '_ss_corrsRunning.nullRhosDark.npy'));
        
        significantPos = find(sum(rhos > nullRhos, 2) / size(nullRhos,2) > 0.975);
        significantNeg = find(sum(rhos < nullRhos, 2) / size(nullRhos,2) > 0.975);
        if isempty([significantPos; significantNeg])
            continue
        end
        
        indStart = find(caData.time > darkData.interval(1), 1);
        indEnd = find(caData.time < darkData.interval(2), 1, 'last');
        t = caData.time(indStart : indEnd);
        running = traces.downSample(runData.running, runData.time, ...
            sigmaXCorr, t);
        running = zscore(running);
        tr = caData.traces(indStart : indEnd, :);
        invalid = all(isnan(tr),2);
        running(invalid) = [];
        tr(invalid,:) = [];
        t(invalid) = [];
        
        trFilt = NaN(size(tr));
        for n = 1:size(tr,2)
            trFilt(:,n) = traces.downSample(tr(:,n), t, sigmaXCorr, t, false);
        end
        trNorm = (trFilt - nanmean(trFilt)) ./ nanstd(trFilt);
        
        bin = median(diff(t));
        lagsBins = floor(lags(1) / bin) : ceil(lags(2) / bin);
        lagsT = lagsBins .* bin;
        runCorr = xcorr(running, lagsBins(end), 'unbiased');
        ccorrsPos = NaN(length(lagsT), length(significantPos));
        acorrsPos = NaN(length(lagsT), length(significantPos));
        for n = 1:length(significantPos)
            ccorrsPos(:,n) = xcorr(trNorm(:,significantPos(n)), running, ...
                lagsBins(end), 'unbiased');
            acorrsPos(:,n) = xcorr(trNorm(:,significantPos(n)), ...
                lagsBins(end), 'unbiased');
        end
        ccorrsNeg = NaN(length(lagsT), length(significantNeg));
        acorrsNeg = NaN(length(lagsT), length(significantNeg));
        for n = 1:length(significantNeg)
            ccorrsNeg(:,n) = xcorr(trNorm(:,significantNeg(n)), running, ...
                lagsBins(end), 'unbiased');
            acorrsNeg(:,n) = xcorr(trNorm(:,significantNeg(n)), ...
                lagsBins(end), 'unbiased');
        end
        
        allXCorrs{count,1} = ccorrsPos;
        allXCorrs{count,2} = ccorrsNeg;
        allLags{count} = lagsT;
        count = count + 1;
        
        ax = zeros(1, 5);
        figure('Position', [680 50 750 950])
        tiledlayout(5, 1)
        nexttile % running autocorrelation
        plot(lagsT, runCorr, 'k')
        xlim(lagsT([1 end]))
        xlabel('Lag (in s)')
        title('Autocorrelation of running')
        ax(1) = gca;
        nexttile % xcorr: pos. corr. neurons with running
        if ~isempty(significantPos)
            plot(lagsT, ccorrsPos)
        end
        if length(significantPos) > 1
            hold on
            plot(lagsT, mean(ccorrsPos,2), 'k', 'LineWidth', 2)
        end
        xlim(lagsT([1 end]))
        xlabel('Shift (in s) of neural activity relative to running')
        title('Cross-corr of positively correlated units')
        ax(2) = gca;
        nexttile % auto-corr: pos. corr. neurons
        if ~isempty(significantPos)
            plot(lagsT, acorrsPos)
        end
        if length(significantPos) > 1
            hold on
            plot(lagsT, mean(acorrsPos,2), 'k', 'LineWidth', 2)
        end
        xlim(lagsT([1 end]))
        xlabel('Lag (in s)')
        title('Auto-corr of positively correlated units')
        ax(3) = gca;
        nexttile % xcorr: neg. corr. neurons with running
        if ~isempty(significantNeg)
            plot(lagsT, ccorrsNeg)
        end
        if length(significantNeg) > 1
            hold on
            plot(lagsT, mean(ccorrsNeg,2), 'k', 'LineWidth', 2)
        end
        xlim(lagsT([1 end]))
        xlabel('Shift (in s) of neural activity relative to running')
        title('Cross-corr of negatively correlated units')
        ax(4) = gca;
        nexttile % auto-corr: neg. corr. neurons
        if ~isempty(significantNeg)
            plot(lagsT, acorrsNeg)
        end
        if length(significantNeg) > 1
            hold on
            plot(lagsT, mean(acorrsNeg,2), 'k', 'LineWidth', 2)
        end
        xlim(lagsT([1 end]))
        xlabel('Lag (in s)')
        title('Auto-corr of negatively correlated units')
        ax(5) = gca;
        linkaxes(ax, 'x')
        sgtitle(sprintf('%s %s', name, date))
        
        saveas(gcf, sprintf('%s %s.png', name, date))
    end
end

%% Fit running kernel
subjDirs = dir(fullfile(folderBase, neuralUnits, 'SS*'));
for subj = 1:length(subjDirs)
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folderBase, neuralUnits, name, '2*'));
    for dt = 1:length(dateDirs)
        date = dateDirs(dt).name;
        folder = fullfile(folderBase, neuralUnits, name, date, '001');
        
        darkData = io.getDarknessData(folder);
        if isempty(darkData.interval)
            continue
        end
        caData = io.getCalciumData(folder);
        runData = io.getRunningData(folder);
        rhos = readNPY(fullfile(folder, '_ss_corrsRunning.rhosDark.npy'));
        nullRhos = readNPY(fullfile(folder, '_ss_corrsRunning.nullRhosDark.npy'));
        
        significantPos = find(sum(rhos > nullRhos, 2) / size(nullRhos,2) > 0.975);
        significantNeg = find(sum(rhos < nullRhos, 2) / size(nullRhos,2) > 0.975);
        if isempty([significantPos; significantNeg])
            continue
        end
        
        indStart = find(caData.time > darkData.interval(1), 1);
        indEnd = find(caData.time < darkData.interval(2), 1, 'last');
        t = caData.time(indStart : indEnd);
        running = traces.downSample(runData.running, runData.time, ...
            sigma, t);
        
        [A, numSamples, winTimes] = krnl.getToeplitz(t, [], [], ...
            {running}, {krnlLims}, false);
        tr = caData.traces(indStart : indEnd, :);
        trFilt = NaN(size(tr));
        for n = 1:size(tr,2)
            trFilt(:,n) = traces.downSample(tr(:,n), t, sigma, t);
        end
        trNorm = (trFilt - nanmean(trFilt)) ./ nanstd(trFilt);
        noNaN = ~any(isnan(trNorm), 1);
        hasNaN = find(~noNaN);
        runKrnls = NaN(numSamples, size(trNorm,2));
        runKrnls(:, noNaN) = A \ trNorm(:, noNaN);
        predictions = NaN(size(trNorm));
        for n = hasNaN
            valid = ~isnan(trNorm(:,n));
            runKrnls(:,n) = A(valid,:) \ trNorm(valid,n);
%             runKrnls(:,n) = ridge(trNorm(valid,n), A(valid,:), 0.05);
            predictions(valid,n) = A(valid,:) * runKrnls(:,n);
        end
        
        figure('Position', [680 130 750 850])
        tiledlayout(2, 1)
%         tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact')
        % plot kernels of positively correlated units
        nexttile
        if ~isempty(significantPos)
            plot(winTimes{1}, runKrnls(:, significantPos))
        end
        ylabel('\DeltaF/F_0')
        title('Positively correlated units')
        xlim(winTimes{1}([1 end]))
        % plot kernels of negatively correlated units
        nexttile
        if ~isempty(significantNeg)
            plot(winTimes{1}, runKrnls(:, significantNeg))
        end
        xlabel('Time (s)')
        ylabel('\DeltaF/F_0')
        title('Negatively correlated units')
        xlim(winTimes{1}([1 end]))
        sgtitle(sprintf('%s %s', name, date))
        
        for n = significantNeg'
            figure
            hold on
            plot(t, trNorm(:,n), 'k')
            plot(t, predictions(:,n))
        end
    end
end