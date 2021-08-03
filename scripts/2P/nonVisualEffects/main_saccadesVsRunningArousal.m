%% Parameters
neuralUnits = 'boutons';
% minimum distance between saccades
minSaccDist = 0.2;
% Smoothing running and neural data
sigma = 1; % in sec
% lags of cross-correlation
lags = [-10 10];
% Length of kernel
krnlLims = [-0.5 0.5];

%% Folders
folderBase = 'C:\STORAGE\OneDrive - University of Sussex\Lab\DATA\DataToPublish\arousal_NYP_matlab';
folderTools = 'C:\STORAGE\workspaces';
folderResults = 'C:\STORAGE\OneDrive - University of Sussex\Lab\RESULTS';

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath('C:\dev\workspace\CortexLab')

%% Plot pupil size, eye position and saccade onsets
fSacc = fullfile(folderResults, 'boutons', 'nonVisualEffects', ...
    'plots_saccades', 'saccades-pupilSize-eyePos_gratings');
if ~isfolder(fSacc)
    mkdir(fSacc)
end
subjDirs = dir(fullfile(folderBase, neuralUnits, 'SS*'));
for subj = 1:length(subjDirs)
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folderBase, neuralUnits, name, '2*'));
    for dt = 1:length(dateDirs)
        date = dateDirs(dt).name;
        folder = fullfile(folderBase, neuralUnits, name, date, '001');
        
        if ~isfile(fullfile(folder, '_ss_recordings.gratings_intervals.npy'))
            continue
        end
        gratingPeriod = readNPY(fullfile(folder, '_ss_recordings.gratings_intervals.npy'));
        eyeData = io.getPupilData(folder);
        
        % find times of saccades (and direction, amplitude)
        saccIntervals = ...
            eye.findSaccades(eyeData.pos(:,1), eyeData.pos(:,2), ...
            round(minSaccDist ./ median(diff(eyeData.time))), 1);
        
        % subtract mean from eye positions and filter
        pos = eyeData.pos - nanmean(eyeData.pos, 1);
        pos = medfilt1(pos, 3, [], 1, 'includenan', 'truncate');
        
        ax = [0 0 0];
        figure('Position', [20 450 1890 530])
        tiledlayout(3,1)
        nexttile
        hold on
        plot(eyeData.time, eyeData.pupilSize, 'k')
        plot(eyeData.time(saccIntervals(:,1)), eyeData.pupilSize(saccIntervals(:,1)), 'r*')
        ax(1) = gca;
        ylabel('Pupil size')
        nexttile
        hold on
        plot(eyeData.time, pos(:,1), 'k')
        plot(eyeData.time(saccIntervals(:,1)), pos(saccIntervals(:,1),1), 'r*')
        ax(2) = gca;
        ylabel('Horizontal eye position')
        nexttile
        hold on
        plot(eyeData.time, pos(:,2), 'k')
        plot(eyeData.time(saccIntervals(:,1)), pos(saccIntervals(:,1),2), 'r*')
        ax(3) = gca;
        ylabel('Vertical eye position')
        linkaxes(ax, 'x')
        xlim(gratingPeriod)
        xlabel('Time (in s)')
        
        lim2 = get(ax(2), 'yLim');
        lim3 = get(ax(3), 'yLim');
        set(ax(3), 'yLim', lim3(1)+0.5*diff(lim3) + [-1 1].*(0.5*diff(lim2)))
        
        sgtitle(sprintf('%s %s: gratings', name, date))
        
        saveas(gcf, fullfile(fSacc, sprintf('%s_%s.png', name, date)))
        savefig(gcf, fullfile(fSacc, sprintf('%s_%s', name, date)), 'compact')
        close gcf
    end
end

%% Plot responses to gratings + saccade times
subjDirs = dir(fullfile(folderBase, neuralUnits, 'SS*'));
fRes = fullfile(folderResults, 'boutons', 'nonVisualEffects', ...
    'plots_saccades', 'response_gratingTrials');
if ~isfolder(fRes)
    mkdir(fRes)
end
fSacc = fullfile(folderResults, 'boutons', 'nonVisualEffects', ...
    'plots_saccades', 'saccades_gratingTrials');
if ~isfolder(fSacc)
    mkdir(fSacc)
end
for subj = 1:length(subjDirs)
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folderBase, neuralUnits, name, '2*'));
    for dt = 1:length(dateDirs)
        date = dateDirs(dt).name;
        folder = fullfile(folderBase, neuralUnits, name, date, '001');
        
        gratingData = io.getGratingInfo(folder);
        if isempty(gratingData)
            continue
        end
        caData = io.getCalciumData(folder);
        eyeData = io.getPupilData(folder);
        largePupilTrials = readNPY(fullfile(folder, '_ss_gratingTrials.largePupil.npy'));
        parsSmall = readNPY(fullfile(folder, '_ss_tuning.parametersSmall.npy'));
        nullParsSmall = readNPY(fullfile(folder, '_ss_tuning.nullParametersSmall.npy'));
        parsLarge = readNPY(fullfile(folder, '_ss_tuning.parametersLarge.npy'));
        nullParsLarge = readNPY(fullfile(folder, '_ss_tuning.nullParametersLarge.npy'));
        curvesSmall = readNPY(fullfile(folder, '_ss_tuning.curvesSmall.npy'));
        curvesLarge = readNPY(fullfile(folder, '_ss_tuning.curvesLarge.npy'));
        curves = cat(3, curvesSmall, curvesLarge);
        isSuppr = readNPY(fullfile(folder, '_ss_tuning.isSuppressed.npy'));
        predictions = readNPY(fullfile(folder, '_ss_gratingPredictions.dff.npy'));
        predTime = readNPY(fullfile(folder, '_ss_gratingPredictions.timestamps.npy'));
        
        % find times of saccades (and direction, amplitude)
        saccIntervals = ...
            eye.findSaccades(eyeData.pos(:,1), eyeData.pos(:,2), ...
            round(minSaccDist ./ median(diff(eyeData.time))), 1, true);
        saccades = false(size(eyeData.time));
        saccades(saccIntervals) = true;
        saveas(gcf, fullfile(fSacc, ...
            sprintf('%s_%s.png', name, date)))
        savefig(gcf, fullfile(fSacc, ...
            sprintf('%s_%s', name, date)), 'compact')
        close gcf
        
        % sort units according to effect of arousal on tuning
        n = length(caData.ids);
        numSh = 200;
        pars = cat(3, parsSmall, parsLarge);
        nullPars = cat(4, nullParsSmall, nullParsLarge);
        maxima = NaN(n,2);
        nullMaxima = NaN(n,2,numSh);
        means = NaN(n,2);
        nullMeans = NaN(n,2,numSh);
        for iUnit = 1:n
            for cond = 1:2
                if ~isnan(pars(iUnit,2,cond)) % tuned
                    pd = pars(iUnit,1,cond);
                    means(iUnit,cond) = mean(curves(iUnit,:,cond));
                    oris = mod(pd + [0 90 180], 360);
                    resp = gratings.orituneWrappedConditions(pars(iUnit,:,cond), oris);
                    maxima(iUnit,cond) = resp(1);
                    
                    resp = NaN(numSh, 3);
                    crv = NaN(numSh, 360);
                    for sh = 1:numSh
                        oris = mod(nullPars(iUnit,1,sh,cond) + [0 90 180], 360);
                        resp(sh,:) = gratings.orituneWrappedConditions( ...
                            nullPars(iUnit,:,sh,cond), oris);
                        crv(sh,:) = gratings.orituneWrappedConditions(...
                            nullPars(iUnit,:,sh,cond), 1:360);
                    end
                    nullMaxima(iUnit,cond,:) = resp(:,1);
                    nullMeans(iUnit,cond,:) = mean(crv,2);
                else
                    means(iUnit,cond) = pars(iUnit,1,cond);
                    nullMeans(iUnit,cond,:) = nullPars(iUnit,1,:,cond);
                end
            end
        end
        % invert sign of responses of suppressed cells
        maxima(isSuppr==1,:) = -maxima(isSuppr==1,:);
        nullMaxima(isSuppr==1,:,:) = -nullMaxima(isSuppr==1,:,:);
        % response modulation and significance
        modFun = @(a,b) (b-a)./((abs(a)+abs(b)) ./ 2) .* 100;
        mx = maxima;
        ind = all(isnan(maxima),2);
        mx(ind,:) = means(ind,:);
        nmx = nullMaxima;
        nmx(ind,:,:) = nullMeans(ind,:,:);
        respMod = modFun(mx(:,1), mx(:,2));
        nullMod = modFun(squeeze(nmx(:,1,:)), squeeze(nmx(:,2,:)));
        confInt = prctile(nullMod, [2.5 97.5], 2);
        sgnfcnt = respMod < confInt(:,1) | respMod > confInt(:,2);
        [~,sorted] = sort(abs(respMod), 'descend');
        sorted(~sgnfcnt(sorted)) = [];
        
        % get responses, predictions and saccade times aligned to stimulus onsets
        [caAligned, timeCaAligned] = traces.getAlignedTraces(caData.traces, ...
            caData.time, gratingData.times(:,1), [-3 5]);
        [predAligned, timePredAligned] = traces.getAlignedTraces(predictions, ...
            predTime, gratingData.times(:,1), [-3 5]);
        [saccAligned, timeSaccAligned] = traces.getAlignedTraces(saccades, ...
            eyeData.time, gratingData.times(:,1), [-3 5]);
        
        % get matrix: trial indices [directions x repetitions]
        trials = zeros(size(largePupilTrials));
        for st = 1:size(largePupilTrials,1)
            trials(st,:) = find(gratingData.IDs == st);
        end
        
        % make plots
        stimDur = median(diff(gratingData.times,1,2));
        fResSession = fullfile(fRes, sprintf('%s_%s', name, date));
        if ~isfolder(fResSession)
            mkdir(fResSession)
        end
        for iUnit = 1:min(20, length(sorted))
            unit = sorted(iUnit);
            [~,prefStim] = min(abs(pars(unit,1,cond) - ...
                gratingData.directions(1:size(largePupilTrials,1))));
            stimResps = caAligned(:, trials(prefStim,:), unit);
            stimPreds = predAligned(:, trials(prefStim,:), unit);
            mini = min(stimResps(:));
            maxi = max(stimResps(:));
            yDist = maxi - mini;
            small = find(~largePupilTrials(prefStim,:));
            large = find(largePupilTrials(prefStim,:));
            mini = -(max(length(small), length(large))-1) * yDist + ...
                mini - 0.1 * yDist;
            maxi = maxi + 0.1 * yDist;
            
            ax = [0 0];
            figure('Position', [680 90 1180 890])
            tiledlayout(1,2)
            nexttile % small pupil
            hold on
            plot([0 0], [mini maxi], 'k')
            plot([1 1] .* stimDur, [mini maxi], 'k')
            for tr = 1:length(small)
                plot(timeCaAligned, stimResps(:, small(tr)) - (tr-1) * yDist, ...
                    'k', 'LineWidth', 2)
                plot(timePredAligned, stimPreds(:, small(tr)) - (tr-1) * yDist, ...
                    'Color', [1 1 1].*0.5, 'LineWidth', 2)
                indSacc = saccAligned(:, trials(prefStim, small(tr))) > 0;
                plot(timeSaccAligned(indSacc), ones(sum(indSacc),1) .* ...
                    (-(tr-1) * yDist), '^r', 'MarkerFaceColor', 'r')
            end
            ax(1) = gca;
            xlabel('Time from stimulus onset (in s)')
            title('Small pupil')
            nexttile % large pupil
            hold on
            plot([0 0], [mini maxi], 'k')
            plot([1 1] .* stimDur, [mini maxi], 'k')
            for tr = 1:length(large)
                plot(timeCaAligned, stimResps(:, large(tr)) - (tr-1) * yDist, ...
                    'k', 'LineWidth', 2)
                plot(timePredAligned, stimPreds(:, large(tr)) - (tr-1) * yDist, ...
                    'Color', [1 1 1].*0.5, 'LineWidth', 2)
                indSacc = saccAligned(:, trials(prefStim, large(tr))) > 0;
                plot(timeSaccAligned(indSacc), ones(sum(indSacc),1) .* ...
                    (-(tr-1) * yDist), '^r', 'MarkerFaceColor', 'r')
            end
            ax(2) = gca;
            linkaxes(ax, 'xy')
            xlim(timeCaAligned([1 end]))
            ylim([mini maxi])
            xlabel('Time from stimulus onset (in s)')
            title('Large pupil')
            sgtitle(sprintf('Unit %d: Responses to pref. dir. (%d°), resp. mod.: %.2f%%', ...
                unit, gratingData.directions(prefStim), respMod(unit)))
            
            saveas(gcf, fullfile(fResSession, ...
                sprintf('%03d_unit%03d.png', iUnit, unit)))
            savefig(gcf, fullfile(fResSession, ...
                sprintf('%03d_unit%03d', iUnit, unit)), 'compact')
            close gcf
        end
    end
end

%% Calculate response modulations on all trials vs only trials without saccades
responseModulations = [];
fRes = fullfile(folderResults, 'boutons', 'nonVisualEffects', ...
    'plots_saccades');
if ~isfolder(fRes)
    mkdir(fRes)
end

subjDirs = dir(fullfile(folderBase, neuralUnits, 'SS*'));
for subj = 1:length(subjDirs)
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folderBase, neuralUnits, name, '2*'));
    for dt = 1:length(dateDirs)
        date = dateDirs(dt).name;
        folder = fullfile(folderBase, neuralUnits, name, date, '001');
        
        gratingData = io.getGratingInfo(folder);
        if isempty(gratingData)
            continue
        end
        caData = io.getCalciumData(folder);
        eyeData = io.getPupilData(folder);
        parsSmall = readNPY(fullfile(folder, '_ss_tuning.parametersSmall.npy'));
        largePupilTrials = readNPY(fullfile(folder, '_ss_gratingTrials.largePupil.npy'));
        isSuppr = readNPY(fullfile(folder, '_ss_tuning.isSuppressed.npy'));
        
        % find times of saccades (and direction, amplitude)
        saccIntervals = eye.findSaccades(eyeData.pos(:,1), eyeData.pos(:,2), ...
            round(minSaccDist ./ median(diff(eyeData.time))), 1);
        saccades = false(size(eyeData.time));
        for s = 1:length(saccIntervals)
            saccades(saccIntervals(s,1) : saccIntervals(s,2)) = true;
        end
        
        % get responses, predictions and saccade times aligned to stimulus onsets
        stimDur = median(diff(gratingData.times,1,2));
        [caAligned, timeCaAligned] = traces.getAlignedTraces(caData.traces, ...
            caData.time, gratingData.times(:,1), [-1 stimDur+1]); % [t x trials x units]
        [saccAligned, timeSaccAligned] = traces.getAlignedTraces(saccades, ...
            eyeData.time, gratingData.times(:,1), [-1 stimDur+1]); % [t x trials]
        
        % get matrix: trial indices [directions x repetitions]
        trials = zeros(size(largePupilTrials));
        for st = 1:size(largePupilTrials,1)
            trials(st,:) = find(gratingData.IDs == st);
        end
        
        % get response modulation at pref. dir. for each unit
        indBaseline = timeCaAligned < 0;
        indStim = timeCaAligned > 0 & timeCaAligned < stimDur;
        baselines = squeeze(mean(caAligned(indBaseline, :, :), 1)); % [trials x units]
        stimResps = squeeze(mean(caAligned(indStim, :, :), 1)); % [trials x units]
        trialSacc = squeeze(any(saccAligned > 0, 1)); % [trials]
        trialSacc = trialSacc(trials); % [stim x rep]
        
        trialResps = reshape(stimResps(trials,:) - baselines(trials,:), ...
            size(trials,1), size(trials,2), []); % [stim x rep x unit]
        indLarge = repmat(largePupilTrials, 1, 1, size(trialResps,3));
        smallResps = trialResps;
        smallResps(indLarge) = NaN;
        smallResps = squeeze(nanmean(smallResps, 2)); % [stim x unit]
        largeResps = trialResps;
        largeResps(~indLarge) = NaN;
        largeResps = squeeze(nanmean(largeResps, 2)); % [stim x unit]
        
        indSacc = repmat(trialSacc, 1, 1, size(trialResps,3));
        smallResps_noSacc = trialResps;
        smallResps_noSacc(indLarge | indSacc) = NaN;
        smallResps_noSacc = squeeze(nanmean(smallResps_noSacc, 2)); % [stim x unit]
        largeResps_noSacc = trialResps;
        largeResps_noSacc(~indLarge | indSacc) = NaN;
        largeResps_noSacc = squeeze(nanmean(largeResps_noSacc, 2)); % [stim x unit]
        
        respMod = NaN(size(caData.traces,2), 1);
        respMod_noSacc = NaN(size(caData.traces,2), 1);
        modFun = @(a,b) (b-a)./((abs(a)+abs(b)) ./ 2) .* 100;
        for iUnit = 1:size(caData.traces,2)
            if isnan(parsSmall(iUnit,1))
                continue
            end
            [~,prefStim] = min(abs(parsSmall(iUnit,1) - gratingData.directions));
            respMod(iUnit) = modFun(smallResps(prefStim, iUnit), ...
                largeResps(prefStim, iUnit));
            respMod_noSacc(iUnit) = modFun(smallResps_noSacc(prefStim, iUnit), ...
                largeResps_noSacc(prefStim, iUnit));
            if smallResps(prefStim, iUnit) < 0 % suppressed
                respMod(iUnit) = - respMod(iUnit);
                respMod_noSacc(iUnit) = - respMod_noSacc(iUnit);
            end
        end
        responseModulations = [responseModulations; [respMod, respMod_noSacc]];
    end
end

figure
hold on
scatter(responseModulations(:,1), responseModulations(:,2), 'k')
plot([-200 200], [-200 200], 'Color', [1 1 1].*0.5, 'LineWidth', 2)
axis square
xlabel('Resp. mod. - all trials')
ylabel('Resp. mod. - trials without saccades')
saveas(gcf, fullfile(fRes, 'responseModulations_scatter_allTrials-vs-noSaccadeTrials.png'))
savefig(gcf, fullfile(fRes, ...
    'responseModulations_scatter_allTrials-vs-noSaccadeTrials'), 'compact')