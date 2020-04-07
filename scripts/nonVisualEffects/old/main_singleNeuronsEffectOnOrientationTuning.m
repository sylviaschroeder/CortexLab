%% Load database
db_driftingGratings

%% Folders
folderROIData = 'C:\Temp2p';

%% Parameters
% smoothing (low-pass filter) of neural traces and non-visual signal before
% correlation between the two
smoothing = 3; %in sec
filtPoly = 3;

% for removal of slow drift in neural responses (high-pass filter)
window = 400; %in sec
percentile = 5;

%% Loop across datasets

for k=1:length(db)
    
    %% Load data
    folder = [folderROIData filesep db(k).subject filesep ...
        db(k).date filesep num2str(db(k).exp)];
    fileStart = [db(k).date '_' num2str(db(k).exp) '_' ...
        db(k).subject];
    file = [fileStart '_2P_plane%03d_ROI.mat'];
    clear info
    neurons = cell(1, length(db(k).planes));
    neuronIDs = [];
    for iPlane=1:length(db(k).planes)
        % load meta
        data=load(fullfile(folder, sprintf(file,db(k).planes(iPlane))));
        info(iPlane) = orderfields(data.meta);
        ind = find(strcmp(data.meta.ROI.CellClasses, 's'));
        neurons{iPlane} = ind;
        neuronIDs = [neuronIDs; [ones(length(ind),1) * ...
            db(k).planes(iPlane), ind']];
    end
    
    %% Preprocess calcium traces
    % get background amplitude (should be calculated when extracting
    % Ca-traces)
    
    % interpolate neural responses to get a single series of time points
    [calciumTraces, calciumTime] = ssLocal.matchSampleTimesAcrossPlanes(info);
    for iPlane = 1:length(db(k).planes)
        calciumTraces{iPlane} = calciumTraces{iPlane}(:,neurons{iPlane});
    end
    calciumTraces = cat(2, calciumTraces{:});
    % caclulate F_0
    [~,F_0] = ssLocal.removeSlowDrift(calciumTraces, calciumTime, window, percentile);
    calciumTraces = (calciumTraces - F_0) ./ F_0;
    % low-pass filter traces (only for correlation with non-visual data)
    filtWindow = ceil(smoothing / median(diff(calciumTime)));
    if mod(filtWindow,2) == 0
        filtWindow = filtWindow-1;
    end
    calciumTraces_smoothed = sgolayfilt(calciumTraces, filtPoly, filtWindow);
    
    %% Partition non-visual data: high vs. low
    % get optimal threshold for running/not running and large/small pupil
    % ((1) running, (2) pupil diameter)
    [thresholds, r_squared, conditionNames] = nonVis.findBestThreshold(info, 1); % remove slow drift == 1
    
    % classify data
    % (1) running
    [negativeCond, positiveCond, conditionTime] = ...
        nonVis.classifyData(info(1), conditionNames{1}, 1, thresholds(1));
    class_run.negCond = negativeCond;
    class_run.posCond = positiveCond;
    class_run.time = conditionTime;
    % (2) pupil diameter
    [negativeCond, positiveCond, conditionTime] = ...
        nonVis.classifyData(info(1), conditionNames{2}, 1, thresholds(2));
    class_pupil.negCond = negativeCond;
    class_pupil.posCond = positiveCond;
    class_pupil.time = conditionTime;
    
    %% Classify stimulus responses according to non-visual data
    % get stimulus information
    stimTimes = ppbox.getStimTimes(info(1));
    stimSequence = ppbox.getStimSequence(info(1));
    stimMatrix = ppbox.buildStimMatrix(stimSequence, stimTimes, calciumTime);
    
    % calculate tuning curves using integral of raw responses to stimuli
    % (1) running(tuningCurves: [stimulus x neuron x condition])
    tuningCurves_run = ...
        gratings.getConditionedOrientationTuning_raw( ...
        calciumTraces, calciumTime, [class_run.negCond, class_run.posCond], ...
        class_run.time, stimMatrix, stimSequence);
    % (2) pupil diameter
    tuningCurves_pupil = ...
        gratings.getConditionedOrientationTuning_raw( ...
        calciumTraces, calciumTime, [class_pupil.negCond, class_pupil.posCond], ...
        class_pupil.time, stimMatrix, stimSequence);
    
    %% Classify neurons according to correlation with non-visual data
    % find neurons most correlated to running
    ballData = nonVis.getRunningSpeed(data.meta);
    if ~isempty(ballData)
        filtWindow = ceil(smoothing / median(diff(ballData.t)));
        if mod(filtWindow,2) == 0
            filtWindow = filtWindow-1;
        end
        % remove high-frequencies from running speed
        running = sgolayfilt(ballData.total, filtPoly, filtWindow);
        [runResults, figHandles] = nonVis.getCorrToNonVisData(calciumTraces_smoothed, ...
            calciumTime, running, ballData.t, 'Running Speed', ...
            neuronIDs, .3, 1);
        for iFig = 1:length(figHandles.traces)
            savefig(figHandles.traces(iFig), fullfile(folderCorrRunning, ...
                [fileStart sprintf('_traces%02d',iFig)]), 'compact');
        end
        savefig(figHandles.rho, fullfile(folderCorrRunning, ...
            [fileStart '_rhos']), 'compact');
    end
    
    % find neurons most correlated to pupil diameter
    [pupilData, pupilTime] = nonVis.loadPupilData(data.meta);
    if ~isempty(pupilData)
        pupilTime(length(pupilData.x)+1:end) = [];
        filtWindow = ceil(smoothing / median(diff(pupilTime)));
        if mod(filtWindow,2) == 0
            filtWindow = filtWindow-1;
        end
        pupilDiam = sqrt(4 * pupilData.area / pi);
        % remove high-frequencies from running speed
        pupilDiam = sgolayfilt(pupilDiam, filtPoly, filtWindow);
        pupilDiam(pupilData.blink | ~pupilData.goodFit) = NaN;
        [pupilResults, figHandles] = nonVis.getCorrToNonVisData(calciumTraces_smoothed, ...
            calciumTime, pupilDiam, pupilTime, 'Pupil Diam.', ...
            neuronIDs, .3, 1);
        for iFig = 1:length(figHandles.traces)
            savefig(figHandles.traces(iFig), fullfile(folderCorrPupil, ...
                [fileStart sprintf('_traces%02d',iFig)]), 'compact');
        end
        savefig(figHandles.rho, fullfile(folderCorrPupil, ...
            [fileStart '_rhos']), 'compact');
    end
    
    %% Determine effect of non-visual components on neural responses
    if ~isempty(ballData)
        % Plot running data
        % (1) All neurons
        figure
        dataX = reshape(tuningCurves_run(:,:,1),[],1);
        dataY = reshape(tuningCurves_run(:,:,2),[],1);
        ind = ~any(isnan([dataX, dataY]), 2);
        dataX = dataX(ind);
        dataY = dataY(ind);
        plot(dataX, dataY, 'k.')
        mini = min([0;dataX;dataY]);
        maxi = max([dataX;dataY]);
        axis([mini 1.2*maxi mini 1.2*maxi])
        hold on
        % fit ax + b
        [f, gof] = fit(dataX, dataY, 'poly1');
        plot([mini maxi], f([mini maxi]), 'm')
        a = coeffvalues(f);
        confidenceInt = confint(f);
        text(mini+0.6*(1.2*maxi-mini), mini+0.9*(1.2*maxi-mini), ...
            sprintf('y = ax + b\na = %.2f (%.2f %.2f)\nb = %.2f (%.2f %.2f)\nadj R^2 = %.3f', ...
            a(1), confidenceInt(1,1), confidenceInt(2,1), ...
            a(2), confidenceInt(1,2), confidenceInt(2,2), gof.adjrsquare), ...
            'Color', 'm')
        % fit ax
        [f, gof] = fit(dataX, dataY, 'a*x');
        plot([mini maxi], f([mini maxi]), 'r')
        a = coeffvalues(f);
        confidenceInt = confint(f);
        text(mini+0.6*(1.2*maxi-mini), mini+0.12*(1.2*maxi-mini), ...
            sprintf('y = ax\na = %.2f (%.2f %.2f)\nadj R^2 = %.3f', ...
            a, confidenceInt(1), confidenceInt(2), gof.adjrsquare), ...
            'Color', 'r')
        % fit x + b
        [f, gof] = fit(dataX, dataY, 'x+b');
        plot([mini maxi], f([mini maxi]), 'b')
        a = coeffvalues(f);
        confidenceInt = confint(f);
        text(mini+0.03*(1.2*maxi-mini), mini+0.9*(1.2*maxi-mini), ...
            sprintf('y = x + b\nb = %.2f (%.2f %.2f)\nadj R^2 = %.3f', ...
            a, confidenceInt(1), confidenceInt(2,1), gof.adjrsquare), ...
            'Color', 'b')
        % labels
        xlabel('Low running speed')
        ylabel('High running speed')
        title('Avg. response per grating and neuron - all neurons')
        savefig(gcf, fullfile(folderEffectRunning, [fileStart '_allNeurons']), 'compact');
        
        % (2) only well correlated neurons
        figure
        ind = runResults.order(runResults.rho >= 0.3 & runResults.p < 0.05);
        dataX = reshape(tuningCurves_run(:,ind,1),[],1);
        dataY = reshape(tuningCurves_run(:,ind,2),[],1);
        ind = ~any(isnan([dataX, dataY]), 2);
        dataX = dataX(ind);
        dataY = dataY(ind);
        plot(dataX, dataY, 'k.')
        mini = min([0;dataX;dataY]);
        maxi = max([dataX;dataY]);
        axis([mini 1.2*maxi mini 1.2*maxi])
        hold on
        % fit ax + b
        [f, gof] = fit(dataX, dataY, 'poly1');
        plot([mini maxi], f([mini maxi]), 'm')
        a = coeffvalues(f);
        confidenceInt = confint(f);
        text(mini+0.6*(1.2*maxi-mini), mini+0.9*(1.2*maxi-mini), ...
            sprintf('y = ax + b\na = %.2f (%.2f %.2f)\nb = %.2f (%.2f %.2f)\nadj R^2 = %.3f', ...
            a(1), confidenceInt(1,1), confidenceInt(2,1), ...
            a(2), confidenceInt(1,2), confidenceInt(2,2), gof.adjrsquare), ...
            'Color', 'm')
        % fit ax
        [f, gof] = fit(dataX, dataY, 'a*x');
        plot([mini maxi], f([mini maxi]), 'r')
        a = coeffvalues(f);
        confidenceInt = confint(f);
        text(mini+0.6*(1.2*maxi-mini), mini+0.12*(1.2*maxi-mini), ...
            sprintf('y = ax\na = %.2f (%.2f %.2f)\nadj R^2 = %.3f', ...
            a, confidenceInt(1), confidenceInt(2), gof.adjrsquare), ...
            'Color', 'r')
        % fit x + b
        [f, gof] = fit(dataX, dataY, 'x+b');
        plot([mini maxi], f([mini maxi]), 'b')
        a = coeffvalues(f);
        confidenceInt = confint(f);
        text(mini+0.03*(1.2*maxi-mini), mini+0.9*(1.2*maxi-mini), ...
            sprintf('y = x + b\nb = %.2f (%.2f %.2f)\nadj R^2 = %.3f', ...
            a, confidenceInt(1), confidenceInt(2), gof.adjrsquare), ...
            'Color', 'b')
        % labels
        xlabel('Low running speed')
        ylabel('High running speed')
        title('Avg. response per grating and neuron - only well corr. neurons')
        savefig(gcf, fullfile(folderEffectRunning, [fileStart '_wellCorrNeurons']), 'compact');
    end
    
    if ~isempty(pupilData)
        % Plot pupil data
        % (1) All neurons
        figure
        dataX = reshape(tuningCurves_pupil(:,:,1),[],1);
        dataY = reshape(tuningCurves_pupil(:,:,2),[],1);
        ind = ~any(isnan([dataX, dataY]), 2);
        dataX = dataX(ind);
        dataY = dataY(ind);
        plot(dataX, dataY, 'k.')
        mini = min([0;dataX;dataY]);
        maxi = max([dataX;dataY]);
        axis([mini 1.2*maxi mini 1.2*maxi])
        hold on
        % fit ax + b
        [f, gof] = fit(dataX, dataY, 'poly1');
        plot([mini maxi], f([mini maxi]), 'm')
        a = coeffvalues(f);
        confidenceInt = confint(f);
        text(mini+0.6*(1.2*maxi-mini), mini+0.9*(1.2*maxi-mini), ...
            sprintf('y = ax + b\na = %.2f (%.2f %.2f)\nb = %.2f (%.2f %.2f)\nadj R^2 = %.3f', ...
            a(1), confidenceInt(1,1), confidenceInt(2,1), ...
            a(2), confidenceInt(1,2), confidenceInt(2,2), gof.adjrsquare), ...
            'Color', 'm')
        % fit ax
        [f, gof] = fit(dataX, dataY, 'a*x');
        plot([mini maxi], f([mini maxi]), 'r')
        a = coeffvalues(f);
        confidenceInt = confint(f);
        text(mini+0.6*(1.2*maxi-mini), mini+0.12*(1.2*maxi-mini), ...
            sprintf('y = ax\na = %.2f (%.2f %.2f)\nadj R^2 = %.3f', ...
            a, confidenceInt(1), confidenceInt(2), gof.adjrsquare), ...
            'Color', 'r')
        % fit x + b
        [f, gof] = fit(dataX, dataY, 'x+b');
        plot([mini maxi], f([mini maxi]), 'b')
        a = coeffvalues(f);
        confidenceInt = confint(f);
        text(mini+0.03*(1.2*maxi-mini), mini+0.9*(1.2*maxi-mini), ...
            sprintf('y = x + b\nb = %.2f (%.2f %.2f)\nadj R^2 = %.3f', ...
            a, confidenceInt(1), confidenceInt(2,1), gof.adjrsquare), ...
            'Color', 'b')
        % labels
        xlabel('Small pupil')
        ylabel('Large pupil')
        title('Avg. response per grating and neuron - all neurons')
        savefig(gcf, fullfile(folderEffectPupil, [fileStart '_allNeurons']), 'compact');
        
        % (2) only well correlated neurons
        figure
        ind = pupilResults.order(pupilResults.rho >= 0.3 & pupilResults.p < 0.05);
        dataX = reshape(tuningCurves_pupil(:,ind,1),[],1);
        dataY = reshape(tuningCurves_pupil(:,ind,2),[],1);
        ind = ~any(isnan([dataX, dataY]), 2);
        dataX = dataX(ind);
        dataY = dataY(ind);
        plot(dataX, dataY, 'k.')
        mini = min([0;dataX;dataY]);
        maxi = max([dataX;dataY]);
        axis([mini 1.2*maxi mini 1.2*maxi])
        hold on
        % fit ax + b
        [f, gof] = fit(dataX, dataY, 'poly1');
        plot([mini maxi], f([mini maxi]), 'm')
        a = coeffvalues(f);
        confidenceInt = confint(f);
        text(mini+0.6*(1.2*maxi-mini), mini+0.9*(1.2*maxi-mini), ...
            sprintf('y = ax + b\na = %.2f (%.2f %.2f)\nb = %.2f (%.2f %.2f)\nadj R^2 = %.3f', ...
            a(1), confidenceInt(1,1), confidenceInt(2,1), ...
            a(2), confidenceInt(1,2), confidenceInt(2,2), gof.adjrsquare), ...
            'Color', 'm')
        % fit ax
        [f, gof] = fit(dataX, dataY, 'a*x');
        plot([mini maxi], f([mini maxi]), 'r')
        a = coeffvalues(f);
        confidenceInt = confint(f);
        text(mini+0.6*(1.2*maxi-mini), mini+0.12*(1.2*maxi-mini), ...
            sprintf('y = ax\na = %.2f (%.2f %.2f)\nadj R^2 = %.3f', ...
            a, confidenceInt(1), confidenceInt(2), gof.adjrsquare), ...
            'Color', 'r')
        % fit x + b
        [f, gof] = fit(dataX, dataY, 'x+b');
        plot([mini maxi], f([mini maxi]), 'b')
        a = coeffvalues(f);
        confidenceInt = confint(f);
        text(mini+0.03*(1.2*maxi-mini), mini+0.9*(1.2*maxi-mini), ...
            sprintf('y = x + b\nb = %.2f (%.2f %.2f)\nadj R^2 = %.3f', ...
            a, confidenceInt(1), confidenceInt(2), gof.adjrsquare), ...
            'Color', 'b')
        % labels
        xlabel('Small pupil')
        ylabel('Large pupil')
        title('Avg. response per grating and neuron - only well corr. neurons')
        savefig(gcf, fullfile(folderEffectPupil, [fileStart '_wellCorrNeurons']), 'compact');
    end
    
    close all
end