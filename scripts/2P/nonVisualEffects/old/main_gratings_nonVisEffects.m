%% Information on dataset
subject='M150410_SS045';
expDate = '2015-05-04';
exp=1;
planes = 1:4;

info=ppbox.infoPopulate(subject, expDate, exp);

folderROIData = 'C:\Temp2p';
%% Set parameters
method={'separated', 'raw'};
lambda = 0.5;

%% Preprocess neural data


%% Get non-neural data and classify
% find optimal threshold for non-visual data to classify population
% responses
folder = [folderROIData filesep subject filesep expDate filesep num2str(exp)];
file = [info.basename2p '_plane%03d_ROI.mat'];
for iPlane=1:length(planes)
    % load meta
    data=load(fullfile(folder, sprintf(file,planes(iPlane))));
    m = orderfields(data.meta);
    infoROI(iPlane) = m;
end
[thresholds, r_squared, conditionNames] = nonVis.findBestThreshold(infoROI, 1); % remove slow drift == 1
[~, ind] = min(r_squared);
threshold = thresholds(ind);
condition = conditionNames{ind};

% classify data
[negativeCond, positiveCond, conditionTime, conditionNames] = ...
    nonVis.classifyData(info, condition, 1, threshold);

%% Determine tuning for both non-neural conditions (e.g. running and not running)

predictionErrors = [];
fittingParameters = [];
adjustedRsquared = [];
baselineActivity = [];

for iPlane = 1:length(planes)
    % load meta
    filePath = fullfile(info.folderProcessed, sprintf('%s_plane%03d_ROI', info.basename2p, planes(iPlane)));
    load(filePath, 'meta')
    % get raw Ca-traces
    neuronInds = strcmp(meta.ROI.CellClasses, 's');
    tracesRaw = meta.F(:,neuronInds);
    % get stimulus information
    [~, stimSequence, stimMatrix, frameTimes] = ssLocal.getStimulusResponseInfo(meta);
    orientations = gratings.getOrientations(stimSequence);
    
    % calculate tuning curves using integral of raw responses to stimuli
    [tuningCurvesRaw, timeAverageResp, trialConditions] = ...
        gratings.getConditionedOrientationTuning_raw( ...
        tracesRaw, frameTimes, [negativeCond, positiveCond], conditionTime, stimMatrix, stimSequence);
    % fit raw tuning curves
    [fitParsRaw, adjRsqRaw] = gratings.fitTuningCurves(tuningCurvesRaw,orientations(:,1));
    
    % calculate tuning curves using rigde regression and SVD
    [timeCourses, tunings, offsets, baselines, filters, predictions, errors, traces, kernelTimes] = ...
        gratings.getConditionedOrientationTuning_separated(tracesRaw, frameTimes, ...
        [negativeCond, positiveCond], conditionTime, stimMatrix, ...
        stimSequence, lambda);
    % curves from ridge regression
    ind = kernelTimes > 0;
    tuningCurvesRidge = squeeze(mean(filters.cond_regr(:,ind,:), 2));
    tuningCurvesRidge = bsxfun(@plus, tuningCurvesRidge, baselines);
    tunCR_nonCond = squeeze(mean(filters.nonCond_regr(:,ind,:), 2));
    tunCR_nonCond = bsxfun(@plus, tunCR_nonCond, baselines);
    tuningCurvesRidge = cat(3, tuningCurvesRidge(1:size(orientations,1),:), ...
        tuningCurvesRidge(size(orientations,1)+1:end,:), ...
        tunCR_nonCond);
    [fitParsRidge, adjRsqRidge] = gratings.fitTuningCurves(tuningCurvesRidge,orientations(:,1));
    % curves from stimulus-time separated filters
    tuningCurvesSep = bsxfun(@plus, tunings, bsxfun(@plus, offsets, baselines));
    [fitParsSep, adjRsqSep] = gratings.fitTuningCurves(tuningCurvesSep,orientations(:,1));
    
    fields = fieldnames(errors);
    for f = 1:length(fields)
        if ~isfield(predictionErrors, fields{f})
            predictionErrors.(fields{f}) = errors.(fields{f});
        else
            predictionErrors.(fields{f}) = [predictionErrors.(fields{f}), errors.(fields{f})];
        end
    end
    if isempty(fittingParameters)
        fittingParameters.raw = fitParsRaw;
        fittingParameters.ridge = fitParsRidge;
        fittingParameters.sep = fitParsSep;
    else
        fittingParameters.raw = [fittingParameters.raw; fitParsRaw];
        fittingParameters.ridge = [fittingParameters.ridge; fitParsRidge];
        fittingParameters.sep = [fittingParameters.sep; fitParsSep];
    end
    if isempty(adjustedRsquared)
        adjustedRsquared.raw = adjRsqRaw;
        adjustedRsquared.ridge = adjRsqRidge;
        adjustedRsquared.sep = adjRsqSep;
    else
        adjustedRsquared.raw = [adjustedRsquared.raw; adjRsqRaw];
        adjustedRsquared.ridge = [adjustedRsquared.ridge; adjRsqRidge];
        adjustedRsquared.sep = [adjustedRsquared.sep; adjRsqSep];
    end
    baselineActivity = [baselineActivity, baselines];
    
    for neuron = 1:size(tracesRaw,2)
        condColors = 'kr';
        modelOris = 0:5:360;
        
        % Plot predictions using classified trials
        figure('Position', [10 680 1900 420])
        hold on
        h = zeros(1,3);
        h(1) = plot(frameTimes, traces(:,neuron), 'k');
        h(2) = plot(frameTimes, predictions.cond_regr(:,neuron), 'b');
        h(3) = plot(frameTimes, predictions.cond_sep(:,neuron), 'r');
        legend(h, 'raw', sprintf('regression (err.: %.2f)',errors.cond_regr(neuron)), ...
            sprintf('separated (err.: %.2f)',errors.cond_sep(neuron)))
        title('Predictions under classification of non-visual data')
        xlabel('Time (in s)')
        ylabel('Calcium response')
        
        % Plot predictions without classifying data according to non-visual data
        figure('Position', [10 680 1900 420])
        hold on
        h = zeros(1,3);
        h(1) = plot(frameTimes, traces(:,neuron), 'k');
        h(2) = plot(frameTimes, predictions.nonCond_regr(:,neuron), 'b');
        h(3) = plot(frameTimes, predictions.nonCond_sep(:,neuron), 'r');
        legend(h, 'raw', sprintf('regression (err.: %.2f)',errors.nonCond_regr(neuron)), ...
            sprintf('separated (err.: %.2f)',errors.nonCond_sep(neuron)));
        title('Predictions without classification of non-visual data')
        xlabel('Time (in s)')
        ylabel('Calcium response')
        
        % Plot tuning curves from integrated raw responses
        maxi = 0;
        figure('Position', [10 680 1900 420])
        subplot(1,3,1)
        hold on
        h = [0 0];
        for cond = 1:2
            resp = squeeze(timeAverageResp(:,neuron,:));
            resp(squeeze(trialConditions(cond,:,:))==0) = NaN;
            plot(orientations(:,1), resp, 'o', 'Color', condColors(cond))
            modelFit = oritune(squeeze(fitParsRaw(neuron,cond,:)), modelOris);
            h(cond) = plot(modelOris, modelFit, 'Color', condColors(cond), 'LineWidth', 2);
            m = max([resp(:);modelFit(:)]);
            if m > maxi, maxi = m; end
        end
        hleg = legend(h, sprintf('%s (aR^2=%.2f)',conditionNames{1},adjRsqRaw(neuron,1)), ...
            sprintf('%s (aR^2=%.2f)',conditionNames{2},adjRsqRaw(neuron,2)));
        set(hleg, 'Box', 'off', 'Color', 'none')
        title('Tuning from raw responses')
        xlabel('Direction')
        ylabel('Calcium response')
        
        % Plot tuning curves from ridge regression filters
        subplot(1,3,2)
        hold on
        h = [0 0];
        for cond = 1:2
            plot(orientations(:,1), tuningCurvesRidge(:,neuron,cond), 'o', 'Color', condColors(cond))
            modelFit = oritune(squeeze(fitParsRidge(neuron,cond,:)), modelOris);
            h(cond) = plot(modelOris, modelFit, 'Color', condColors(cond), 'LineWidth', 2);
            m = max([resp(:);modelFit(:)]);
            if m > maxi, maxi = m; end
        end
        hleg = legend(h, sprintf('%s (aR^2=%.2f)',conditionNames{1},adjRsqRidge(neuron,1)), ...
            sprintf('%s (aR^2=%.2f)',conditionNames{2},adjRsqRidge(neuron,2)));
        set(hleg, 'Box', 'off', 'Color', 'none')
        title('Tuning from ridge regression filters')
        xlabel('Direction')
        ylabel('Calcium response')
        
        % Plot tuning curves from stimulus-time separated filters
        subplot(1,3,3)
        hold on
        h = [0 0];
        for cond = 1:2
            plot(orientations(:,1), tuningCurvesSep(:,neuron,cond), 'o', 'Color', condColors(cond))
            modelFit = oritune(squeeze(fitParsSep(neuron,cond,:)), modelOris);
            h(cond) = plot(modelOris, modelFit, 'Color', condColors(cond), 'LineWidth', 2);
            m = max([resp(:);modelFit(:)]);
            if m > maxi, maxi = m; end
        end
        hleg = legend(h, sprintf('%s (aR^2=%.2f)',conditionNames{1},adjRsqSep(neuron,1)), ...
            sprintf('%s (aR^2=%.2f)',conditionNames{2},adjRsqSep(neuron,2)));
        set(hleg, 'Box', 'off', 'Color', 'none')
        title('Tuning from stimulus-time separated filters')
        xlabel('Direction')
        ylabel('Calcium response')
        for pl = 1:3
            subplot(1,3,pl)
            xlim([-5 365])
            ylim([0 1.05*maxi])
        end
    end
end

%% Summary statistics
% Prediction errors
gratings.plotPredictionErrors(predictionErrors);
% Fitting parameters
%  Raw data
gratings.compareConditionedTuningParameters(fittingParameters.raw, baselineActivity, ...
   conditionNames, adjustedRsquared.raw);
%  Ridge regression
gratings.compareConditionedTuningParameters(fittingParameters.ridge, baselineActivity, ...
   conditionNames, adjustedRsquared.ridge);
%  Stimulus-time separated kernels
gratings.compareConditionedTuningParameters(fittingParameters.sep, baselineActivity, ...
   conditionNames, adjustedRsquared.sep);
% Goodness of fits
r2 = cat(2, adjustedRsquared.raw(:), adjustedRsquared.ridge(:), adjustedRsquared.sep(:));
gratings.compareGoodnessOfFits(r2, {'raw', 'ridge regression', 'separated'});
