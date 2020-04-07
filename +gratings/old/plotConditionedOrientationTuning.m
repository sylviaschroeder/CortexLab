function [fitParameters, adjustedRsquared, conditionNames] =  ...
    plotConditionedOrientationTuning(meta, condition, plotData, ...
    method, planes)
% meta          data information; has to contain field .ROI
% condition     if empty, condition resulting in best partitioning of
%               population response will be chosen, other supply one of: 
%               'running', 'pupilDiam', 'pupilDiamDeriv', 'eyeMove'
% plotConditionData     if 1, condition data and classification will be
%                       plotted
% method        either 'separated' or 'raw'
% (planes)      image planes that will be considered in determining
%               population response, if not supplied only plane of meta
%               structure will be used

% fitParameters     [neuron x condition x parameter]; contains parameters
%                   resulting from fit for each neuron and each condition;
%                   parameters: preferred orientation (Dp), ampl. at pref.
%                   dir. (Rp), ampl. at opp. dir. (Rn), offset of tuning 
%                   curve from zero (Ro), tuning width (sigma)
% adjustedRsquared  [neuron x condition+1]; contains adjusted R_squared of
%                   fit for each neuron and each condition (plus all data 
%                   not divided into conditions)

lambda = 0.5; % for ridge regression, if empty optimal value will be 
              % determined via cross-validation
rSquaredThreshold = 0.3; % threshold 

parameters = {'Pref. direction', 'Amplitude at preferred direction', ...
    'Amplitude at null direction', 'Offset', 'Tuning width'};

if nargin < 5 || isempty(planes)
    planes = meta.iPlane;
end

% find optimal threshold for non-visual data to classify population
% responses
[thresholds, r_squared, conditions] = nonVis.findBestThreshold(meta, planes);
% if condition not supplied, choose the best condition of non-visual data
if isempty(condition)
    [~, ind] = min(r_squared);
    threshold = thresholds(ind);
    condition = conditions{ind};
else
    ind = strcmp(condition, conditions);
    threshold = thresholds(ind);
end

% classify data
[negativeCond, positiveCond, conditionTime, conditionNames] = ...
    nonVis.classifyData(meta, condition, plotData, threshold);

%% Do for each plane!
fitParameters = [];
adjustedRsquared = [];
for iPlane = planes
    % load meta
    filePath = fullfile(info.folderProcessed, ...
        sprintf('%s_plane%03d_ROI', info.basename2p, iPlane));
    load(filePath, 'meta')
    % get stimulus information
    [~, stimSequence, stimMatrix, frameTimes] = ...
        ssLocal.getStimulusResponseInfo(meta);
    orientations = gratings.getOrientations(stimSequence);
    
    switch method
        case 'separated'
            [~, tunings, offsets, ~, ~, errors] = ...
                gratings.getConditionedOrientationTuning_separated(meta.F, frameTimes, ...
                [negativeCond, positiveCond], conditionTime, stimMatrix, ...
                stimSequence, lambda);
            tuningCurves = bsxfun(@plus, tunings, offsets);
            if plotData == 1
                figure
                plot(errors.cond_regr, errors.nonCond_regr, 'k.')
                maxi = max([errors.cond_regr, errors.nonCond_regr]);
                axis([0 maxi 0 maxi])
                hold on
                plot([0 maxi], [0 maxi], 'k:')
                title('Prediction errors of regression filters')
                xlabel('Trials divided into conditions')
                ylabel('Trials not divided into conditions')
                
                figure
                plot(errors.cond_sep, errors.nonCond_sep, 'k.')
                maxi = max([errors.cond_regr, errors.nonCond_regr]);
                axis([0 maxi 0 maxi])
                hold on
                plot([0 maxi], [0 maxi], 'k:')
                title('Prediction errors of stimulus-time separated filters')
                xlabel('Trials divided into conditions')
                ylabel('Trials not divided into conditions')
            end
        case 'raw'
            tuningCurves = gratings.getConditionedOrientationTuning_raw(meta.F, ...
                frameTimes, [negativeCond, positiveCond], conditionTime, stimMatrix, ...
                stimSequence);
    end
    % tuningCurves: [stimulus x neuron x condition]
    
    [fitPars, adjRsqu] = gratings.fitTuningCurves(tuningCurves, ...
        orientations(:,1));
    % fitPars: [neuron x condition x parameter]
    
    fitParameters = [fitParameters; fitPars];
    adjustedRsquared = [adjustedRsquared; adjRsqu];
end

for par = 2:4
    figure('Position', [860 680 1040 420])
    
    subplot(1,2,1)
    [f, gof] = fit(fitParameters(:,1,par), fitParameters(:,2,par), 'a*x');
    plot(fitParameters(:,1,par), fitParameters(:,2,par), 'k.')
    hold on
    mini = min(reshape(fitParameters(:,[1 2],par),[],1));
    maxi = max(reshape(fitParameters(:,[1 2],par),[],1));
    plot([mini maxi], f([mini maxi]), 'r')
    axis([mini maxi mini maxi])
    a = coeffvalues(f);
    confidenceInt = confint(f);
    text(mini+0.05*(maxi-mini), mini+0.9*(maxi-mini), ...
        sprintf('y = a * x\na = %.2f (%.2f %.2f)\nadj R^2 = %.3f', ...
        a, confidenceInt(1), confidenceInt(2), gof.adjrsquare))
    title(sprintf('%s (all data used)', parameters{par}))
    xlabel(['Not ' condition])
    ylabel(condition)
    
    subplot(1,2,2)
    ind = all(adjustedRsquared(:,1:2) >= rSquaredThreshold, 2);
    [f, gof] = fit(fitParameters(ind,1,par), fitParameters(ind,2,par), 'a*x');
    plot(fitParameters(ind,1,par), fitParameters(ind,2,par), 'k.')
    hold on
    mini = min(reshape(fitParameters(:,[1 2],par),[],1));
    maxi = 1.05 * max(reshape(fitParameters(:,[1 2],par),[],1));
    plot([mini maxi], f([mini maxi]), 'r')
    axis([mini maxi mini maxi])
    a = coeffvalues(f);
    confidenceInt = confint(f);
    text(mini+0.05*(maxi-mini), mini+0.9*(maxi-mini), ...
        sprintf('y = a * x\na = %.2f (%.2f %.2f)\nadj R^2 = %.3f', ...
        a, confidenceInt(1), confidenceInt(2), gof.adjrsquare))
    title(sprintf('%s (only good fits used)', parameters{par}))
    xlabel(['Not ' condition])
    ylabel(condition)
end

% compare various combinations of parameters
newParameters = cat(3, sum(fitParameters(:,:,[2 4]), 3), ... % Rp + R0
    fitParameters(:,:,2) ./ fitParameters(:,:,3), ...        % Rp / Rn
    (fitParameters(:,:,2)-fitParameters(:,:,3)) ./ ...
    sum(fitParameters(:,:,[2 3]),3), ...                     % (Rp-Rn) / (Rp+Rn)
    fitParameters(:,:,2) ./ fitParameters(:,:,4));           % Rp / R0
parNames = {'Max. responses (Rp+R0)', ...
    'Direction selectivity (Rp/Rn)', ...
    'Direction selectivity (Rp-Rn)/(Rp+Rn)', ...
    'Tuning strength (Rp/R0)'};
for par = 1:size(newParameters,3)
    figure('Position', [860 680 1040 420])
    
    subplot(1,2,1)
    [f, gof] = fit(newParameters(:,1,par), newParameters(:,2,par), 'a*x');
    plot(newParameters(:,1,par), newParameters(:,2,par), 'k.')
    hold on
    mini = min(reshape(newParameters(:,[1 2],par),[],1));
    maxi = max(reshape(newParameters(:,[1 2],par),[],1));
    plot([mini maxi], f([mini maxi]), 'r')
    axis([mini maxi mini maxi])
    a = coeffvalues(f);
    confidenceInt = confint(f);
    text(mini+0.05*(maxi-mini), mini+0.9*(maxi-mini), ...
        sprintf('y = a * x\na = %.2f (%.2f %.2f)\nadj R^2 = %.3f', ...
        a, confidenceInt(1), confidenceInt(2), gof.adjrsquare))
    title(sprintf('%s (all data used)', parNames{par}))
    xlabel(['Not ' condition])
    ylabel(condition)
    
    subplot(1,2,2)
    ind = all(adjustedRsquared(:,1:2) >= rSquaredThreshold, 2);
    [f, gof] = fit(newParameters(ind,1,par), newParameters(ind,2,par), 'a*x');
    plot(newParameters(ind,1,par), newParameters(ind,2,par), 'k.')
    hold on
    mini = min(reshape(newParameters(:,[1 2],par),[],1));
    maxi = max(reshape(newParameters(:,[1 2],par),[],1));
    plot([mini maxi], f([mini maxi]), 'r')
    axis([mini maxi mini maxi])
    a = coeffvalues(f);
    confidenceInt = confint(f);
    text(mini+0.05*(maxi-mini), mini+0.9*(maxi-mini), ...
        sprintf('y = a * x\na = %.2f (%.2f %.2f)\nadj R^2 = %.3f', ...
        a, confidenceInt(1), confidenceInt(2), gof.adjrsquare))
    title(sprintf('%s (only good fits used)', parNames{par}))
    xlabel(['Not ' condition])
    ylabel(condition)
end