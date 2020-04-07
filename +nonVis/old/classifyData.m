function [negativeCondition, positiveCondition, conditionTime, labels, ...
    conditionSignal] = ...
    classifyData(meta, condition, doPlot, threshold)

% meta                  info structure describing imaging data set
% condition             String: 'running', 'pupilDiam', 'pupilDiamDeriv', or 'eyeMove'
% doPlot                if 1, classified data are plotted

% negativeCondition     [1 x time]; each entry is 0 or 1 specifying at each
%                       time point whether the opposite of the condition is
%                       true (e.g. case running: is 1 if mouse is 
%                       considered to be stationary, is 0 otherwise)
% positiveCondition     [1 x time]; each entry is 0 or 1 specifying at each
%                       time point whether the condition is true (e.g. case
%                       running: is 1 if mouse is considered to be running,
%                       0 otherwise)
% conditionTime         [1 x time]; time of condition data (e.g. sample
%                       times for running)
% labels                {<label of negative condition>, <label of positive
%                       condition>}
% conditionSignal       raw data used for classification (e.g. ball
%                       movements in case of running)

conditions = {'running', 'pupilDiam', 'pupilDiamDeriv', 'eyeMove'};
condLabels = {'Ball movement', 'Pupil diamter', 'Deriv. of pupil diameter', ...
    'Eye movement'};

switch condition
    case conditions{1}
        ballData = nonVis.getRunningSpeed(meta);
        conditionSignal = ballData.total;
        conditionTime = ballData.t;
        dataValues = unique(conditionSignal);
        nValues = hist(conditionSignal, dataValues);
        if nargin < 4
            thresh = find(nValues > ...
                max(1, 0.005 * length(conditionSignal)), 10); % should be around 1.0
            thresholdPositive = dataValues(thresh(5)); % running threshold
                                                   % 2: ~0.5, used before 24.04.2015
                                                   % 3: around 1.0
                                                   % 4: around 1.5, ...
            thresholdNegative = thresholdPositive; % non-running threshold
        end
        minIntervalPositive = 0.5; % in sec; min. time of threshold excess
                                   % to count as running
        minIntervalNegative = 0.5; % in sec; min. time of threshold excess
                                   % to count as stationary
        smoothingWindow = 0.3; % in sec
        labels = {'not running', 'running'};
        condName = condLabels{1};
    case conditions{2}
        [pupilData, pupilTimes] = nonVis.loadPupilData(meta);
        conditionSignal = sqrt(4 * pupilData.area / pi);
        conditionSignal(pupilData.blink | ~pupilData.goodFit) = NaN;
        condMean = nanmean(conditionSignal);
        condStd = nanstd(conditionSignal);
        conditionSignal(conditionSignal > condMean+5*condStd | ...
            conditionSignal < condMean-5*condStd) = NaN;
        conditionTime = pupilTimes;
        conditionTime(length(pupilData.area)+1 : end) = [];
        if nargin < 4
            thresholdPositive = prctile(conditionSignal, 50);
            thresholdNegative = thresholdPositive;
        end
        minIntervalPositive = 0.1;
        minIntervalNegative = 0.1;
        smoothingWindow = 2;
        labels = {'pupil constricted', 'pupil dilated'};
        condName = condLabels{2};
    case conditions{3}
        [pupilData, pupilTimes] = nonVis.loadPupilData(meta);
        conditionTime = pupilTimes;
        conditionTime(length(pupilData.area)+1 : end) = [];
        conditionSignal = sqrt(4 * pupilData.area / pi);
        ind = isnan(conditionSignal);
        pup = interp1(conditionTime(~ind), conditionSignal(~ind), pupilTimes(ind), 'pchip', 0);
        conditionSignal(ind) = pup;
        sr = 1 / median(diff(conditionTime));
        % [b, a] = butter(4, (2/sr)); % low-pass butterworth filter of order 4 and
        % cut-off frequency at 1 Hz (as in Reimer,...,
        % Tolias, 2014)
        [b, a] = butter(4, (2/sr) / 2); % cut-off now at 0.5 Hz
        conditionSignal = filtfilt(b, a, conditionSignal);
        conditionSignal = diff(conditionSignal);
        conditionTime = conditionTime((1:end-1)) + 1/sr/2;
        conditionSignal(pupilData.blink | ~pupilData.goodFit) = NaN;
        conditionSignal = interp1(conditionTime, conditionSignal, popTimes, 'pchip', NaN)';
        if nargin < 4
            thresholdPositive = 0;
            thresholdNegative = 0;
        end
        minIntervalPositive = 0.1;
        minIntervalNegative = 0.1;
        smoothingWindow = 0;
        labels = {'pupil constricting', 'pupil dilating'};
        condName = condLabels{3};
    case conditions{4}
        [pupilData, pupilTimes] = nonVis.loadPupilData(meta);
        conditionTime = pupilTimes;
        conditionTime(length(pupilData.x)+1 : end) = [];
        condSamplingRate = 1 / median(diff(conditionTime));
        x = pupilData.x;
        x(pupilData.blink | ~pupilData.goodFit) = NaN;
        x = smooth(x, 0.5 * condSamplingRate, 'lowess');
        [binCounts, bins] = hist(x, 50);
        [~, maxInd] = max(binCounts);
        mostCommon = bins(maxInd);
        x = x - mostCommon;
        y = pupilData.y;
        y(pupilData.blink | ~pupilData.goodFit) = NaN;
        y = smooth(y, 0.5 * condSamplingRate, 'lowess');
        [binCounts, bins] = hist(y, 50);
        [~, maxInd] = max(binCounts);
        mostCommon = bins(maxInd);
        y = y - mostCommon;
        conditionSignal = sqrt(x.^2 + y.^2);
        
%         conditionSignal = sqrt(pupilData.x.^2 + pupilData.y.^2);
%         conditionSignal(pupilData.blink | ~pupilData.goodFit) = NaN;
%         conditionSignal = smooth(conditionSignal, ...
%             0.5 * condSamplingRate, 'lowess');
%         [binCounts, bins] = hist(conditionSignal, 50);
%         [~, maxInd] = max(binCounts);
%         mostCommon = bins(maxInd);
%         conditionSignal = abs(conditionSignal - mostCommon);  

        if nargin < 4
            thresholdPositive = prctile(conditionSignal, 66);
            thresholdNegative = thresholdPositive;
        end
        minIntervalPositive = 0;
        minIntervalNegative = 1;
        smoothingWindow = 0;
        labels = {'eye stationary', 'eye moved away'};
        condName = condLabels{4};
end

if nargin > 3
    thresholdPositive = threshold;
    thresholdNegative = threshold;
end

condSamplingRate = 1 / median(diff(conditionTime));
if smoothingWindow > 0
    smoothedCondSignal = smooth(conditionSignal, ...
        smoothingWindow * condSamplingRate, 'lowess');
else
    smoothedCondSignal = conditionSignal;
end
positiveCondition = nonVis.getThresholdedPeriods(smoothedCondSignal, ...
    conditionTime, thresholdPositive, 'above', minIntervalPositive);
negativeCondition = nonVis.getThresholdedPeriods(smoothedCondSignal, ...
    conditionTime, thresholdNegative, 'below', minIntervalNegative);

if doPlot == 1
    plotData(conditionSignal, conditionTime, positiveCondition, ...
        negativeCondition, condName, labels)
end

end

function plotData(data, time, positiveCond, negativeCond, yName, labels)
colors = [0 0 0; 1 0 0];

figure('Position', [565 680 1320 420])
handles = [0 0 0 ];
handles(3) = plot(time, data, 'Color', [0.5 0.5 0.5]);
hold on
% highValue = prctile(data, 95);
% lowValue = prctile(data, 5);
% dataRange = highValue - lowValue;
neg = negativeCond;
neg(neg == 0) = NaN;
handles(1) = plot(time, data .* neg, 'Color', colors(1,:), 'LineWidth', 1);
% handles(1) = plot(time, neg * highValue, ...
%     'Color', colors(1,:), 'LineWidth', 2);
pos = positiveCond;
pos(pos == 0) = NaN;
handles(2) = plot(time, data .* pos, 'Color', colors(2,:), 'LineWidth', 1);
% handles(2) = plot(time, pos .* (highValue + 0.05*dataRange), ...
%     'Color', colors(2,:), 'LineWidth', 2);
xlabel('Time (in s)')
ylabel(yName)
legend(handles, [labels, {'unclassified'}], 'Location', 'NorthWest')
axis tight

end