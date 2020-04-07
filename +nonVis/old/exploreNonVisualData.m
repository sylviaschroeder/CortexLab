function [popResp, popTimes, running, pupilDiam, pupilDiamDeriv, eyeMove] = ...
    exploreNonVisualData(infoROI, titleLine, doPlot, removeDrift, window)

popResp = [];
popTimes = [];
running =[];
pupilDiam = [];
pupilDiamDeriv = [];
eyeMove = [];

if nargin < 4
    removeDrift = 0;
end
if nargin < 5
    window = 400; % in sec to remove slow drift
end

[traces, popTimes] = ssLocal.matchSampleTimesAcrossPlanes(infoROI);

for iPlane = 1:length(infoROI)
    
    % get population response
    neuronInds = strcmp(infoROI(iPlane).ROI.CellClasses, 's');
    pResp = traces{iPlane}(:,neuronInds);
    if removeDrift == 1
        pResp = ssLocal.removeSlowDrift(pResp, popTimes, window);
    end
%     if isfield(infoROI(iPlane), 'movingFrames')
%         pResp(infoROI(iPlane).movingFrames,:) = NaN;
%     end
    pResp = bsxfun(@rdivide, bsxfun(@minus, pResp, nanmean(pResp)), nanstd(pResp));
    if ~isempty(popResp) && size(pResp,1) ~= size(popResp,1)
        mini = min(size(pResp,1), size(popResp,1));
        pResp(mini+1:end,:) = [];
        popResp(mini+1:end,:) = [];
        popTimes(mini+1:end) = [];
    end
    popResp = [popResp, pResp];
end
popResp = mean(popResp, 2);

% get running data
ballData = nonVis.getRunningSpeed(infoROI(iPlane));
if isempty(ballData)
    running = NaN(length(popTimes), 1);
else
    running = ballData.total;
    times = ballData.t;
    running = interp1(times, running, popTimes, 'pchip', NaN)';
end

% get pupil data
[pupilData, times] = nonVis.loadPupilData(infoROI(iPlane));
if ~isempty(pupilData)
    times(length(pupilData.x)+1:end) = [];
    % pupil diameter
    pupilDiam = sqrt(4 * pupilData.area / pi);
    pupilDiam(pupilData.blink | ~pupilData.goodFit) = NaN;
    pupilDiam = interp1(times, pupilDiam, popTimes, 'pchip', NaN)';
    % derivative of pupil diameter
    pupilDiamDeriv = sqrt(4 * pupilData.area / pi);
    ind = isnan(pupilDiamDeriv);
    pup = interp1(times(~ind), pupilDiamDeriv(~ind), times(ind), 'pchip', 0);
    pupilDiamDeriv(ind) = pup;
    sr = 1 / median(diff(times));
    % [b, a] = butter(4, (2/sr)); % low-pass butterworth filter of order 4 and
    % cut-off frequency at 1 Hz (as in Reimer,...,
    % Tolias, 2014)
    [b, a] = butter(4, (2/sr) / 2);
    pupilDiamDeriv = filtfilt(b, a, pupilDiamDeriv);
    pupilDiamDeriv = diff(pupilDiamDeriv);
    pupilDiamDeriv(pupilData.blink(1:end-1) | ~pupilData.goodFit(1:end-1)) = NaN;
    pupilDiamDeriv = interp1(times(1:end-1), pupilDiamDeriv, popTimes, 'pchip', NaN)';
    % eye movement
    sr = median(diff(times));
    x = pupilData.x;
    x(pupilData.blink | ~pupilData.goodFit) = NaN;
    x = smooth(x, 0.5 * sr, 'lowess');
    [binCounts, bins] = hist(x, 50);
    [~, maxInd] = max(binCounts);
    mostCommon = bins(maxInd);
    x = x - mostCommon;
    y = pupilData.y;
    y(pupilData.blink | ~pupilData.goodFit) = NaN;
    y = smooth(y, 0.5 * sr, 'lowess');
    [binCounts, bins] = hist(y, 50);
    [~, maxInd] = max(binCounts);
    mostCommon = bins(maxInd);
    y = y - mostCommon;
    eyeMove = sqrt(x.^2 + y.^2);
    eyeMove = interp1(times, eyeMove, popTimes, 'pchip', NaN)';
end

if doPlot == 1
    figure('Position', [730 60 1130 1040])
    labels = {'Population response', 'Running'};
    matrix = [popResp, running];
    if ~isempty(pupilData)
        labels = [labels, {'Pupil diameter', ...
            'Deriv. of pupil diam.', 'Eye movement'}];
        matrix = [matrix, pupilDiam, pupilDiamDeriv, eyeMove];
    end
    
    [~,ax,bigAx] = gplotmatrix(matrix, [], [], [], [], [], [], 'variable', ...
        labels, labels);
    
    indNaN = isnan(popResp);
    matrix(indNaN,:) = [];
    [rho, p] = corr(matrix);
    for row = 1:size(ax,2)-1
        for col = row+1:size(ax,2)
            axes(ax(row, col));
            xlim = get(ax(row, col), 'xlim');
            ylim = get(ax(row, col), 'ylim');
            text(0.05*range(xlim)+xlim(1), 0.9*range(ylim)+ylim(1), ...
                sprintf('rho=%.2f\np=%.3f', rho(row, col), p(row, col)))
        end
    end
    if nargin > 2 && ~isempty(titleLine)
        title(bigAx, titleLine, 'FontWeight', 'bold', 'Interpreter', 'none')
    end
end