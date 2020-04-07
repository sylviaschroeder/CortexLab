function plotPupilRunningAndCaTraces(traces, traceTimes, pupilData, ...
    pupilTimes, ballData, neuronIDs)

% normalize traces to min of zero and max of one
normTraces = traces - repmat(min(traces, [], 1), size(traces,1), 1);
normTraces = normTraces ./ repmat(max(normTraces, [], 1), size(normTraces,1), 1);
numFigs = ceil(size(traces,2) / 30);
numTracesPerFig = ceil(size(traces,2) / numFigs);
if isempty(pupilTimes)
    pupilTimes = traceTimes(end);
end
maxTime = min([traceTimes(end), pupilTimes(end), ballData.t(end)]);

if nargin < 6
    neuronIDs = 1:size(traces,2);
end

for fig = 1:numFigs
    figure('Position', [854 49 910 1068])
    ax = zeros(1,4);
    
    if ~isempty(pupilData)
        ax(1) = subplot(5,1,1);
        set(gca, 'Position', [0.1 0.94 0.85 0.05])
        xData = pupilData.x;
        xData(pupilData.blink | pupilData.goodFit == 0) = NaN;
        mx = nanmean(xData);
        sx = nanstd(xData);
        xData(xData > mx + 5*sx | xData < mx - 5*sx) = NaN;
        plot(pupilTimes, xData, 'k')
        hold on
        plot(pupilTimes(pupilData.goodFit == 0), nanmean(xData), 'b.')
        plot(pupilTimes(pupilData.blink), nanmean(xData), 'r.')
        set(gca, 'XTickLabel', [], 'box', 'off')
        xlim([0 maxTime])
        rng = max(xData) - min(xData);
        ylim([min(xData) - 0.05*rng, max(xData) + 0.05*rng])
        ylabel('Pupil: X')
        
        ax(2) = subplot(5,1,2);
        set(gca, 'Position', [0.1 0.87 0.85 0.05])
        yData = pupilData.y;
        yData(pupilData.blink | pupilData.goodFit == 0) = NaN;
        my = nanmean(yData);
        sy = nanstd(yData);
        yData(yData > my + 5*sy | yData < my - 5*sy) = NaN;
        plot(pupilTimes, yData, 'k')
        hold on
        plot(pupilTimes(pupilData.goodFit == 0), nanmean(yData), 'b.')
        plot(pupilTimes(pupilData.blink), nanmean(yData), 'r.')
        set(gca, 'XTickLabel', [], 'box', 'off')
        xlim([0 maxTime])
        rng = max(yData) - min(yData);
        ylim([min(yData) - 0.05*rng, max(yData) + 0.05*rng])
        ylabel('Pupil: Y')
        
        ax(3) = subplot(5,1,3);
        set(gca, 'Position', [0.1 0.8 0.85 0.05])
        area = pupilData.area;
        area(pupilData.blink | pupilData.goodFit == 0) = NaN;
        ma = nanmean(area);
        sa = nanstd(area);
        area(area > ma + 5*sa | area < ma - 5*sa) = NaN;
        diam = sqrt(4 * area / pi);
        plot(pupilTimes, diam, 'k')
        hold on
        plot(pupilTimes(pupilData.goodFit == 0), nanmean(diam), 'b.')
        plot(pupilTimes(pupilData.blink), nanmean(diam), 'r.')
        set(gca, 'XTickLabel', [], 'box', 'off')
        xlim([0 maxTime])
        ylim([0, 1.05*max(diam)])
        ylabel('Pupil: Diam.')
    end
    
    ax(4) = subplot(5,1,4);
    set(gca, 'Position', [0.1 0.73 0.85 0.05])
    plot(ballData.t, ballData.total, 'k')
    set(gca, 'XTickLabel', [], 'box', 'off')
    xlim([0 maxTime])
    ylim([0 1.05*max(ballData.total)])
    ylabel('Running')
    
    ax(5) = subplot(4,1,4);
    set(gca, 'Position', [0.1 0.06 0.85 0.65])
    ind = (fig-1)*numTracesPerFig+1 : min(fig*numTracesPerFig, size(traces,2));
    tr = normTraces(:, ind);
    plot(traceTimes, tr - repmat(1:size(tr,2), size(tr,1), 1), 'k')
    set(gca, 'YTick', (-size(tr,2):-1)+0.5, 'YTickLabel', flip(neuronIDs(ind)), ...
        'box', 'off')
    xlim([0 maxTime])
    ylim([-size(tr,2) 0])
    xlabel('Time (in s)')
    ylabel('Normalized Ca-traces (neuron IDs)')
    
    linkaxes(ax, 'x')
end