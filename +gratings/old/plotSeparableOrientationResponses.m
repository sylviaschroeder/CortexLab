function plotSeparableOrientationResponses(calciumTraces, traceTimes, ...
    filterPred, filterSepPred)

for iCell = 1:size(calciumTraces, 2)
    h = zeros(1, 3);
    figure('Position', [10 680 1900 420])
    h(1) = plot(traceTimes, calciumTraces(:,iCell), 'k');
    hold on
    h(2) = plot(traceTimes, filterPred(:,iCell), 'b');
    h(3) = plot(traceTimes, filterSepPred(:,iCell), 'r');
    
    set(gca, 'Position', [0.035 0.11 0.95 0.82])
    axis tight
    title(['ROI ' num2str(iCell)])
    xlabel('Time (in s)')
    ylabel('Calcium response')
    legend(h, 'Raw', 'Filters', 'Separated')
end