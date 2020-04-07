function cleanTraces = regressRedFromGreenTraces(tracesGreen, tracesRed, ...
    frameTimes, doPlot, folder)
% Performs a linear fit of tracesRed (traces from red channel) to low
% values (estimate of baseline) of tracesGreen (traces from green channel)
% and subtracts the fit. 
% Assumption: tracesGreen = cleanTraces + r * tracesRed

% tracesGreen   [time x ROIs]
% tracesRed     [time x ROIs]
% frameTimes    [1 x time]
% doPlot        Plots are generated if set to 1
% folder        str; where to save figures, nothing is saved if kept empty
%
% cleanTraces   [time x ROIs]

window = 400;
percentile = 5;

if nargin<4
    doPlot = 0;
end
if nargin<5
    folder = [];
end
if ~isempty(folder) && ~isdir(folder)
    mkdir(folder)
end
    

[highPassFilteredGreen, smoothedGreen] = ssLocal.removeSlowDrift( ...
    tracesGreen, frameTimes, window, percentile);
[highPassFilteredRed, smoothedRed] = ssLocal.removeSlowDrift( ...
    tracesRed, frameTimes, window, percentile);

% assume: greenTrace = realSignal - r * redTrace
% now estimate r and return realSignal
options.numN = 20; % number of values to base fit on 
options.minNp = 5; % minimum percentile of red trace to be considered
options.maxNP = 90; % maximum percentile of red trace to be considered
options.pCell = 10; % values in green trace up to this percentile are 
                    % considered for fit
options.noNeg = 0; % allow negative correlation between green and red traces
indNaN = any(isnan([highPassFilteredGreen; highPassFilteredRed]), 1);
hpfGreen = highPassFilteredGreen;
hpfGreen(:, indNaN) = [];
hpfRed = highPassFilteredRed;
hpfRed(:, indNaN) = [];
[~, fittingResults] = estimateNeuropil(hpfGreen', ...
    hpfRed', options);

ct = tracesGreen(:,~indNaN) - bsxfun(@times, fittingResults.corrFactor(:,2)', ...
    hpfRed);
cleanTraces = NaN(size(ct,1), size(indNaN,2));
cleanTraces(:,~indNaN) = ct;

corrFactors = NaN(size(tracesGreen,2),2);
corrFactors(~indNaN,:) = fittingResults.corrFactor;
if doPlot == 1
    for iCell = 1:size(tracesGreen,2)
        if any(isnan(tracesGreen(:,iCell)))
            continue
        end
        
        figure('Position',[262 45 1653 932])
        ax1 = subplot(2,4,1:3);
        hold on
        plot(frameTimes, tracesGreen(:,iCell), 'k');
        plot(frameTimes, smoothedGreen(:,iCell), 'b');
        plot(frameTimes, highPassFilteredGreen(:,iCell), 'r');
        plot(frameTimes, cleanTraces(:,iCell), 'c');
        title(sprintf('Trace %d - Green trace', iCell))
        
        ax2 = subplot(2,4,5:7);
        hold on
        plot(frameTimes, tracesRed(:,iCell), 'k');
        plot(frameTimes, smoothedRed(:,iCell), 'b');
        plot(frameTimes, highPassFilteredRed(:,iCell), 'r')
        title('Red trace')
        linkaxes([ax1 ax2], 'x')
        xlim([0 frameTimes(end)])
        
        subplot(2,4,8)
        hold on
        plot(highPassFilteredRed(:,iCell), highPassFilteredGreen(:,iCell), ...
            'k.', 'MarkerSize', 2)
        xmin = min(highPassFilteredRed(:,iCell));
        xmax = max(highPassFilteredRed(:,iCell));
        ymin = min(highPassFilteredGreen(:,iCell));
        ymax = max(highPassFilteredGreen(:,iCell));
        plot([xmin xmax], corrFactors(iCell,2).*[xmin xmax]+ ...
            corrFactors(iCell,1),'r')
        axis([xmin xmax ymin ymax])
        xlabel('Red trace')
        ylabel('Green trace')
        
        if ~isempty(folder)
            savefig(gcf, fullfile(folder, sprintf('Trace%04d',iCell)), 'compact')
            saveas(gcf, fullfile(folder, sprintf('Trace%04d.jpg',iCell)))
            close(gcf)
        end
    end
end