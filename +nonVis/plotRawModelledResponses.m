function plotRawModelledResponses(respRaw, respMod, baseline, time, ...
    stimDuration, titleStr)

% respRaw   [time x trial x stimulus]
% respMod   [time x trial x stimulus]
% baseline  [time x trial x stimulus]

numStim = size(respRaw, 3);
numTrial = size(respRaw, 2);
mini = min([respRaw(:); respMod(:)]);
maxi = max([respRaw(:); respMod(:)]);
yRange = maxi - mini;
yGap = 0.05*yRange;
xRange = diff(time([1 end]));
xGap = 0.1*xRange;
labels = {'raw', 'model', 'model baseline'};

figure('Position', [75 1 1846 1003])
hold on
h = NaN(1,3);
for stim = 1:numStim
%     yOffset = -(yRange+yGap) * (stim-1);
    xOffset = (xRange+xGap) * (stim-1);
    for trial = 1:numTrial
%         subplot(numStim, numTrial, (stim-1)*numStim+trial)
%         xOffset = (xRange+xGap) * (trial-1);
        yOffset = -(yRange+yGap) * (trial-1);
        plot([0 0]+xOffset, [mini maxi]+yOffset, 'k:')
        plot([1 1]*stimDuration+xOffset, [mini maxi]+yOffset, 'k:')
        plot(time([1 end])+xOffset, [0 0]+yOffset, 'k:')
        h(1)=plot(time+xOffset, respRaw(:,trial,stim)+yOffset, 'k', 'LineWidth', 2);
        if ~isempty(respMod)
            h(2)=plot(time+xOffset, respMod(:,trial,stim)+yOffset, 'r', 'LineWidth', 2);
        end
        if ~isempty(baseline)
            h(3)=plot(time+xOffset, baseline(:,trial,stim)+yOffset, ...
                'b', 'LineWidth', 2);
        end
    end
end
set(gca, 'XTick', round([time(1) 0 stimDuration time(end)].*10)./10, ...
    'YTick', round([mini maxi].*10)./10)
xlim([time(1) time(end)+xOffset])
ylim([mini+yOffset round(maxi*10)/10])
xlabel('Stimulus (and Time in s)')
ylabel('Trial')
ind = ~isnan(h);
legend(h(ind), labels(ind))
if nargin>5
    title(titleStr)
end