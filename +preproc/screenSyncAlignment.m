function s = screenSyncAlignment(block, pdt, pdchan)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

pdt = pdt(:); % ensure it's a row vector

[~, vals] = kmeans(pdchan, 2); % find best 2 photodiode value clusters (i.e. it is triggered or untriggered)
thresh = mean(vals); % threshold is mean their mean
pdFlips = abs(diff(pdchan > thresh)) > 0; % look for the flips
% assume flip time is halfway between the threshold crossings
pdFlipTimes = mean([pdt([false ; pdFlips]) pdt([pdFlips ; false])], 2);

diffs = diff(pdFlipTimes) > .01;
% diffs = [true; or(diffs, [false; diffs(1:end-1)])];
diffs = or(diffs, [false; diffs(1:end-1)]);
pdFlipTimes =  pdFlipTimes(diffs); %sort([pdFlipTimes(1); pdFlipTimes([diffs; false]); ...
    %pdFlipTimes([false; diffs])]);

blockFlipsTimes = block.stimWindowUpdateTimes;

p = pdFlipTimes;
b = blockFlipsTimes;
ax = zeros(1,3);
f1 = figure('Position', [4 50 502 946]);
ax(1) = gca;
f2 = figure('Position', [516 558 1402 420]);
ax(2) = gca;
f3 = figure('Position', [516 55 1402 420]);
ax(3) = gca;
linkaxes(ax, 'x');
for k = 1:3
    hold(ax(k), 'on');
end
diodeHeight = 1;
blockHeight = 1.5;
while true
    for k = 1:3
        cla(ax(k))
    end
    figure(f2)
    plot(pdt, pdchan)
    h1 = plot(p, ones(1, length(p)).* diodeHeight, 'o');
    h2 = plot(b, ones(1, length(b)).* blockHeight, 'x');
    xlim(pdt([1 end]))
    legend([h1 h2], {'diode','block'})
    xlabel('Time (s)')
    title(sprintf('Diode: %d events; block: %d events', length(p), length(b)))
    
    L = min(length(p), length(b));
    co = robustfit(b(1:L), p(1:L));
    figure(f1)
    plot(b(1:L), p(1:L), '.')
    plot(b([1 end]), b([1 end]) * co(2) + co(1))
    xlabel('Diode events')
    ylabel('Block events')
    
    figure(f3)
    diffs = p(1:L) - (b(1:L) * co(2) + co(1));
    plot(p(1:L), diffs)
    ylim([min(diffs) max(diffs)])
    xlabel('Time (s)')
    ylabel('Residuals')
    
%     rect = getrect(ax(2));
    rect = drawrectangle(ax(2)).Position;
    if diodeHeight>rect(2) && diodeHeight<rect(2)+rect(4) % diode events
        ind = p>rect(1) & p<rect(1)+rect(3);
        p(ind) = [];
    elseif blockHeight>rect(2) && blockHeight<rect(2)+rect(4) %block events
        ind = b>rect(1) & b<rect(1)+rect(3);
        b(ind) = [];
    else % no events selected -> finished
        break
    end
end

pdFlipTimes = p;
blockFlipsTimes = b;

close(f1)
close(f2)
close(f3)

% if numel(pdFlipTimes) - 1 == numel(blockFlipsTimes)
%   % in most datasets the first photodiode flip occurs without a
%   % corresponding stimWindowUpdate element, so just ignore that one
%   pdFlipTimes = pdFlipTimes(2:end);
% end
% 
% % On b2 scope there is also an extra photodiode flip at the end of the experiment
% if numel(pdFlipTimes) - 2 == numel(blockFlipsTimes)
%     pdFlipTimes([1 end]) = [];
% end
% 
% % photodiode now picks up two additional flips at the end of sessions for
% % unknown reasons. 2016/06/22 SF
% if numel(pdFlipTimes) - 3 == numel(blockFlipsTimes)
%     pdFlipTimes([1 end-1 end]) = [];
% end

[co, stats] = robustfit(blockFlipsTimes, pdFlipTimes);
assert(stats.ols_s < 0.02, 'Significant error in fit')

s.coeff = co';
s.blockToPdTimeFrame = @(t)t*co(2) + co(1);
s.pdToBlockTImeFrame = @(t)(t - co(1))/co(2);
%offset to place every block flip after corresponding photodiode flip
lag = -max(blockFlipsTimes*co(2) - pdFlipTimes);
toPDTimeFrameLag = @(t)t*co(2) + lag;

if isfield(block.trial, 'stimulusCueStartedTime')
    s.stimOnTimes = follows(...
        toPDTimeFrameLag([block.trial.stimulusCueStartedTime]), pdFlipTimes);
end

if isfield(block.trial, 'stimulusCueEndedTime')
    s.stimOffTimes = follows(...
        toPDTimeFrameLag([block.trial.stimulusCueEndedTime]), pdFlipTimes);
end

    function t = follows(a, b)
        n = numel(a);
        t = zeros(size(a));
        ti = t;
        for ii = 1:n
            ti(ii) = find(b > a(ii), 1);
            t(ii) = b(ti(ii));
        end
        
        d = t - a;
        range = max(d) - min(d);
        assert((range/mean(d)) < 4.2, 'delta range is much larger than the mean');
    end

end

