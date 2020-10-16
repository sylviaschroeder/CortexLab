function b = alignBlock2TL(block, tl, doPlot)

% upsample wheel signal in block to match sampling rate of timeline
dtTL = median(diff(tl.rawDAQTimestamps));
t_bl_up = (block.inputs.wheelTimes(1) : dtTL : block.inputs.wheelTimes(end))';
wh_bl_up = interp1(block.inputs.wheelTimes, block.inputs.wheelValues, ...
    t_bl_up);

% extract wheel signal from timeline
wh_tl_raw = wheel.correctCounterDiscont( ...
    tl.rawDAQData(:,strcmp('rotaryEncoder', {tl.hw.inputs.name})))';

% smooth first derivative of both wheel signals
% we can't just use the raw signals because they are on different scales
% (additive shift); derivates are too noisy so they need smoothing
win = normpdf(-20:20, 0, 4);
wh_bl = conv(diff(wh_bl_up), win, 'same');
wh_tl = conv(diff(wh_tl_raw), win, 'same');

% find delay between both wheel signals
d = finddelay(wh_bl, wh_tl);

% determine difference between times to contruct linear coefficients b
% (assume the times are only shifted, not scaled -> slope is set to 1)
b = [1; 0];
if d >= 0
    b(2) = tl.rawDAQTimestamps(1+d) - t_bl_up(1);
else
    b(2) = tl.rawDAQTimestamps(1) - t_bl_up(1-d);
end

if doPlot
    figure('Position', [80 555 1715 420])
    tiledlayout(1,2);
    nexttile
    hold on
    plot(wh_tl, 'k')
    plot((1:length(wh_bl)) + d, wh_bl, 'r')
    
    nexttile
    hold on
    plot(wh_tl_raw, 'k')
    plot((1:length(wh_bl_up)) + d, wh_bl_up, 'r')
end