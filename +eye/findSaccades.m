function [saccadeTimes, velocity, velocityAligned, bins] = ...
    findSaccades(x, y, minDist, scaleThresh, doPlot)

if nargin < 4
    scaleThresh = 1;
end
if nargin < 5
    doPlot = 0;
end

diffX = diff(x);
diffY = diff(y);

velocity = sqrt(diffX.^2 + diffY.^2);
acceleration = diff(velocity);

vel_log = log(velocity);
vel_log(isinf(vel_log)) = [];
b_vel = floor(min(vel_log)*10)/10 : 0.1 : ceil(max(vel_log)*10)/10;
n = hist(vel_log, b_vel);
n = log(n);
ind = ~isinf(n);
f = fit(b_vel(ind)', n(ind)', 'gauss1');
coeffs = coeffvalues(f);
thresh_vel = exp(sum(coeffs(2:3))) * scaleThresh; % threshold = mean + STD of velocity distribution (in log-log histogram)

reverse1 = acceleration(1:end-1) .* acceleration(2:end);
reverse2 = acceleration(1:end-2) .* acceleration(3:end);

rev1_ln = log(reverse1(reverse1 < 0) .* -1);
bins1 = -.5:0.5:max(rev1_ln)+1;
n1 = hist(rev1_ln, bins1);
f1 = fit(bins1(2:end)', n1(2:end)', 'poly5');
b = 0:0.1:max(rev1_ln);
mins = islocalmin(f1(b));
ind = find(mins,1);
thresh_1 = -exp(b(ind));

rev2_ln = log(reverse2(reverse2 < 0) .* -1);
bins2 = -.5:0.5:max(rev2_ln)+1;
n2 = hist(rev2_ln, bins2);
f2 = fit(bins2(2:end)', n2(2:end)', 'poly5');
b = 0:0.1:max(rev2_ln);
mins = islocalmin(f2(b));
ind = find(mins,1);
thresh_2 = -exp(b(ind));
if isempty(thresh_2)
    thresh_2 = min(reverse2)/2;
end

sacc = false(length(x),1);
sacc(1:end-1) = sacc(1:end-1) | velocity > thresh_vel;
% sacc(2:end-2) = sacc(2:end-2) | reverse1<thresh_1;
% sacc(2:end-3) = sacc(2:end-3) | reverse2<thresh_2;
% sacc(2:end-4) = sacc(2:end-4) | reverse3<thresh_3;
sacc = conv(double(sacc), ones(1, minDist));
sacc(length(x)+1:end) = [];
saccadeTimes = find(diff(sacc>0)>0) + 1;
velPrev = velocity(saccadeTimes-1) > thresh_vel/2;
saccadeTimes(velPrev) = saccadeTimes(velPrev) - 1;

bins = -minDist : 3*minDist;
ind = saccadeTimes + bins;
nonVal = ind<1 | ind>length(velocity);
ind(nonVal) = 1;
velocityAligned = velocity(ind);
velocityAligned(nonVal) = NaN;

if doPlot > 0
    t = 1:length(x);
    figure('Position', [1 41 1920 1083])
    ax = zeros(1,6);
    subplot(6,6,1:5) % 1
    plot(x,'k')
    hold on
    plot(t(saccadeTimes), x(saccadeTimes), '*r')
    ylabel('eye pos in x')
    ax(1) = gca;
    subplot(6,6,7:11) % 2
    plot(y,'k')
    hold on
    plot(t(saccadeTimes), y(saccadeTimes), '*r')
    ylabel('eye pos in y')
    ax(2) = gca;
    subplot(6,6,13:17) % 3
    plot(2:length(x),velocity,'k')
    hold on
    plot([1 length(x)], [1 1].*thresh_vel, 'r')
    ylabel('velocity')
    ax(3) = gca;
    subplot(6,6,19:23) % 4
    plot(2:length(x)-1,acceleration,'k')
    ylabel('acceleration')
    ax(4) = gca;
    subplot(6,6,25:29) % 5
    plot(2:length(x)-2,reverse1,'k')
    ylabel('rev in t+1')
    hold on
%     plot([1 length(x)],[1 1].*thresh_1, 'r')
    ax(5) = gca;
    subplot(6,6,31:35) % 6
    plot(2:length(x)-3,reverse2,'k')
    ylabel('rev in t+2')
    hold on
%     plot([1 length(x)],[1 1].*thresh_2, 'r')
    ax(6) = gca;
    linkaxes(ax, 'x')
    xlim([1 length(x)])
    xlabel('Time (in samples)')
    
    subplot(6,6,6)
    plot(bins, nanmean(velocityAligned), 'k')
    xlim(bins([1 end]))
    xlabel('Time from saccade')
    ylabel('Velocity')
    subplot(6,6,18)
    plot(b_vel, n, 'k')
    hold on
    plot(f, 'r')
    plot(log(thresh_vel), f(log(thresh_vel)), 'ro')
    xlim(b_vel([1 end]))
    xlabel('velocity (e^x)')
    legend off
    subplot(6,6,30) % 5
    plot(bins1(2:end), n1(2:end), 'k')
    hold on
    plot(f1)
    xlim(bins1([2 end]))
    legend off
    subplot(6,6,36) % 6
    plot(bins2(2:end), n2(2:end), 'k')
    hold on
    plot(f2)
    xlim(bins2([2 end]))
    legend off
    xlabel('rev magn in e^x')
end