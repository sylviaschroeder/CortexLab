% load q structure for dataset
% load pupil data: results + eyeLog
pupil = results;
pupil = nonVis.getPupilDiam(pupil);
pupilTime = [eyeLog.TriggerData.Time];
ind = isnan(pupil);
pupil = pupil(~ind);
pupilTime = pupilTime(~ind);
time = q.frameTimes;
pupil = interp1(pupilTime, pupil, time)';

F = q.cellMat;
F = bsxfun(@rdivide, F - mean(F), std(F));

perf = zeros(length(time), 1);
for tr = 1:size(q.stimTimes,2)
    n = find(hist(q.stimTimes(:,tr), q.frameTimes));
    perf(n(1):n(2)) = q.trialConditions(tr,5);
end
perf = (perf+1) ./ 2;
perf_smooth = smooth(perf, 200);

r1 = corr(pupil, F); % pupil vs. calcium
[~,ind_sort] = sort(r1, 'descend');

figure('Position', [1 41 1920 1083])
subplot(8,1,1)
plot(time, pupil)
xlim(time([1 end]))
ylabel('Pupil diameter')
axis tight
set(gca, 'box', 'off')

subplot(8,1,2)
plot(time, perf_smooth, 'k')
axis tight
xlim(time([1 end]))
ylim([0 1])
ylabel('Performance')
set(gca, 'box', 'off')

subplot(8,1,3:5)
hold on
for k = 1:10
    plot(time,F(:,ind_sort(k))-(k-1)*5,'r')
end
xlim(time([1 end]))
title('Most positively correlated with pupil')
ylabel('\DeltaF/F')
axis tight
set(gca, 'box', 'off')

subplot(8,1,6:8)
hold on
for k = 1:10
    plot(time,F(:,ind_sort(end-k+1))-(k-1)*5,'b')
end
xlim(time([1 end]))
xlabel('Time (s)')
title('Most negatively correlated with pupil')
ylabel('\DeltaF/F')
axis tight
set(gca, 'box', 'off')

r2 = corr(perf_smooth, F); % performance vs. calcium

binSize = .05;
mini = round(min([r1 r2]) / binSize) * binSize;
maxi = round(max([r1 r2]) / binSize) * binSize;
bins = mini : binSize : maxi;
figure('Position', [215 670 1620 420])
subplot(1,3,1)
hist(r1, bins)
xlim([mini-binSize maxi+binSize])
xlabel('Corr. coeff.')
title('Pupil and calcium')
subplot(1,3,2)
hist(r2, bins)
xlim([mini-binSize maxi+binSize])
xlabel('Corr. coeff.')
title('Performance and calcium')
subplot(1,3,3)
plot(pupil, perf_smooth, 'k.')
xlabel('Pupil')
ylabel('Performance')
title(sprintf('rho = %.3f', corr(pupil,perf_smooth)))