folderData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\electrophys\';
folderEye = '\\zserver.cortexlab.net\Data\EyeCamera';
% folderEye = 'J:\Eye';
folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\Electrophys_SC\nonVisualEffects\';

%% Define datasets
% db_ephys_driftingGratings
db_ephys_darkScreens

%% Parameters
stabilityThreshold = 1; % for drift in spike rate for correlation analysis
windowStd = 0.5; % in sec; STD of Gaussian window to convolve spike rates
% nonVisualSignal = 'pupil'; %'running' or 'pupil';
% smoothing (low-pass filter) of non-visual signal before
smoothStd = 0.25; %in sec
ignoreLaserTrials = 1;

%% Correlation with behaviour (continuous traces)
draws = 200;
corrs = struct([]);
for k = 1:length(db)
    fprintf('Dataset %d of %d\n', k, length(db))
    corrs(k).subject = db(k).subject;
    corrs(k).date = db(k).date;
    
    if isfield(db, 'exp')
        corrs(k).exp = db(k).exp;
        % load data (time, spikeCounts, cellIDs, depths, stimMatrix,
        % directions, blanks, stimTimes, laserOn, nSpikes, sampleRate, 
        % spikeWidths, timelineToEphys, waveforms)
        load(fullfile(folderData, db(k).subject, db(k).date, ...
            sprintf('%02d_data.mat', db(k).exp)));
        runningFile = fullfile(folderData, db(k).subject, db(k).date, ...
            sprintf('%02d_running.mat', db(k).exp));
        pupilFile = fullfile(folderEye, db(k).subject, db(k).date, ...
            sprintf('%02d_eye_processed.mat', db(k).exp));
        pupilTimeFile = fullfile(folderEye, db(k).subject, db(k).date, ...
            sprintf('%02d_eyeTime.mat', db(k).exp));
    else
        corrs(k).exp = 'dark';
        % load data (cellIDs, depths, nSpikes, sampleRate, spikeCounts, 
        % spikeWidths, time, timelineToEphys, waveforms)
        load(fullfile(folderData, db(k).subject, db(k).date, 'dark_data.mat'));
        runningFile = fullfile(folderData, db(k).subject, db(k).date, ...
            'dark_running.mat');
        pupilFile = fullfile(folderEye, db(k).subject, db(k).date, ...
            'eye_processed.mat');
        pupilTimeFile = fullfile(folderEye, db(k).subject, db(k).date, ...
            'eyeTime.mat');
    end
    corrs(k).cellIDs = cellIDs';
    
    % (1) running
    if exist(runningFile, 'file') > 0
        % load running
        load(runningFile);
    else
        % load Timeline
        tlFile = dat.expFilePath(db(k).subject, db(k).date, db(k).tlExp, ...
            'timeline', 'master');
        load(tlFile)
        wheel = double(Timeline.rawDAQData(:, ...
            strcmp({Timeline.hw.inputs.name}, 'rotaryEncoder')));
        tlTime = Timeline.rawDAQTimestamps .* timelineToEphys(1) + ...
            timelineToEphys(2);
        wheelData = nonVis.getRunningSpeed_wheel(wheel, tlTime, smoothStd); % units per sec
        running = interp1(wheelData.t, wheelData.total, time, 'pchip');
        if all(running == 0)
            running = NaN;
        end
        save(runningFile, 'running');
    end
    % (2) pupil
    if exist(pupilFile, 'file') > 0
        data = load(pupilFile, 'results');
        pupil = data.results;
        pupil = nonVis.getPupilDiam(pupil, 16);
        if exist(pupilTimeFile, 'file') > 0
            data = load(pupilTimeFile);
            pupilTime = data.eyeTime;
        else
            data = load(fullfile(folderEye, db(k).subject, db(k).date, ...
                num2str(db(k).tlExp), 'eye_timeStamps.mat'));
            pupilTime = data.tVid .* timelineToEphys(1) + timelineToEphys(2);
        end
        nanTimes = pupilTime(isnan(pupil));
        nanTimes = hist(nanTimes, time)>0;
        pupil = interp1(pupilTime, pupil, time, 'pchip');
        pupil(nanTimes) = NaN;
    else
        pupil = [];
    end
    
    % only consider "stable" neurons (no big drift in terms of spike rate)
    % collapse spikes into 1 s bins, divide time series into 10 epochs of
    % equal size, use mean and std of all epochs to evaluate stability
    sr = 1/diff(time([1 2]));
    toAdd = round(ceil(length(time)/sr)*sr - length(time));
    sc = [spikeCounts; NaN(toAdd, size(spikeCounts,2))];
    sc = reshape(sc, round(sr), ceil(length(time)/sr), size(spikeCounts,2));
    sc = squeeze(nansum(sc,1));
    toAdd = ceil(size(sc,1)/10)*10 - size(sc,1);
    sc = [sc; NaN(toAdd, size(sc,2))];
    sc = reshape(sc, size(sc,1)/10, 10, size(sc,2));
    means = squeeze(nanmean(sc, 1));
    stds = squeeze(nanstd(sc, 0, 1));
    ratios = mean(stds,1) ./ std(means, 0, 1);
    corrs(k).isStable = ratios' > stabilityThreshold;
    indStable = find(ratios > stabilityThreshold);
    spikeCounts = spikeCounts(:, indStable);
    
    winSize = round(windowStd * sr);
    win = normpdf(-5*winSize:5*winSize, 0, winSize);
    spikeRates = zeros(size(spikeCounts));
    for iCell = 1:size(spikeRates,2)
        spikeRates(:,iCell) = conv(spikeCounts(:,iCell), win, 'same');
    end
    if exist('stimMatrix', 'var') && ignoreLaserTrials == 1
        ind = sum(stimMatrix(laserOn,:),1) > 0;
        spikeRates(ind,:) = NaN;
    end
    
    ind = isnan(running) | any(isnan(spikeRates),2);
    if ~isempty(pupil)
        ind = ind | isnan(pupil);
    end
    spikeRates(ind,:) = [];
    running(ind) = [];
    [r,p] = corr(running, spikeRates);
    corrs(k).running.rhos = NaN(length(cellIDs),1);
    corrs(k).running.rhos(indStable) = r;
    corrs(k).running.pValues = NaN(length(cellIDs),1);
    corrs(k).running.pValues(indStable) = p;
    
    corrs(k).pupil.rhos = NaN(length(cellIDs),1);
    corrs(k).pupil.pValues = NaN(length(cellIDs),1);
    if ~isempty(pupil)
        pupil(ind) = [];
        [r,p] = corr(pupil, spikeRates);
        corrs(k).pupil.rhos(indStable) = r;
        corrs(k).pupil.pValues(indStable) = p;
    end
    corrs(k).cellDepths = db(k).yPosRange(2) - depths';
    
    % generate null distribution by shifting non-visual signal in time
    shifts = randi(length(running), draws, 1);
    runningShifted = NaN(length(running), draws);
    for j = 1:draws
        runningShifted(:,j) = circshift(running, shifts(j));
    end
    r = corr(runningShifted, spikeRates);
    corrs(k).running.nullRhos = NaN(length(cellIDs), draws);
    corrs(k).running.nullRhos(indStable,:) = r';
    corrs(k).pupil.nullRhos = NaN(length(cellIDs), draws);
    if ~isempty(pupil)
        pupilShifted = NaN(length(pupil), draws);
        for j = 1:draws
            pupilShifted(:,j) = circshift(pupil, shifts(j));
        end
        r = corr(pupilShifted, spikeRates);
        corrs(k).pupil.nullRhos(indStable,:) = r';
    end
end
pars.stabilityThreshold = stabilityThreshold;
pars.windowStd = windowStd;
pars.smoothStd = smoothStd;
pars.ignoreLaserTrials = ignoreLaserTrials;
if isfield(db, 'exp')
    save(fullfile(folderResults, ...
        'corrsDuringGratings_pupilAndRunning.mat'), 'corrs', 'pars')
else
    save(fullfile(folderResults, ...
        'corrsDuringDarkScreens_pupilAndRunning.mat'), 'corrs', 'pars')
end

%% Plot population results
if isfield(db, 'exp')
    data = load(fullfile(folderResults, ...
        'corrsDuringGratings_pupilAndRunning.mat'));
    stimulus = 'Gratings';
else
    data = load(fullfile(folderResults, ...
        'corrsDuringDarkScreens_pupilAndRunning.mat'));
    stimulus = 'Black screens';
end
corrs = data.corrs;
fields = {'running', 'pupil'};
rhos = cell(1,length(fields));
isSgnf = cell(1,length(fields));
null = cell(1,length(fields));
depths = [];
for k = 1:length(corrs)
    for f = 1:length(fields)
        r = corrs(k).(fields{f}).rhos;
        confIntv = prctile(corrs(k).(fields{f}).nullRhos, [2.5 97.5], 2);
        rhos{f} = [rhos{f}; r];
        isSgnf{f} = [isSgnf{f}; r>confIntv(:,2) | r<confIntv(:,1)];
        null{f} = [null{f}; corrs(k).(fields{f}).nullRhos];
    end
    depths = [depths; corrs(k).cellDepths];
end

labels = {'running speed', 'pupil diameter'};
% Histogram of correlation coefficients
for f = 1:length(fields)
    figure
%     bins = -.5:.05:.8;
    bins = -1:.05:1;
    n1 = hist(rhos{f}(isSgnf{f} == 1),bins)';
    n2 = hist(rhos{f}(isSgnf{f} == 0),bins)';
    bar(bins, [n1,n2], 'stacked')
    colormap([0 0 0; 1 1 1])
    xlabel(['Correlation: firing rate and ' labels{f}])
    ylabel('#Neurons')
    title([stimulus ': corr. with ' labels{f}])
    xlim([bins(1)-.05 bins(end)+.05])
    set(gca, 'box', 'off')
    legend('p < 0.05', 'p >= 0.05')
    legend('boxoff')
    
    % Depth plot
    figure
    hold on
    scatter(rhos{f}(isSgnf{f}==1), depths(isSgnf{f}==1), 'ko', ...
        'MarkerFaceColor', 'k')
    scatter(rhos{f}(isSgnf{f}==0), depths(isSgnf{f}==0), 'ko')
%     scatter(rhos{f}(rhos{f}>=0), depths(rhos{f}>=0), 'o', 'MarkerEdgeColor', 'none', ...
%         'MarkerFaceColor', 'r', 'MarkerFaceAlpha', .5);
%     scatter(-rhos{f}(rhos{f}<0), depths(rhos{f}<0), 'o', 'MarkerEdgeColor', 'none', ...
%         'MarkerFaceColor', 'b', 'MarkerFaceAlpha', .5)
    set(gca, 'YDir','reverse')
    xlim([bins(1)-.05 bins(end)+.05])
    xlabel(['Correlation: firing rate and ' labels{f}])
    ylabel('Depth from surface of SC')
    title([stimulus ': corr. with ' labels{f}])
    legend('pos corr','neg corr')
    legend('boxoff')
    
    % Cumulative distributions: data vs. null
    x = sort(rhos{f}(~isnan(rhos{f})), 'ascend');
    y = (1:length(x)) ./ length(x);
    xNull = sort(null{f}(~isnan(null{f})), 'ascend');
    yNull = (1:length(xNull)) ./ (length(xNull));
    xNull = [x(1); xNull; x(end)];
    yNull = [0 yNull 1];
    figure
    hold on
    plot(x, y, 'k')
    plot(xNull, yNull, ':', 'Color', [.5 .5 .5])
    xlim(x([1 end]))
    xlabel(['Correlation with ' labels{f}])
    ylabel('Proportion of neurons')
    title([stimulus ': corr. with ' labels{f}])
    legend('data', 'null distribution', 'Location', 'NorthWest')
end

%% Plot traces of most and least correlated neurons
% num neurons for each condition
n = 10;
% cols = [.7 0 0; .3 .3 .3; 0 0 .7];
cols = [.7 0 0; 0 0 .7];
minSamples = 200;
labels = {'running speed', 'pupil diameter'};
reference = 'running'; % 'running' or 'pupil'
f = 1; % 1 or 2

% load correlation results (corrs)
if isfield(db, 'exp')
    load(fullfile(folderResults, 'corrsDuringGratings_pupilAndRunning.mat'))
else
    load(fullfile(folderResults, 'corrsDuringDarkScreens_pupilAndRunning.mat'))
end
depths = [];
for k = 1:length(corrs)
    depths = [depths; corrs(k).cellDepths];
end

for k = 1:length(db)
    fprintf('Dataset %d of %d\n', k, length(db))
    if isfield(db, 'exp')
        prefix = sprintf('%02d', db(k).exp);
    else
        prefix = 'dark';
    end
    % load data (time, spikeCounts, cellIDs, dephts, stimMatrix,
    % directions, blanks, stimTimes, laserOn)
    load(fullfile(folderData, db(k).subject, db(k).date, ...
        [prefix '_data.mat']));
    % load running (in units per sec)
    load(fullfile(folderData, db(k).subject, db(k).date, ...
        [prefix '_running.mat']));
    cmPerUnit = 2*pi * 8.75 / (4 * 1024);
    running = running * cmPerUnit; % now in cm/s
    % load pupil
    if isfield(db, 'exp')
        pupilFile = sprintf('%02d_eye_processed.mat', db(k).exp);
        pupilTimeFile = sprintf('%02d_eyeTime.mat', db(k).exp);
    else
        pupilFile = 'eye_processed.mat';
        pupilTimeFile = 'eyeTime.mat';
    end
    if exist(fullfile(folderEye, db(k).subject, db(k).date, ...
            pupilFile), 'file')
        data = load(fullfile(folderEye, db(k).subject, db(k).date, ...
            pupilFile), 'results');
        pupil = data.results;
        pupil = nonVis.getPupilDiam(pupil, 16);
        if exist(fullfile(folderEye, db(k).subject, db(k).date, ...
                pupilTimeFile), 'file') > 0
            data = load(fullfile(folderEye, db(k).subject, db(k).date, ...
                pupilTimeFile));
            pupilTime = data.eyeTime;
        else
            data = load(fullfile(folderEye, db(k).subject, db(k).date, ...
                num2str(db(k).tlExp), 'eye_timeStamps.mat'));
            pupilTime = data.tVid .* timelineToEphys(1) + timelineToEphys(2);
        end
        nanTimes = pupilTime(isnan(pupil));
        nanTimes = hist(nanTimes, time)>0;
        pupil = interp1(pupilTime, pupil, time, 'pchip');
        pupil(nanTimes) = NaN;
    else
        pupil = NaN(length(time), 1);
    end
    % smooth spike counts, get spike rates
    sr = 1 / diff(time([1 2]));
    winSize = round(windowStd * sr);
    win = normpdf(-5*winSize:5*winSize, 0, winSize);
    indStable = find(corrs(k).isStable);
    spikeRates = zeros(size(spikeCounts,1), length(indStable));
    
    for iCell = 1:size(spikeRates,2)
        spikeRates(:,iCell) = conv(spikeCounts(:,indStable(iCell)), ...
            win, 'same') .* sr;
    end
    avgRates = mean(spikeRates,1);
    [rSigned,indSigned] = sort(corrs(k).(reference).rhos(indStable),'ascend');
    [rAbs,indAbs] = sort(abs(corrs(k).(reference).rhos(indStable)),'ascend');
    sets = cell(1,2);
    ind = find(avgRates(indSigned) > 1, n, 'last');
    sets{1} = [flip(indSigned(ind)), flip(rSigned(ind))];
    ind = find(avgRates(indSigned) > 1, n);
    sets{2} = [indSigned(ind), rSigned(ind)];
%     ind = find(avgRates(indAbs) > 1, n);
%     sets{2} = [indAbs(ind), rAbs(ind)];
    
    if strcmp(reference, 'pupil')
        ind = pupil > nanmedian(pupil);
        nonVisual = pupil;
        nonVisual2 = running;
        label2 = 'running speed';
    else
        ind = running > 5;
        nonVisual = running;
        nonVisual2 = pupil;
        label2 = 'pupil diameter';
    end
    if all(isnan(nonVisual))
        continue
    end
    starts = find(diff(ind)==1);
    stops = find(diff(ind)==-1);
    if starts(1)>stops(1)
        starts = [1; starts];
    end
    if starts(end)>stops(end)
        stops(end+1) = length(nonVisual);
    end
    ind = (stops - starts) >= minSamples;
    starts = starts(ind);
    stops = stops(ind);
    
    dist = 5;
    t = time-time(1);
    ax = zeros(1,3);
    figure('Position', [1 41 1920 1083])
    subplot(10,1,1)
    hold on
    mini = min(nonVisual);
    maxi = max(nonVisual);
    fill(t([starts';stops';stops';starts']), [mini mini maxi maxi]', 'k', ...
        'EdgeColor', 'none', 'FaceColor', [1 1 1].*.9)
    plot(t,nonVisual,'Color', [0 0 0])
    ylabel(reference)
    axis tight
    set(gca,'box','off','XTickLabel',[])
    title('dark')
    ax(1) = gca;
%     title('Correlation with behaviour (gratings and optogenetics)')

    subplot(10,1,2)
    hold on
    mini = min(nonVisual2);
    maxi = max(nonVisual2);
    fill(t([starts';stops';stops';starts']), [mini mini maxi maxi]', 'k', ...
        'EdgeColor', 'none', 'FaceColor', [1 1 1].*.9)
    plot(t,nonVisual2,'Color', [0 0 0])
    ylabel(label2)
    axis tight
    set(gca,'box','off','FontName','Arial','XTickLabel',[])
    ax(2) = gca;
    
    subplot(10,1,3:10)
    hold on
    mini = -(length(sets)*n + (length(sets)-1)*2 - 0.5) * dist;
    maxi = dist;
    fill(t([starts';stops';stops';starts']), [mini mini maxi maxi]', 'k', ...
        'EdgeColor', 'none', 'FaceColor', [1 1 1].*.9)
    y = 0;
    for s = 1:length(sets)
        for j = 1:n
            tr = spikeRates(:,sets{s}(j,1));
            tr = (tr-mean(tr))/std(tr);
            plot(t,tr+y, 'Color', cols(s,:))
            y = y-dist;
        end
        y = y - 2*dist;
    end
    axis tight
    ylim([mini-.2 maxi])
    xlabel('Time (s)')
    ylabel('Firing rate (normalised)')
    ax(3) = gca;
    linkaxes(ax, 'x')
%     set(gca,'FontName','Arial','YTick',[])
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(fullfile(folderResults, ['plots_' reference], 'continuousTraceCorrWithBehaviour', ...
        sprintf('traces_%s_%s_%s.jpg', db(k).subject, db(k).date, prefix)), ...
        '-djpeg','-r0')
    close gcf
    
    % Run previous code cell first for the following code!
    % Cumulative distributions: data vs. null
    x = sort(rhos{f}(~isnan(rhos{f})), 'ascend');
    y = (1:length(x)) ./ length(x);
    xNull = sort(null{f}(~isnan(null{f})), 'ascend');
    yNull = (1:length(xNull))' ./ (length(xNull));
    xNull = [min(x); xNull; max(x)];
    yNull = [0; yNull; 1];
    figure
    hold on
    plot(x, y, 'k', 'LineWidth', 2)
    plot(xNull, yNull, ':', 'Color', [.5 .5 .5], 'LineWidth', 2)
    for s = 1:length(sets)
        heights = interp1(x, y, sets{s}(:,2));
        plot(sets{s}(:,2), heights, 'o', 'MarkerEdgeColor', 'none', ...
            'MarkerFaceColor', cols(s,:))
    end
    xlim(x([1 end]))
    xlabel(['Correlation with ' reference])
    ylabel('Proportion of neurons')
    legend('data', 'null distribution', 'Location', 'NorthWest')
    legend('boxoff')
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(fullfile(folderResults, ['plots_' reference], 'continuousTraceCorrWithBehaviour', ...
        sprintf('distribution_%s_%s_%s.jpg', db(k).subject, db(k).date, prefix)), ...
        '-djpeg','-r0')
    close gcf
    
    % Depth plot
    figure
    hold on
    scatter(rhos{f}(isSgnf{f}==1), depths(isSgnf{f}==1), 'ko', ...
        'MarkerFaceColor', 'k')
    scatter(rhos{f}(isSgnf{f}==0), depths(isSgnf{f}==0), 'ko')
    for s = 1:length(sets)
        scatter(sets{s}(:,2), corrs(k).cellDepths(indStable(sets{s}(:,1))), ...
            'o', 'MarkerFaceColor', cols(s,:), 'MarkerEdgeColor', 'none')
    end
    set(gca, 'YDir','reverse')
    xlim([bins(1)-.05 bins(end)+.05])
    xlabel(['Correlation: firing rate and ' labels{f}])
    ylabel('Depth from surface of SC')
    title([stimulus ': corr. with ' labels{f}])
    legend('pos corr','neg corr')
    legend('boxoff')
end
