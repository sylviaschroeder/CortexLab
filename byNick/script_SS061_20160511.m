ksDir = 'J:\Electrophys\SS061\20160511';
[spikeTimes,spikeAmps,spikeDepths] = ksDriftmap(ksDir);
pars.loadPCs = true;
sp = loadKSdir(ksDir,pars);
cids = sp.cids(sp.cgs==2);
ind = ismember(sp.clu, cids);
cluYpos = clusterAverage(sp.clu(ind), spikeDepths(ind));
[~,order] = sort(cluYpos,'descend');

% load stimulus information and running data

figure('Position',[1921 1 1920 1123])
subplot(7,1,1)
plot(wheelData.t,wheelData.total,'k')
title('Running speed')
ax1=gca;
subplot(7,1,2:7)
colors = prism(7);
nColBins = 10;
hold on
for iCell=1:length(order)
    id = cids(order(iCell));
    ind = sp.clu == id;
    ampRange = quantile(spikeAmps(ind),[.1 .9]);
    colorBins = [min(spikeAmps(ind))-1, linspace(ampRange(1), ampRange(2), ...
        nColBins-1), max(spikeAmps(ind))];
    grayness = linspace(0,0.9,nColBins);
    hueInd = mod(iCell,size(colors,1));
    if hueInd == 0
        hueInd = size(colors,1);
    end
    cols = bsxfun(@times,colors(hueInd,:),(1-grayness)') + ...
        bsxfun(@times,[1 1 1].*0.5,grayness');
    for b = 1:nColBins-1
        s = sp.clu==cids(iCell) & spikeAmps>colorBins(b) & spikeAmps<=colorBins(b+1);
        if sum(s)>0
            plot(spikeTimes(s), spikeDepths(s), '.', 'Color', cols(b,:));
        end
    end
end
xlabel('time');
ylabel('y position');
colormap prism
colorbar('Limits',[0 0.093],'Ticks',[],'Position',[.932 .111 .02 .691])
linkaxes([ax1 gca],'x')
xlim([min(spikeTimes),max(spikeTimes)])


%% load spike data

Fs = 30000;
clu = readNPY('spike_clusters.npy');
ss = readNPY('spike_times.npy');
st = double(ss)/Fs;
spikeTemplates = readNPY('spike_templates.npy'); % note: zero-indexed
% clu = spikeTemplates;
tempScalingAmps = readNPY('amplitudes.npy');

[cids, cgs] = readClusterGroupsCSV('cluster_groups.csv');

% noiseClusters = cids(cgs==0);
% st = st(~ismember(clu, noiseClusters));
% spikeTemplates = spikeTemplates(~ismember(clu, noiseClusters));
% tempScalingAmps = tempScalingAmps(~ismember(clu, noiseClusters));
% clu = clu(~ismember(clu, noiseClusters));
% cgs = cgs(~ismember(cids, noiseClusters));
% cids = cids(~ismember(cids, noiseClusters));

goodClusters = cids(cgs==2);
st = st(ismember(clu, goodClusters));
spikeTemplates = spikeTemplates(ismember(clu, goodClusters));
tempScalingAmps = tempScalingAmps(ismember(clu, goodClusters));
clu = clu(ismember(clu, goodClusters));
cgs = cgs(ismember(cids, goodClusters));
% cids = cids(ismember(cids, goodClusters));

%% compute template and spike Amplitudes and Depths
load forPRBimecP3opt3
yc = ycoords(connected); xc = xcoords(connected);
temps = readNPY('templates.npy');

winv = readNPY('whitening_mat_inv.npy');

[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW] = ...
    templatePositionsAmplitudes(temps, winv, yc, spikeTemplates, tempScalingAmps);

%% get Ypos and Amp for each cluster
cids = unique(clu);

% without manual sorting
% cluYpos = templateYpos(cids+1);
% cluAmps = tempAmps(cids+1);

% with manual sorting
cluYpos = clusterAverage(clu, spikeDepths);
cluAmps = clusterAverage(clu, spikeAmps);

%% plot y-pos, time and amplitude of each spike
figure
a = spikeAmps-min(spikeAmps);
a = a ./ max(a) .* 20 + .1;
c = spikeAmps-min(spikeAmps);
c = 1 - (c ./ max(c) .* .8 + .2);
scatter(st, spikeDepths, a);

%% compute firing rates

cids = unique(clu);
recDur = st(end);
FRs = zeros(size(cids));
for c = 1:length(cids);
    FRs(c) = sum(clu==cids(c))./recDur;
end

%% plot templates by depth

yy = cluYpos; xx = rand(size(yy))*max(yy);

figure;

subplot(2,2,1);

scatter(xx,yy,cluAmps);
title('amplitudes');

subplot(2,2,3);
hist(cluAmps, 25);

subplot(2,2,2);

scatter(xx,yy,FRs);
title('firing rates');

subplot(2,2,4);
hist(FRs, 0:60);

%% stimulus info for gratings

fid = fopen('SS061_20160511_g0_t0.imec_sync.bin');
syncDat = fread(fid, Inf, '*int16');
fclose(fid);

syncFs = 3000;
tSync = (0:length(syncDat)-1)/syncFs;

flipsDown = tSync(syncDat(1:end-1)==-1 & syncDat(2:end)==-2);
flipsUp = tSync(syncDat(1:end-1)==-2 & syncDat(2:end)==-1);


stimStarts_1 = flipsDown(diff([0 flipsDown])>0.2 & flipsDown<2700); %  empirically determined transition time
stimStops_1 = stimStarts_1+2;
load('\\ZSERVER\Data\trodes\M160426_SS061\20160511\1\Protocol.mat')
stimIDs_1 = zeros(size(stimStarts_1));
for q = 1:size(Protocol.seqnums,1)
    stimIDs_1(Protocol.seqnums(q,:)) = q;
end

stimSequence = ppbox.getStimSequence('M160426_SS061',20160511,1);
[directions, blanks] = gratings.getOrientations(stimSequence,'ori1','c1');

spikeRes = 0.005;
timeBeforeAfterStimuli = 10; %in sec

% t = round((stimStarts_1(1)-3)/spikeRes)*spikeRes : spikeRes : stimStops_1(end)+3;
t = round((stimStarts_1(1)-timeBeforeAfterStimuli)/spikeRes)*spikeRes : ...
    spikeRes : (stimStops_1(end)+timeBeforeAfterStimuli);
stimTimes.onset = stimStarts_1;
stimTimes.offset = stimStops_1;
stimMatrix = ppbox.buildStimMatrix(stimSequence, stimTimes, t);

laserOn = [ones(13,1); zeros(13,1)];

%% running data
runWinStd = 0.5; % in sec
% load Timeline
load('\\ZSERVER\Data\expInfo\SS061\2016-05-11\1\2016-05-11_1_SS061_Timeline.mat')

% get conversion from Timeline time into spike data time
phd=Timeline.rawDAQData(:, 2);
phd=(phd-min(phd))/(max(phd)-min(phd));
thr=0.5; % using one threshold here
above=phd>thr;
deltas=[0; diff(above)];
goingDownTimes=Timeline.rawDAQTimestamps(deltas==-1);
starts_1 = goingDownTimes(diff([0 goingDownTimes])>0.2 & goingDownTimes<2700);
b = regress(stimStarts_1', [starts_1',ones(length(starts_1),1)]);

wheel = double(Timeline.rawDAQData(:,3));
tlTime = Timeline.rawDAQTimestamps .* b(1) + b(2);
wheelData = nonVis.getRunningSpeed_wheel(wheel, tlTime, runWinStd);

%% Prepare data used for all analyses of this dataset
yposWin = [640 2390]; %[1500 2390];
% time vector
time = t';
% spike rates
cluInds = cluYpos>yposWin(1) & cluYpos<yposWin(2); % & FRs>.1;
inclClu = cids(cluInds);
spikeCounts = zeros(length(t), length(inclClu));
for iCell = 1:length(inclClu)
    spikeInds = clu==inclClu(iCell) & st>=time(1) & st<=time(end);
    sc = hist(st(spikeInds), t);
    spikeCounts(:,iCell) = sc;
end
indResp = sum(spikeCounts,1) > 0;
spikeCounts = spikeCounts(:,indResp);
% cell IDS
cellIDs = inclClu(indResp);
% depths
depths = cluYpos(cluInds);
depths = depths(indResp);
% running speed
running = interp1(wheelData.t, wheelData.total, t, 'pchip')';

save(fullfile('C:\DATA\electrophys\SS061\20160511', 'data.mat'), ...
    'time', 'spikeRates', 'cellIDs', 'running', 'stimMatrix', 'directions', ...
    'blanks', 'stimTimes', 'laserOn');

%% PSTHs 

ampThresh = 2.5;
% yposWin = [2100 2400];
% yposWin = [1000 2100];
yposWin = [1500 2390];
inclClu = cids(cluAmps>ampThresh & cluYpos>yposWin(1) & cluYpos<yposWin(2) & FRs>0.5);
tempIDs = stimIDs_1;
tempIDs(tempIDs==13) = 27;
psthViewer(st(ismember(clu, inclClu)), clu(ismember(clu, inclClu)), stimStarts_1, [-0.5 3.0], tempIDs);

%% Plot correlation with running
spikeWinStd = 0.1; % in sec
spikeRes = 0.005; % in sec
filtPoly = 3;
smoothing = 3;
ampThresh = 2.5;

% yposWin = [2100 2400]; % intermediate layers
yposWin = [1500 2390]; %[1500 2390]; % visually responsive part of SC (estimated during recordings)
% cluInds = cluAmps>ampThresh & cluYpos>yposWin(1) & cluYpos<yposWin(2) & FRs>0.5;
cluInds = cluYpos>yposWin(1) & cluYpos<yposWin(2);
inclClu = cids(cluInds);
% smooth all spikes from selected neurons and plot together with running
% trace
% (1) grating experiment
t = round(stimStarts_1(1)/spikeRes)*spikeRes : spikeRes : stimStops_1(end);
spikeWinSamples = spikeWinStd / spikeRes;
win = pdf('Normal', -5*spikeWinSamples:5*spikeWinSamples, 0, spikeWinSamples);
filtWindow = ceil(smoothing / spikeRes);
if mod(filtWindow,2) == 0
    filtWindow = filtWindow-1;
end

spikeInds = ismember(clu, inclClu) & st>=stimStarts_1(1) & st<=stimStops_1(end);
spikeCount = hist(st(spikeInds), t);
rate = conv(spikeCount, win, 'same');

wheel = interp1(wheelData.t, wheelData.total, t, 'pchip');

figure
plot(t, wheel)
ax1 = gca;
figure
plot(t, rate-max(rate))
ax2 = gca;
linkaxes([ax1 ax2], 'x')

% for each neuron
spikeCounts = zeros(length(t), length(inclClu));
for iCell = 1:length(inclClu)
    spikeInds = clu==inclClu(iCell) & st>=stimStarts_1(1) & st<=stimStops_1(end);
    sc = hist(st(spikeInds), t);
%     spikeRates(:,iCell) = conv(sc, win, 'same');
    spikeCounts(:,iCell) = sc;
end
spikeCounts = sgolayfilt(spikeCounts, filtPoly, filtWindow);

[runResults, figHandles] = nonVis.getCorrToNonVisData(spikeCounts, ...
    t, wheel, t, 'Running Speed', inclClu, -1, 1, stimStarts_1);

% plot y-pos vs. corr. with running
rho = runResults.rho;
ind = ~isnan(rho);
rho = rho(ind);
pos = cluYpos(cluInds);
pos = pos(runResults.order(ind));
amp = cluAmps(cluInds);
amp = amp(runResults.order(ind));
binSize = .05;
colBinCenters = -1:binSize:1;
colBinEdges = [colBinCenters-binSize/2, colBinCenters(end)+binSize/2];
[~,~,binPerCell] = histcounts(rho,colBinEdges);
colMap = jet(length(colBinCenters));
cols = colMap(binPerCell,:);
figure
scatter(abs(rho), pos, amp, cols, 'LineWidth', 2)
colormap jet
b = colorbar('Ticks',0:.1:1,'TickLabels',-1:.2:1);
b.Label.String = 'Corr. w running';
xlabel('Corr. with running (absolute value)')
ylabel('Depth (microns)')
set(gca,'YTick',1500:300:2400,'YTickLabel',900:-300:0)

%% Plot tuning curves (single neurons): V1 inactivation and running

runningThreshold = 0;

yposWin = [1500 2390]; %[1500 2390]; % visually responsive part of SC (estimated during recordings)
inclClu = cids(cluYpos>yposWin(1) & cluYpos<yposWin(2));

wheel = interp1(wheelData.t, wheelData.total, t);
isRunning = double(wheel > runningThreshold);
runPerTrial = ssLocal.getTracesPerStimulus(isRunning', stimMatrix, [0 0]);
running = squeeze(mean(runPerTrial,4)); % [stimuli x repetition]
running = double(running >= 0.4);

gratings.tuningCurveViewer(st(ismember(clu, inclClu)), clu(ismember(clu, inclClu)), ...
    cids, cluYpos, stimMatrix, t, directions, blanks, laserOn, running);

%% Make scatter plot matrix and grand average tuning curves

runningThreshold = 0;
stimOffsets = [0 0];
yposWin = [1500 2400]; %[1500 2390]; % visually responsive part of SC (estimated during recordings)

offsetSamples = round(stimOffsets ./ spikeRes);
inclClu = cids(cluYpos>yposWin(1) & cluYpos<yposWin(2));

wheel = interp1(wheelData.t, wheelData.total, t);
isRunning = double(wheel > runningThreshold);
runPerTrial = ssLocal.getTracesPerStimulus(isRunning', stimMatrix, [0 0]);
running = squeeze(mean(runPerTrial,4)); % [stimuli x repetition]
running = double(running >= 0.4);

conditions = bsxfun(@plus, laserOn .* 2, running) + 1;
%1: laser off & not running
%2: laser off & running
%3: laser on & not running
%4: laser on & running

respSummary = NaN(length(inclClu), 4, length(unique(conditions(:))));
% 1st dim.: neurons
% 2nd dim.: response at (1) baseline, (2) pref. dir., (3) 90 deg. from 
% pref., (4) 180 deg. from pref.
% 3rd dim.: conditions (see 'conditions' above)
labelsResponses = {'basel.','pref. dir.','orthog. dir.','oppos. dir.'};
labelsConditions = {'laser off, no running','laser off, running', ...
    'laser on, no running','laser on, running'};
R2s = NaN(length(inclClu),1);

x = 0:359;
tuningCurves = NaN(length(inclClu), length(x), length(unique(conditions(:))));
blankResponses = NaN(length(inclClu), length(unique(conditions(:))));

timePlus = [t(1)-spikeRes t t(end)+spikeRes];
for iCell = 1:length(inclClu)
    nchars = fprintf('Neuron %d of %d', iCell, length(inclClu));
    iSt = st(clu==inclClu(iCell));
    if isempty(iSt)
        continue
    end
    spikeCount = hist(iSt, timePlus);
    spikeCount([1 end]) = [];
    respPerTrial = ssLocal.getTracesPerStimulus(spikeCount(:), ...
        stimMatrix, offsetSamples);
    response = squeeze(mean(respPerTrial,4)) ./ spikeRes; % [stimulus x repetition]
    if nansum(response(:)) == 0
        continue
    end
    
    [parameters, blankResp, predictions] = gratings.fitTuningCurveConditions( ...
        response, directions, blanks, conditions, [1 3 5]);
    for c = 1:length(unique(conditions(:)))
        prefDirs = parameters(1,c) + [0 90 180];
        respPref = orituneWrappedConditions(parameters(:,c),prefDirs);
        respSummary(iCell,:,c) = [blankResp(c), respPref];
        
        tCurve = orituneWrappedConditions(parameters(:,c),x);
        tuningCurves(iCell,:,c) = tCurve;
    end
    blankResponses(iCell,:) = blankResp;
    
    resp = response;
    resp(blanks,:) = [];
    pred = predictions;
    pred(blanks,:) = [];
    R2 = 1 - nansum((resp(:)-pred(:)).^2) / ...
        nansum((resp(:)-blankResp(1)).^2);
    R2s(iCell) = R2;
    fprintf(repmat('\b', 1, nchars));
end

% Plot average tuning curves for various conditions
% 1. align tuning curves to pref. direction + divide by max. response when
% laser off and not running
xNorm = -180:180;
tcNorm = NaN(size(tuningCurves,1),length(xNorm),size(tuningCurves,3));
brNorm = NaN(size(blankResponses));
for iCell = 1:size(tuningCurves,1)
    if R2s(iCell)<.2
        continue
    end
    [maxResp,prefInd] = max(squeeze(tuningCurves(iCell,:,1)));
    inds = mod(prefInd + xNorm, 360);
    inds(inds == 0) = 360;
    for c = 1:size(tuningCurves,3)
        tcNorm(iCell,:,c) = tuningCurves(iCell,inds,c) ./ maxResp;
    end
    brNorm(iCell,:) = blankResponses(iCell,:) ./ maxResp;
end
% 2. plot aligned and normalized tuning curves
colors = [0 0 0; 1 0 0; 0 0 0; 1 0 0];
lins = {'-','-','--','--'};
minR2 = 0.2;
% only consider certain neurons
% subsets = {R2s>=minR2 & respSummary(:,2,1)>respSummary(:,1,1) & ...
%     respSummary(:,2,2)>respSummary(:,2,1); ...
%     R2s>=minR2 & respSummary(:,2,1)>respSummary(:,1,1) & ...
%     respSummary(:,2,2)<respSummary(:,2,1);
%     R2s>=minR2 & respSummary(:,2,1)<respSummary(:,1,1) & ...
%     respSummary(:,2,2)>respSummary(:,2,1);
%     R2s>=minR2 & respSummary(:,2,1)<respSummary(:,1,1) & ...
%     respSummary(:,2,2)<respSummary(:,2,1)};
% descrips = {'Contrast+ & Running+', 'Contrast+ & Running-', ...
%     'Cotrast- & Running+', 'Contrast- & Running-'};
subsets = {R2s>=minR2 & respSummary(:,2,2)>respSummary(:,2,1); ...
    R2s>=minR2 & respSummary(:,2,2)<respSummary(:,2,1)};
descrips = {'Running+', 'Running-'};
for s = 1:length(subsets)
    f = figure;
    hold on
    h = zeros(1,size(tcNorm,3));
%     maxi = max([reshape(tcNorm(subsets{s},:,:),[],1); ...
%         reshape(brNorm(subsets{s},:),[],1)]);
%     mini = min([reshape(tcNorm(subsets{s},:,:),[],1); ...
%         reshape(brNorm(subsets{s},:),[],1)]);
    maxi = max([reshape(tcNorm(subsets{s},:,:),[],1)]);
    mini = min([reshape(tcNorm(subsets{s},:,:),[],1)]);
    for c = 1:size(tcNorm,3)
%         figure
%         plot(xNorm,tcNorm(subsets{s},:,c)')
%         hold on
%         plot(210,brNorm(subsets{s},c),'o')
%         xlim([-180 220])
%         ylim([mini-.05*(maxi-mini), maxi+.05*(maxi-mini)])
%         set(gca,'XTick',[-180:90:180,210],'XTickLabel',{-180:90:180,'blank'})
%         xlabel('Direction')
%         ylabel('Normalized response')
%         title(labelsConditions{c})
        
        ind = all(~isnan(tcNorm(:,:,c)),2)&all(~isinf(tcNorm(:,:,c)),2) & ...
            subsets{s};
        m = smooth(median(tcNorm(ind,:,c),1),30,'lowess')';
        se = std(tcNorm(ind,:,c),0,1) / sqrt(sum(ind));
        figure(f)
        fill(xNorm([1:end,end:-1:1]),[m+se,flip(m)-flip(se)],'k', ...
            'EdgeColor','none','FaceColor',colors(c,:),'FaceAlpha',.2)
        h(c) = plot(xNorm,m,lins{c},'Color',colors(c,:),'LineWidth',2);
%         m = median(brNorm(ind,c));
%         se = std(brNorm(ind,c)) / sqrt(sum(ind));
%         errorbar(200+c*10,m,se,'o','Color',colors(c,:),'LineWidth',2)
%         set(gca,'XTick',[-180:90:180,225],'XTickLabel',{-180:90:180,'blank'})
        xlabel('Direction')
        ylabel('Normalized response')
        title(sprintf('%s (n=%d)',descrips{s},sum(ind)))
    end
    legend(h,labelsConditions)
end

% Plot all variables against each other
valid = R2s >= 0.3;
for r1 = 1:4
    for c1 = 1:4
        for r2 = r1+1:4
            figure
            scatter(respSummary(valid,r1,c1), respSummary(valid,r2,c1))
            xlabel(labelsResponses{r1})
            ylabel(labelsResponses{r2})
            title(labelsConditions{c1})
            mini = min(reshape(respSummary(valid,[r1 r2],c1),[],1));
            maxi = max(reshape(respSummary(valid,[r1 r2],c1),[],1));
            range = maxi-mini;
            mini = mini-.05*range;
            maxi = maxi+.05*range;
            hold on
            plot([mini maxi],[mini maxi])
            axis([mini maxi mini maxi])
            axis square
        end
        for c2 = c1+1:4
            figure
            scatter(respSummary(valid,r1,c1), respSummary(valid,r1,c2))
            xlabel(labelsConditions{c1})
            ylabel(labelsConditions{c2})
            title(labelsResponses{r1})
            mini = min(reshape(respSummary(valid,r1,[c1 c2]),[],1));
            maxi = max(reshape(respSummary(valid,r1,[c1 c2]),[],1));
            range = maxi-mini;
            mini = mini-.05*range;
            maxi = maxi+.05*range;
            hold on
            plot([mini maxi],[mini maxi])
            axis([mini maxi mini maxi])
            axis square
        end
    end
end

% 1st: response difference: pref. dir. - opp. dir.
% 2nd: pref. dir. - orth. dir.s
% 3rd: pref. dir. - min(opp. dir., orth. dirs)
diffs{1} = squeeze(respSummary(:,2,:) - respSummary(:,4,:));
diffs{2} = squeeze(respSummary(:,2,:) - respSummary(:,3,:));
diffs{3} = squeeze(respSummary(:,2,:) - min(respSummary(:,[3 4],:),[],2));
labelsDiffs = {'Pref. - opp.', 'Pref. - orth', 'Pref. - min(opp.,orth.)'};
condPairs = [1 2; 3 4; 1 3; 2 4];
valid = R2s >= 0.3;
for d = 1:length(diffs)
    for c = 1:size(condPairs,1)
        figure
        scatter(diffs{d}(valid,condPairs(c,1)), diffs{d}(valid,condPairs(c,2)))
        xlabel(labelsConditions{condPairs(c,1)})
        ylabel(labelsConditions{condPairs(c,2)})
        title(labelsDiffs{d})
        mini = min(reshape(diffs{d}(valid,condPairs(c,:)),[],1));
        maxi = max(reshape(diffs{d}(valid,condPairs(c,:)),[],1));
        range = maxi-mini;
        mini = mini-.05*range;
        maxi = maxi+.05*range;
        hold on
        plot([mini maxi],[mini maxi])
        axis([mini maxi mini maxi])
        axis square
    end
end

% 1st: response difference: pref. during running - pref. when not running
% 2nd: response ratio: pref. during running / pref. when not running
diffs = cell(1,4);
diffs{1} = squeeze(respSummary(:,2,[2 4]) - respSummary(:,2,[1 3]));
diffs{2} = squeeze(respSummary(:,2,[2 4]) ./ respSummary(:,2,[1 3]));
diffs{3} = squeeze(respSummary(:,1,[2 4]) - respSummary(:,1,[1 3]));
diffs{4} = squeeze(respSummary(:,1,[2 4]) ./ respSummary(:,1,[1 3]));
labelsDiffs = {'Pref.: running - not running', 'Pref.: running / not running', ...
    'Baseline: running - not running', 'Baseline: running / not running'};
valid = R2s >= 0.3
for d = 1:length(diffs)
    figure
    scatter(diffs{d}(valid,1), diffs{d}(valid,2))
    xlabel('laser off')
    ylabel('laser on')
    title(labelsDiffs{d})
    tmp = reshape(diffs{d}(valid,:),[],1);
    tmp(isinf(tmp)) = [];
    mini = min(tmp);
    maxi = max(tmp);
    range = maxi-mini;
    mini = mini-.05*range;
    maxi = maxi+.05*range;
    hold on
    if mod(d,2)==1
        plot([0 0],[mini maxi],'k:')
        plot([mini maxi],[0 0],'k:')
    else
        plot([1 1],[mini maxi],'k:')
        plot([mini maxi],[1 1],'k:')
    end
    axis([mini maxi mini maxi])
    axis square
end

% compare running effect across laser conditions
laserOff = (respSummary(valid,2,2)-min(respSummary(valid,[3 4],2),[],2)) ./ ...
    (respSummary(valid,2,1)-min(respSummary(valid,[3 4],1),[],2));
laserOn = (respSummary(valid,2,4)-min(respSummary(valid,[3 4],4),[],2)) ./ ...
    (respSummary(valid,2,3)-min(respSummary(valid,[3 4],3),[],2));
ind = any([laserOff, laserOn] > 30, 2);
laserOff(ind) = [];
laserOn(ind) = [];
figure
scatter(laserOff, laserOn)
xlabel('Laser off')
ylabel('Laser on')
title('pref-min: run / not run')
mini = min([laserOff; laserOn]);
maxi = max([laserOff; laserOn]);
range = maxi-mini;
mini = mini-.05*range;
maxi = maxi+.05*range;
hold on
plot([mini maxi],[mini maxi])
axis([mini maxi mini maxi])
axis square

laserOff = (respSummary(valid,2,2)-min(respSummary(valid,[3 4],2),[],2)) - ...
    (respSummary(valid,2,1)-min(respSummary(valid,[3 4],1),[],2));
laserOn = (respSummary(valid,2,4)-min(respSummary(valid,[3 4],4),[],2)) - ...
    (respSummary(valid,2,3)-min(respSummary(valid,[3 4],3),[],2));
figure
scatter(laserOff, laserOn)
xlabel('Laser off')
ylabel('Laser on')
title('pref-min: run - not run')
mini = min([laserOff; laserOn]);
maxi = max([laserOff; laserOn]);
range = maxi-mini;
mini = mini-.05*range;
maxi = maxi+.05*range;
hold on
plot([0 0],[mini maxi],'k:')
plot([mini maxi],[0 0],'k:')
axis([mini maxi mini maxi])
axis square

% normalize responses
respSumNorm = bsxfun(@minus, respSummary, min(respSummary(:,[3 4],1),[],2));
respSumNorm = bsxfun(@rdivide, respSumNorm, respSumNorm(:,2,1));
figure
boxplot([squeeze(respSumNorm(valid,2,[1 2 3 4])), ...
    squeeze(min(respSumNorm(valid,[3 4],:),[],2))])
ylim([-1.5 4.5])
set(gca,'XTickLabel',{'Pref,L off,R off','Pref,L off,R on','Pref,L on,R off',...
    'Pref,L on,R on','Min,L off,R off','Min,L off,R on','Min,L on,R off',...
    'Min,L on,R on'},'XTickLabelRotation',45);

%% Plot modulation index due to inactivation
% exclude low amplitude and low firing rate neurons
minRate = 0.5;
durTrial = 2;
tempWindow = [0.4 1.7];
laserOn_1 = [ones(13,1); zeros(13,1)];
ampThresh = 2.5;

clu_1 = clu(st>=stimStarts_1(1) & st<=stimStops_1(end));
dur_1 = stimStops_1(end) - stimStarts_1(1);
[cells_1, spikeC_1] = countUnique(clu_1);
cells_1(spikeC_1/dur_1 < minRate) = [];
firingRateInds_1 = ismember(cids, cells_1);
goodNeurons_1 = cids(firingRateInds_1 & cluAmps>ampThresh);

modIndices_1 = ephys.getInactivationIndices(st(ismember(clu, goodNeurons_1)), ...
    clu(ismember(clu, goodNeurons_1)), stimStarts_1, tempWindow, ...
    stimIDs_1, laserOn_1);
modIndices_1(modIndices_1<-3) = -3;
modIndices_1(modIndices_1>3) = 3;
figure, hist(modIndices_1)
xlabel('Modulation index')
ylabel('# Units')
colMap = jet(201);
valid = ~isnan(modIndices_1);
cols = zeros(length(modIndices_1),3);
cols(valid,:) = colMap(round(modIndices_1(valid)*33)+101,:);
% areas = round(abs(modIndices_1)*10);
figure
scatter(xx(ismember(cids, goodNeurons_1(valid))), ...
    yy(ismember(cids, goodNeurons_1(valid))), 50, cols(valid,:))
set(gca, 'Color', [0.5 0.5 0.5])
colormap jet
colorbar('Ticks',0:0.1:1, 'TickLabels', -3:3/5:3)
ylim([0 max(yy)])
ylabel('Depth in microns')
title('Modulation of SC (and others) due to inactivation of V1')

%% load stimulus info for sparse noise

load('\\ZSERVER\Data\trodes\M160426_SS061\20160511\2\Protocol.mat')

fid = fopen('SS061_20160511_g0_t0.imec_sync.bin');
syncDat = fread(fid, Inf, '*int16');
fclose(fid);

syncFs = 3000;
tSync = (0:length(syncDat)-1)/syncFs;

flipsDown = tSync(syncDat(1:end-1)==-1 & syncDat(2:end)==-2);
flipsUp = tSync(syncDat(1:end-1)==-2 & syncDat(2:end)==-1);

flipsDown = flipsDown(flipsDown>2835 & flipsDown<3763); % empirical
flipsUp = flipsUp(flipsUp>2835 & flipsUp<3763);

photodiodeFlips = sort([flipsUp flipsDown]);

%% compute sparse noise frames

load(dat.expFilePath('M160426_SS061', '2016-05-11', 2, 'hw-info', 'master'));
myScreenInfo.windowPtr = NaN;
% [allFrames, frameTimes] = computeSparseNoiseFrames(Protocol, photodiodeFlips, myScreenInfo);
stimNum = 1;
ss = eval([Protocol.xfile(1:end-2) '(myScreenInfo, Protocol.pars(:,stimNum));']);

% convert stimulus info to array
imTextSeq = ss.ImageTextures(ss.ImageSequence(1:end-1)); % excluding the last ImageSequence here as a hack to make the right number of photodiode events (?)
q = reshape([imTextSeq{:}], size(ss.ImageTextures{1},1), size(ss.ImageTextures{1},2), []);
stimArray = repmat(q, [1 1 Protocol.nrepeats]); clear q; 

nX = size(stimArray,1);
nY = size(stimArray,2);
stimArrayZeroPad = cat(3,zeros(size(stimArray,1), size(stimArray,2),1), stimArray);
for x = 1:nX
    for y = 1:nY
        stimEventTimes{x,y,1} = photodiodeFlips(stimArrayZeroPad(x,y,1:end-1)==0 & ...
            stimArrayZeroPad(x,y,2:end)==1); % going from grey to white
        stimEventTimes{x,y,2} = photodiodeFlips(stimArrayZeroPad(x,y,1:end-1)==0 & ...
            stimArrayZeroPad(x,y,2:end)==-1); % going from grey to black
    end
end

%% compute RF for some cluster

stimEventsToUse = 2; % 1 is grey-to-white, 2 is grey-to-black. 
spikeCountingWindow = [0 0.1]; % sec

useCIDs = [292 878 344 375 508 570 585 599 621 754 842];

for c = 1:length(useCIDs)
    thisCID = useCIDs(c);
    theseST = st(clu==thisCID);
    
    thisRF = zeros(nX,nY);
    
    for x = 1:nX
        for y = 1:nY            
            [psth, bins, rasterX, rasterY, spikeCounts] = ...
                psthRasterAndCounts(theseST, stimEventTimes{x,y,stimEventsToUse}(:), ...
                spikeCountingWindow, 0.001);
            thisRF(x,y) = mean(spikeCounts);
        end
    end
    
    figure;
    subplot(2,1,1)
    imagesc(thisRF)
    title(sprintf('cluster %d', thisCID));
    
    subplot(2,1,2);
    imagesc(conv2(thisRF, ones(3,3)/9, 'same'));
    
    drawnow
end
