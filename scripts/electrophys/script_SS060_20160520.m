%% Load spike data
ksDir = 'J:\Electrophys\SS060\20160520';
pars.loadPCs = true;
sp = loadKSdir(ksDir,pars);
Fs = sp.sample_rate;
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

%% Plot drift map
ksDir = 'J:\Electrophys\SS060_20160520';
[spikeTimes,spikeAmps,spikeDepths] = ksDriftmap(ksDir);
pars.loadPCs = true;
sp = loadKSdir(ksDir,pars);
cids = sp.cids(sp.cgs==2);

figure('Position',[1921 1 1920 1123])
inclu = ismember(sp.clu,cids);
plotDriftmap(spikeTimes(inclu), spikeAmps(inclu), spikeDepths(inclu))
for iCell=1:length(cids)
    nColBins = 10;
    ind = sp.clu == cids(iCell);
    ampRange = quantile(spikeAmps(ind),[.1 .9]);
    colorBins = [min(spikeAmps(ind))-1, linspace(ampRange(1), ampRange(2), ...
        nColBins-1), max(spikeAmps(ind))];
    redness = linspace(.1,1,nColBins);
    cols = bsxfun(@times,redness',[1 0 0]) + bsxfun(@times,1-redness',[1 1 1]);
    h = zeros(1,nColBins);
    for b = 1:nColBins-1
        s = sp.clu==cids(iCell) & spikeAmps>colorBins(b) & spikeAmps<=colorBins(b+1);
        if sum(s)>0
            h(b) = plot(spikeTimes(s), spikeDepths(s), '.', 'Color', cols(b,:));
        end
    end
    xlabel('time');
    ylabel('y position');
    
    title(sprintf('Cell cluster %d', cids(iCell)))
    pause
    delete(h)
end

%% compute template and spike Amplitudes and Depths
yc = sp.ycoords;
xc = sp.xcoords;
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
hist(FRs, 25);

%% stimulus info for gratings

fid = fopen(fullfile(ksDir, 'SS060_20160520_SC_g0_t0.imec_sync.bin'));
syncDat = fread(fid, Inf, '*int16');
fclose(fid);

syncFs = 3000;
tSync = (0:length(syncDat)-1)/syncFs;

flipsDown = tSync(syncDat(1:end-1)==-1 & syncDat(2:end)==-2);
flipsUp = tSync(syncDat(1:end-1)==-2 & syncDat(2:end)==-1);

stimStarts_1 = flipsDown(diff([0 flipsDown])>0.2 & flipsDown<3100); %  empirically determined transition time
stimStops_1 = stimStarts_1+2;

load('\\ZSERVER\Data\trodes\M160425_SS060\20160520\1\Protocol.mat')
stimIDs_1 = zeros(size(stimStarts_1));
for q = 1:size(Protocol.seqnums,1)
    stimIDs_1(Protocol.seqnums(q,:)) = q;
end

stimSequence = ppbox.getStimSequence('M160425_SS060',20160520,1);
[directions, blanks] = gratings.getOrientations(stimSequence,'ori1','c1');

spikeRes = 0.005;
timeBeforeAfterStimuli = 10; %in sec

t = round((stimStarts_1(1)-timeBeforeAfterStimuli)/spikeRes)*spikeRes : ...
    spikeRes : (stimStops_1(end)+timeBeforeAfterStimuli);
stimTimes.onset = stimStarts_1;
stimTimes.offset = stimStops_1;
stimMatrix = ppbox.buildStimMatrix(stimSequence, stimTimes, t);

laserOn = [ones(13,1); zeros(13,1)];

%% running data
runWinStd = 0.5; % in sec
% load Timeline (Timeline crashed, therefore we need to read the binary
% file)
fid = fopen('\\ZSERVER\Data\expInfo\SS060\2016-05-20\98\2016-05-20_98_SS060_Timeline.dat', 'r');
dat = fread(fid, [13, Inf], '*double')';
fclose(fid);
fs = 2500;
datTime = (1:size(dat,1)) ./ fs;

% get conversation from Timeline time into spike data time
phd=dat(:, 2);
phd=(phd-min(phd))/(max(phd)-min(phd));
thr=0.5; % using one threshold here
above=phd>thr;
deltas=[0; diff(above)];
goingDownTimes=datTime(deltas==-1);
starts_1 = goingDownTimes(diff([0 goingDownTimes])>0.2 & goingDownTimes<2860);
b = regress(stimStarts_1', [starts_1',ones(length(starts_1),1)]);

wheel = double(dat(:,3));
tlTime = datTime .* b(1) + b(2);
wheelData = nonVis.getRunningSpeed_wheel(wheel, tlTime, runWinStd);

%% Prepare data used for all analyses of this dataset
yposWin = [850 2100];
% time vector
time = t'; % (includes 10 sec before and after stimuli)
% spike rates
cluInds = cluYpos>yposWin(1) & cluYpos<yposWin(2);
inclClu = cids(cluInds);
spikeRates = zeros(length(t), length(inclClu));
for iCell = 1:length(inclClu)
    spikeInds = clu==inclClu(iCell) & st>=stimStarts_1(1) & st<=stimStops_1(end);
    sc = hist(st(spikeInds), t);
    spikeRates(:,iCell) = sc;
end
% cell IDS
cellIDs = inclClu;
% running speed
running = interp1(wheelData.t, wheelData.total, t, 'pchip')';

save(fullfile('C:\DATA\electrophys\SS060\20160520', 'data.mat'), ...
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
yposWin = [0 2390]; %[1500 2390]; % visually responsive part of SC (estimated during recordings)
cluInds = cluAmps>ampThresh & cluYpos>yposWin(1) & cluYpos<yposWin(2) & FRs>0.5;
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

wheel = interp1(wheelData.t, wheelData.total, t);

figure
plot(t, wheel)
ax1 = gca;
figure
plot(t, rate-max(rate))
ax2 = gca;
linkaxes([ax1 ax2], 'x')

% for each neuron
spikeRates = zeros(length(t), length(inclClu));
for iCell = 1:length(inclClu)
    spikeInds = clu==inclClu(iCell) & st>=stimStarts_1(1) & st<=stimStops_1(end);
    sc = hist(st(spikeInds), t);
    spikeRates(:,iCell) = conv(sc, win, 'same');
end
spikeRates = sgolayfilt(spikeRates, filtPoly, filtWindow);

[runResults, figHandles] = nonVis.getCorrToNonVisData(spikeRates, ...
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
ylabel('Depth')

%% Plot tuning curves: V1 inactivation and running

runningThreshold = 0;
spikeRes = 0.005;

t = round((stimStarts_1(1)-3)/spikeRes)*spikeRes : spikeRes : stimStops_1(end)+3;

yposWin = [0 2400]; %[1500 2390]; % visually responsive part of SC (estimated during recordings)
inclClu = cids(cluYpos>yposWin(1) & cluYpos<yposWin(2));

stimTimes.onset = stimStarts_1;
stimTimes.offset = stimStops_1;
stimSequence.seq = stimIDs_1;
stimSequence.labels = unique(stimIDs_1);
stimMatrix = ppbox.buildStimMatrix(stimSequence, stimTimes, t);

laserOn = [ones(13,1); zeros(13,1)];

wheel = interp1(wheelData.t, wheelData.total, t);
isRunning = double(wheel > runningThreshold);
runPerTrial = ssLocal.getTracesPerStimulus(isRunning', stimMatrix, [0 0]);
running = squeeze(mean(runPerTrial,4)); % [stimuli x repetition]
running = double(running >= 0.4);

gratings.tuningCurveViewer(st(ismember(clu, inclClu)), clu(ismember(clu, inclClu)), ...
    cids, cluYpos, stimMatrix, t, directions, blanks, laserOn, running);

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
