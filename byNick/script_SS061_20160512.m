





%% load spike data

Fs = 30000;
% clu = readNPY('spike_clusters.npy');
ss = readNPY('spike_times.npy');
st = double(ss)/Fs;
spikeTemplates = readNPY('spike_templates.npy'); % note: zero-indexed
clu = spikeTemplates;
tempScalingAmps = readNPY('amplitudes.npy');

% [cids, cgs] = readClusterGroupsCSV('cluster_groups.csv');
% 
% noiseClusters = cids(cgs==0);
% 
% st = st(~ismember(clu, noiseClusters));
% spikeTemplates = spikeTemplates(~ismember(clu, noiseClusters));
% tempScalingAmps = tempScalingAmps(~ismember(clu, noiseClusters));
% clu = clu(~ismember(clu, noiseClusters));
% cgs = cgs(~ismember(cids, noiseClusters));
% cids = cids(~ismember(cids, noiseClusters));

%% compute template and spike Amplitudes and Depths
load forPRBimecP3opt3
yc = ycoords(connected); xc = xcoords(connected);
temps = readNPY('templates.npy');

winv = readNPY('whitening_mat_inv.npy');

[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW] = ...
    templatePositionsAmplitudes(temps, winv, yc, spikeTemplates, tempScalingAmps);

%% without manual sorting, cluster Ypos and Amp are same as templates

cids = unique(clu);
cluYpos = templateYpos(cids+1);
cluAmps = tempAmps(cids+1);

%% compute firing rates
[cids, spikeCounts] = countUnique(clu);
recDur = st(end);
FRs = spikeCounts./recDur;



%% plot templates by depth

yy = cluYpos; xx = rand(size(yy))*max(yy);

figure;

subplot(2,2,1);

scatter(xx,yy,cluAmps*4);
title('amplitudes');

subplot(2,2,3);
hist(cluAmps, 25);
xlabel('amplitude')

subplot(2,2,2);

scatter(xx,yy,FRs*3);
title('firing rates');

subplot(2,2,4);
hist(FRs, 0:0.5:20);
xlabel('firing rate');
xlim([0 20]);


%% stimulus info for gratings (1st and 3rd exp.)

fid = fopen('SS061_20160512_g0_t0.imec_sync.bin');
syncDat = fread(fid, Inf, '*int16');
fclose(fid);

syncFs = 3000;
tSync = (0:length(syncDat)-1)/syncFs;

flipsDown = tSync(syncDat(1:end-1)==-1 & syncDat(2:end)==-2);
flipsUp = tSync(syncDat(1:end-1)==-2 & syncDat(2:end)==-1);


stimStarts_1 = flipsDown(diff([0 flipsDown])>0.2 & flipsDown>50 & flipsDown<2550); %  empirically determined transition time
stimStops_1 = stimStarts_1+2;
load('\\ZSERVER\Data\trodes\M160426_SS061\20160512\1\Protocol.mat')
stimIDs_1 = zeros(size(stimStarts_1));
for q = 1:size(Protocol.seqnums,1)
    stimIDs_1(Protocol.seqnums(q,:)) = q;
end

stimStarts_3 = flipsDown(diff([0 flipsDown])>0.2 & flipsDown>3665 & flipsDown<4890); %  empirically determined transition time
stimStops_3 = stimStarts_3+2;
load('\\ZSERVER\Data\trodes\M160426_SS061\20160512\3\Protocol.mat')
stimIDs_3 = zeros(size(stimStarts_3));
for q = 1:size(Protocol.seqnums,1)
    stimIDs_3(Protocol.seqnums(q,:)) = q;
end

%% PSTHs 

ampThresh = 2.5;
yposWin = [2100 3200]; % V1
inclClu = cids(cluAmps>ampThresh & cluYpos>yposWin(1) & cluYpos<yposWin(2));

psthViewer(st(ismember(clu, inclClu)), clu(ismember(clu, inclClu)), ...
    stimStarts_1, [-0.5 3.0], stimIDs_1);

psthViewer(st(ismember(clu, inclClu)), clu(ismember(clu, inclClu)), ...
    stimStarts_3, [-0.5 3.0], stimIDs_3);

%% Plot modulation index due to inactivation
% (1) first experiment (laser at electrode location)
% exclude low amplitude and low firing rate neurons (measure firing rate
% only during first 10 repeats of experiment as electrode moved afterwards)
minRate = 0.5;
maxTrial_1 = 26 * 10;
durTrial = 2;
tempWindow = [0.4 1.7];
laserOn_1 = [ones(13,1); zeros(13,1)];
ampThresh = 2.5;

clu_1 = clu(st>=stimStarts_1(1) & st<=stimStarts_1(maxTrial_1)+durTrial);
dur_1 = stimStarts_1(maxTrial_1)+2 - stimStarts_1(1);
[cells_1, spikeC_1] = countUnique(clu_1);
cells_1(spikeC_1/dur_1 < minRate) = [];
firingRateInds_1 = ismember(cids, cells_1);
goodNeurons_1 = cids(firingRateInds_1 & cluAmps>ampThresh);

modIndices_1 = ephys.getInactivationIndices(st(ismember(clu, goodNeurons_1)), ...
    clu(ismember(clu, goodNeurons_1)), stimStarts_1(1:maxTrial_1), tempWindow, ...
    stimIDs_1(1:maxTrial_1), laserOn_1);
modIndices_1(modIndices_1<-10) = -10;
modIndices_1(modIndices_1>10) = 10;
figure, hist(modIndices_1)
xlabel('Modulation index')
ylabel('# Units')
colMap = jet(201);
valid = ~isnan(modIndices_1);
cols = zeros(length(modIndices_1),3);
cols(valid,:) = colMap(round(modIndices_1(valid)*10)+101,:);
% areas = round(abs(modIndices_1)*10);
figure
scatter(xx(ismember(cids, goodNeurons_1(valid))), ...
    yy(ismember(cids, goodNeurons_1(valid))), 50, cols(valid,:))
set(gca, 'Color', [0.5 0.5 0.5])
colormap jet
colorbar('Ticks',0:0.1:1, 'TickLabels', -10:2:10)
ylim([0 max(yy)])
ylabel('Depth in microns')
title('Modulation of V1 and below due to inactivation of V1 (same spot)')

% (2) laser 1 mm anterior of electrode position (entry site)
minTrial_3 = 26*11+11;
dur_3 = stimStops_1(end) - stimStarts_1(minTrial_3) + ...
    (stimStops_3(end) - stimStarts_3(1));
laserOn_3 = [NaN(13,1); zeros(13,1); ones(13,1)];
clu_3 = clu((st>=stimStarts_1(minTrial_3) & st <=stimStops_1(end)) | ...
    (st>=stimStarts_3(1) & st<=stimStops_3(end)));
[cells_3, spikeC_3] = countUnique(clu_3);
cells_3(spikeC_3/dur_3 < minRate) = [];
firingRateInds_3 = ismember(cids, cells_3);
goodNeurons_3 = cids(firingRateInds_3 & cluAmps>ampThresh);

modIndices_3 = ephys.getInactivationIndices(st(ismember(clu, goodNeurons_3)), ...
    clu(ismember(clu, goodNeurons_3)), [stimStarts_1(minTrial_3:end), ...
    stimStarts_3], tempWindow, [stimIDs_1(minTrial_3:end), ...
    stimIDs_3+max(stimIDs_1)], laserOn_3);
modIndices_3(modIndices_3<-10) = -10;
modIndices_3(modIndices_3>10) = 10;
figure, hist(modIndices_3)
xlabel('Modulation index')
ylabel('# Units')
valid_3 = ~isnan(modIndices_3);
cols_3 = zeros(length(modIndices_3),3);
cols_3(valid_3,:) = colMap(round(modIndices_3(valid_3)*10)+101,:);
% areas = round(abs(modIndices_3)*10);
figure
scatter(xx(ismember(cids, goodNeurons_3(valid_3))), ...
    yy(ismember(cids, goodNeurons_3(valid_3))), 50, cols_3(valid_3,:))
set(gca, 'Color', [0.5 0.5 0.5])
colormap jet
colorbar('Ticks',0:0.1:1, 'TickLabels', -10:2:10)
ylim([0 max(yy)])
ylabel('Depth in microns')
title('Modulation of V1 and below due to inactivation of V1 (1mm apart)')

%% load stimulus info for sparse noise

load('\\ZSERVER\Data\trodes\M160426_SS061\20160512\4\Protocol.mat')

expTimeWindow = [4929 5863];

fid = fopen('SS061_20160512_g0_t0.imec_sync.bin');
syncDat = fread(fid, Inf, '*int16');
fclose(fid);

syncFs = 3000;
tSync = (0:length(syncDat)-1)/syncFs;

flipsDown = tSync(syncDat(1:end-1)==-1 & syncDat(2:end)==-2);
flipsUp = tSync(syncDat(1:end-1)==-2 & syncDat(2:end)==-1);

flipsDown = flipsDown(flipsDown>expTimeWindow(1) & flipsDown<expTimeWindow(2)); % empirical
flipsUp = flipsUp(flipsUp>expTimeWindow(1) & flipsUp<expTimeWindow(2));

photodiodeFlips = sort([flipsUp flipsDown]);

%% compute sparse noise frames

load(dat.expFilePath('M160426_SS061', '2016-05-12', 4, 'hw-info', 'master'));
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

% compute firing rates just for this experiment
cids = unique(clu);
[~, spikeCounts] = countUnique([clu(st>expTimeWindow(1) & st<expTimeWindow(2)); cids]);

recDur = diff(expTimeWindow);
FRsThisRec = spikeCounts./recDur;

ampThresh = 2.5;
yposWin = [2100 3200]; % V1
useCIDs = cids(cluAmps>ampThresh & cluYpos>yposWin(1) & cluYpos<yposWin(2) & FRsThisRec>0.2);
whos useCIDs


% useCIDs = [445];

for c = 1:length(useCIDs)
    thisCID = useCIDs(c);
    theseST = st(clu==thisCID);
    
    thisRF = zeros(nX,nY);
    
    for x = 1:nX
        for y = 1:nY            
            [psth, bins, rasterX, rasterY, spikeCounts] = psthRasterAndCounts(theseST, stimEventTimes{x,y,stimEventsToUse}(:), spikeCountingWindow, 0.001);
            thisRF(x,y) = mean(spikeCounts);
        end
    end
    
    figure;
    subplot(2,1,1)
    imagesc(thisRF)
    title(sprintf('cluster %d', thisCID));
    
    subplot(2,1,2);
    imagesc(conv2(thisRF, ones(3,3)/9, 'same'));
    
    pause
end
