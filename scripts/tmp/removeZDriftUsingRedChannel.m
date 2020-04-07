iPl = 1;
planes = [2 3 4 5];
iPlane = planes(iPl);
nFrames1_2= [3622 3622 NaN 3622]+[6009 6008 NaN 6008];

channels = {'green','red'};
detrended = cell(1,2);
diffs = cell(1,2);
avgTraces = cell(1,2);
for iCh = 1:2
    fid = fopen(fullfile('J:\DATA\', ...
        sprintf('SS041_2015-04-23_3_plane%d_%s.bin',iPlane,channels{iCh})), 'r');
    fread(fid,  512*512*nFrames1_2(iPl), '*int16');
    data = fread(fid,  512*512*(7779), '*int16');
    fclose(fid);
    data = reshape(single(data),512,512,[]);
    
    numPix = 8;
    newSize = 512/numPix;
    dataSmall = mean(reshape(data, [512, numPix, newSize, size(data,3)]),2);
    dataSmall = mean(reshape(dataSmall, numPix,newSize,newSize,size(data,3)),1);
    dataSmall = squeeze(dataSmall);
    
    dataDetrend = detrend(reshape(dataSmall,newSize^2,[])')';
    dataDetrend = reshape(dataDetrend,newSize,newSize,[]);
    detrended{iCh} = dataDetrend;
    
    timeAvg = median(dataDetrend,3);
    allDiffs = bsxfun(@minus,dataDetrend,timeAvg);
    meanDiffs = mean(reshape(allDiffs,[],size(allDiffs,3)),1);
    diffs{iCh} = meanDiffs;
    rho = corr(meanDiffs', reshape(dataDetrend,[],size(dataDetrend,3))');
    indNeg = rho<0;
    diffPixAvg = reshape(dataDetrend,[],size(dataDetrend,3));
    diffPixAvg(indNeg,:) = -diffPixAvg(indNeg,:);
    diffPixAvg = mean(diffPixAvg,1);
    avgTraces{iCh} = diffPixAvg;
end
clear data dataSmall dataDetrend

d = load(fullfile('C:\DATA\InfoStructs\M150305_SS041\2015-04-23\3\', ...
    sprintf('2015-04-23_3_M150305_SS041_2P_plane%03d_ROI.mat',iPlane)));
meta = d.meta;
[~,stimSeq,stimMatrix,frameTimes] = ssLocal.getStimulusResponseInfo(meta);
[~,blankStim] = gratings.getOrientations(stimSeq);

d = load('C:\RESULTS\nonVisualEffects\modelGratingResp\running\kernelResponses_crossval.mat');
results = d.results;

smoothing = 3; %in sec
filtPoly = 3;
ballData = nonVis.getRunningSpeed(meta);
if ~isempty(ballData)
    filtWindow = ceil(smoothing / median(diff(ballData.t)));
    if mod(filtWindow,2) == 0
        filtWindow = filtWindow-1;
    end
    nonVisData = sgolayfilt(ballData.total, filtPoly, filtWindow);
    nonVisTime = ballData.t;
end
runInt = interp1(nonVisTime,nonVisData,frameTimes);

figure('Position',[1 41 1920 1083])
subplot(3,1,1)
plot(nonVisTime,nonVisData,'k')
axis tight
title('Running')
ylabel('Running speed')
ax1 = gca;
subplot(3,1,2)
plot(frameTimes,avgTraces{2},'r')
axis tight
ylabel('Flourescence change in red channel')
title('Red channel')
ax2 = gca;
subplot(3,1,3)
plot(frameTimes,avgTraces{1},'g')
axis tight
ylabel('Flourescence change in green channel')
title('Green channel')
ax3 = gca;
linkaxes([ax1 ax2 ax3],'x')
xlim(frameTimes([1 end]))

figure('Position',[1 41 1920 1083])
subplot(3,1,1)
plot(runInt,'k')
axis tight
title('Running')
ylabel('Running speed')
ax1 = gca;
subplot(3,1,2)
plot(avgTraces{2},'r')
axis tight
ylabel('Flourescence change in red channel')
title('Red channel')
ax2 = gca;
subplot(3,1,3)
plot(avgTraces{1},'g')
axis tight
ylabel('Flourescence change in green channel')
title('Green channel')
ax3 = gca;
linkaxes([ax1 ax2 ax3],'x')
xlim([1 length(runInt)])

figure
plot(avgTraces{2},avgTraces{1},'k.')
xlabel('Red trace')
ylabel('Green trace')

figure('Position',[1 41 1920 1083])
subplot(4,1,1)
plot(nonVisTime,nonVisData,'k')
axis tight
title('Running')
ylabel('Running speed')
ax1 = gca;
subplot(4,1,2)
plot(frameTimes,avgTraces{2},'k')
axis tight
ylabel('Flourescence change in red channel')
title('Red channel')
ax2 = gca;

for iCell = 1:length(results(2).plane(iPl).cellIDs)
    k = results(2).plane(iPl).cellIDs(iCell);
    pars = results(2).plane(iPl).modelPars(iCell);
    kernel = pars.kernel;
    if isempty(kernel)
        prediction = [];
    else
        prediction = models.getStimulusResponsePrediction(stimMatrix, ...
            blankStim, pars.alphaEachTrial, kernel);
        prediction = prediction * max(1,mean(meta.F0(:,k)));
    end   
    baseline = pars.baseline;
    baseline = baseline * max(1,mean(meta.F0(:,k))) + meta.F0(:,k); 
    
    subplot(4,1,3)
    plot(frameTimes,meta.Fcorr(:,k), 'k')
    hold on
%     plot(frameTimes,meta.F0(:,k), 'b', 'LineWidth', 2)
    plot(frameTimes,baseline, 'b', 'LineWidth', 2)
    labels = {'F (neuropil corr.)','Fitted baseline'};
    if ~isempty(prediction)
        plot(frameTimes,baseline+prediction, 'g')
        labels = [labels, {'Prediction'}];
    end
    linkaxes([ax1 ax2 gca],'x')
    xlabel('Time (in s)')
    ylabel('Green fluorescence')
    legend(labels)
    axis tight
    title(sprintf('Cell %d',k))
    xlim(frameTimes([1 end]))
    hold off
    
    subplot(4,4,15)
    plot(avgTraces{2}, baseline, 'k.')
    xlabel('Red channel')
    ylabel('Fitted baseline')
    axis tight square
    subplot(4,4,16)
    plot(avgTraces{2}, meta.Fcorr(:,k), 'k.')
    xlabel('Red channel')
    ylabel('F (neuropil corr.)')
    axis tight square
    pause
end

