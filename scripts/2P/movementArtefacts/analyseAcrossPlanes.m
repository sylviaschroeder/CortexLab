dataset = 10;
planes = 2:4; % 2:5;
exp = 1;
folderMeta = 'C:\DATA\InfoStructs\M150611_SS048\2015-12-02\1'; % 'C:\DATA\InfoStructs\M150305_SS041\2015-04-23\3';
fileMeta = '2015-12-02_1_M150611_SS048_2P_plane%03d_ROI.mat';
folderMetaRed = 'C:\DATA\InfoStructs\red\M150611_SS048\2015-12-02\1'; % 'C:\DATA\InfoStructs\red\M150305_SS041\2015-04-23\3'
folderSVD = 'C:\DATA\F\M150611_SS048\2015-12-02\1_2_3_4_6';
folderSVDRed = 'C:\DATA\F\M150611_SS048\2015-12-02\1_2_3_4_6\red';
fileSVD = 'SVD_M150611_SS048_2015-12-02_plane%d.mat';
prctileNeuron = 8;
prctileWindow = 180; % in sec
prctileWindowRed = 180; % in sec
redFilter = 10; % in sec (median filter)

smoothing = 3; % sec
filtPoly = 3;

pixels = 50;

imRed = cell(1,max(planes));
imGreen = cell(1,max(planes));
mergedRed = cell(1,max(planes));
mergedGreen = cell(1,max(planes));
for p = 2  %planes
    data = load(fullfile(folderSVDRed, sprintf(fileSVD, p)));
    U_r = data.U;
    x = size(U_r,2);
    y = size(U_r,1);
    U_r = double(reshape(U_r,x*y,[]));
    V_r = data.Vcell;
    V_r = double(V_r{exp});
    r = U_r * sum(V_r,2);
    r = (r - min(r(:))) ./ (max(r(:)) - min(r(:)));
    imRed{p} = reshape(r,y,x);
    
    data = load(fullfile(folderSVD, sprintf(fileSVD, p)));
    U_g = data.U;
    x = size(U_g,2);
    y = size(U_g,1);
    U_g = double(reshape(U_g,x*y,[]));
    V_g = data.Vcell;
    V_g = double(V_g{exp});
    g = U_g * sum(V_g,2);
    g = (g - min(g(:))) ./ (max(g(:)) - min(g(:)));
    imGreen{p} = reshape(g,y,x);
    
    data = load(fullfile(folderMeta, sprintf(fileMeta,p)));
    m = data.meta;
    [mask,cols] = spatial.plotCellMasks(m,0,1);
    mergedRed{p} = spatial.mergeImageWithMask(imRed{p},mask,cols);
    mergedGreen{p} = spatial.mergeImageWithMask(imGreen{p},mask,cols);
    
    data = load(fullfile(folderMetaRed, sprintf(fileMeta,p)));
    m.R = data.meta.F;
    m.Rcorr = data.meta.Fcorr;
    m.R0 = data.meta.F0;
    meta(p) = m;
end
% data = load('C:\RESULTS\Preprocessing\globalRedTrace_SS041_2015-04-23_exp3_plane2.mat');
% globalRed = data.globalRed;

for p = planes
    figure('position',[5 42 3828 1074])
    subplot(1,4,1)
    imagesc(imRed{p})
    colormap gray
    axis image
    title(sprintf('Plane %d: Red channel',p))
    subplot(1,4,2)
    imshow(mergedRed{p})
    axis image
    subplot(1,4,3)
    imagesc(imGreen{p})
    colormap gray
    axis image
    title('Green channel')
    subplot(1,4,4)
    imshow(mergedGreen{p})
    axis image
end

frameTimes = ppbox.getFrameTimes(meta(planes(1)));
samplingRate = 1/median(diff(frameTimes));
green = cell(1,max(planes));
greenPercentiles = cell(1,max(planes));
deltaF = cell(1,max(planes));
red = cell(1,max(planes));
redDetrended = cell(1,max(planes));
redPercentiles = cell(1,max(planes));
for p = 2 %planes
    t = ppbox.getFrameTimes(meta(p));
    green{p} = interp1(t,meta(p).Fcorr,frameTimes, 'pchip');
    [~, prctlTraces] = preproc.removeSlowDrift(meta(p).Fcorr, ...
        samplingRate, prctileWindow, prctileNeuron);
    [~, prctlTracesRed] = preproc.removeSlowDrift(meta(p).R, ...
        samplingRate, prctileWindowRed, prctileNeuron);
    greenPercentiles{p} = interp1(t, prctlTraces, frameTimes, 'pchip');
    redPercentiles{p} = interp1(t, prctlTracesRed, frameTimes, 'pchip');
    tr = bsxfun(@rdivide, meta(p).Fcorr-prctlTraces, bsxfun(@max, 1, ...
        mean(prctlTraces,1)));
    deltaF{p} = interp1(t, tr, frameTimes, 'pchip');
    red{p} = interp1(t,meta(p).R,frameTimes,'pchip');
    redDetrended{p} = interp1(t,meta(p).R-prctlTracesRed,frameTimes,'pchip');
end

ballData = nonVis.getRunningSpeed(meta(2));
filtWindow = ceil(smoothing / median(diff(ballData.t)));
if mod(filtWindow,2) == 0
    filtWindow = filtWindow-1;
end
total = ballData.total ./ median(diff(ballData.t)) ./ 53;
running = sgolayfilt(total, filtPoly, filtWindow);
running = interp1(ballData.t, running, frameTimes, 'pchip')';

[~, ~, stimMatrix] = ssLocal.getStimulusResponseInfo(meta(planes(1)));

data = load('C:\RESULTS\Preprocessing\data\globalRedTraces_SS048_2015-12-02_1.mat');
globalRed = data.globalRed;
globalRedPerStim = ssLocal.getTracesPerStimulus(globalRed.percentiled, ...
    stimMatrix,[0 0]);
globalRedPerStim = squeeze(mean(globalRedPerStim,4)); % [stim x trial]
data = load('C:\RESULTS\nonVisualEffects\modelGratingResp\kernelFit\results.mat');
results = data.results;

nRF = round(samplingRate * redFilter);
orange = lines(2);
orange = orange(2,:);
for p = 2 %planes
    deltaFPerStim = ssLocal.getTracesPerStimulus(deltaF{p},stimMatrix,[0 0]);
    deltaFPerStim = mean(deltaFPerStim,4); % [neuron x stim x trial]
    redPerStim = ssLocal.getTracesPerStimulus(redDetrended{p},stimMatrix,[0 0]);
    redPerStim = mean(redPerStim,4);
    for c = 1:size(meta(p).F,2)
        if meta(p).ROI.isDuplicate(c)==1 %|| isempty(meta(p).ROI.duplicates(c).plane)
            continue
        end
        
        % show red with moving percentiles + 
        figure('position',[1 41 1920 1083])
        subplot(2,1,1)
        hold on
        plot(frameTimes, green{p}(:,c), 'Color', [0 .5 0])
        plot(frameTimes, greenPercentiles{p}(:,c), 'Color', [0 .5 .5])
        subplot(2,1,2)
        hold on
        plot(frameTimes, red{p}(:,c), 'Color', [.5 0 0])
        plot(frameTimes, redPercentiles{p}(:,c), 'Color', [.5 0 .5])
        
        pause
        close
    end
        
%         % figure (1): Traces of running, global red, local red (detrended),
%         % delta F/F, and raw F
%         ax = zeros(1,7);
%         f1 = figure('position',[1 41 1920 1083]);
%         
%         subplot(12,1,1)
%         plot(frameTimes, running, 'k')
%         ax(1) = gca;
%         [rho,pVal] = corr(running, deltaF{p}(:,c));
%         title(sprintf('Cell %d of plane %d (rho=%.2f, p=%.4f)', ...
%             c, p, rho, pVal))
%         ylabel('Running')
%         axis tight
%         set(gca,'box','off')
%         
%         subplot(12,1,2)
%         plot(frameTimes, globalRed.percentiled, 'r')
%         ax(2) = gca;
%         ylabel('Global red')
%         axis tight
%         set(gca,'box','off')
%         
%         trR = NaN(length(frameTimes), max(planes));
%         trDF = NaN(length(frameTimes), max(planes));
%         trF = NaN(length(frameTimes), max(planes));
%         trR(:,p) = redDetrended{p}(:,c);
%         trDF(:,p) = deltaF{p}(:,c);
%         trF(:,p) = green{p}(:,c);
% %         for d = 1:length(meta(p).ROI.duplicates(c).plane)
% %             plane = meta(p).ROI.duplicates(c).plane(d);
% %             id = meta(p).ROI.duplicates(c).ID(d);
% %             trR(:,plane) = tracesRed{plane}(:,id);
% %             trDF(:,plane) = deltaF{plane}(:,id);
% %             trF(:,plane) = traces{plane}(:,id);
% %         end
%         filteredRed = medfilt1(trR(:,p),nRF);
%         mdl = fitlm(filteredRed, trDF(:,p), 'RobustOpts', 'on');
%         slope = mdl.Coefficients.Estimate(2);
%         correctedGreen = trDF(:,p) - slope.*filteredRed;
%         
%         subplot(12,1,3:4)
%         hold on
%         plot(frameTimes,red{p}(:,c))
%         plot(frameTimes,redPercentiles{p}(:,c))
%         plot(frameTimes,imopenedRed{p}(:,c))
%         legend('raw','prctls','minmax')
%         ax(7) = gca;
%         ylabel('Raw red')
%         axis tight
%         set(gca,'box','off')
%         
%         subplot(12,1,5:6)
%         hold on
% %         tr = bsxfun(@rdivide, bsxfun(@minus, tr, nanmean(tr,1)), nanstd(tr,1));
%         offset = 0;
%         pl = [];
%         depths = [];
%         for k = 1:size(trR,2)
%             if all(isnan(trR(:,k)))
%                 continue
%             end
%             col = [.5 0 0];
%             if k ~= p
%                 col = .5*col + .5;
%             end
%             maxi = max(trR(:,k));
%             mini = min(trR(:,k));
%             plot(frameTimes,trR(:,k)+offset-maxi, 'Color', col)
%             if k == p
%                 plot(frameTimes,filteredRed+offset-maxi, ...
%                     'Color', orange, 'LineWidth', 2)
%             end
%             depths(end+1) = offset - .5*(maxi-mini);
%             pl(end+1) = k;
%             offset = offset - (maxi-mini);
%         end
%         depths = flip(depths);
%         pl = flip(pl);
%         set(gca, 'YTick', depths, 'YTickLabel', pl,'box','off')
%         ylabel('Detrended red')
%         ax(3) = gca;
%         axis tight
%         
%         subplot(12,1,7:8)
%         plot(frameTimes,correctedGreen, 'Color', orange)
%         set(gca,'box','off')
%         ylabel('Corrected \DeltaF/F')
%         ax(4) = gca;
%         axis tight
%         
%         subplot(12,1,9:10)
%         hold on
% %         tr = bsxfun(@rdivide, bsxfun(@minus, tr, nanmean(tr,1)), nanstd(tr,1));
%         offset = 0;
%         pl = [];
%         depths = [];
%         for k = 1:size(trDF,2)
%             if all(isnan(trDF(:,k)))
%                 continue
%             end
%             col = [0 .5 0];
%             if k ~= p
%                 col = .5*col + .5;
%             end
%             maxi = max(trDF(:,k));
%             mini = min(trDF(:,k));
%             plot(frameTimes,trDF(:,k)+offset-maxi, 'Color', col)
%             depths(end+1) = offset - .5*(maxi-mini);
%             pl(end+1) = k;
%             offset = offset - (maxi-mini);
%         end
%         depths = flip(depths);
%         pl = flip(pl);
%         set(gca, 'YTick', depths, 'YTickLabel', pl,'box','off')
%         ylabel('\DeltaF/F')
%         ax(5) = gca;
%         axis tight
%         
%         subplot(12,1,11:12)
%         hold on
% %         tr = bsxfun(@rdivide, bsxfun(@minus, tr, nanmean(tr,1)), nanstd(tr,1));
%         offset = 0;
%         pl = [];
%         depths = [];
%         for k = 1:size(trF,2)
%             if all(isnan(trF(:,k)))
%                 continue
%             end
%             col = [0 0 .5];
%             if k ~= p
%                 col = .5*col + .5;
%             end
%             maxi = max(trF(:,k));
%             mini = min(trF(:,k));
%             plot(frameTimes,trF(:,k)+offset-maxi, 'Color', col)
%             plot(frameTimes,greenPercentiles{p}(:,c)+offset-maxi)
%             plot(frameTimes,imopenedGreen{p}(:,c)+offset-maxi)
%             depths(end+1) = offset - .5*(maxi-mini);
%             pl(end+1) = k;
%             offset = offset - (maxi-mini);
%         end
%         depths = flip(depths);
%         pl = flip(pl);
%         set(gca, 'YTick', depths, 'YTickLabel', pl,'box','off')
%         ylabel('F')
%         xlabel('Time (s)')
%         ax(6) = gca;
%         axis tight
%         
%         linkaxes(ax,'x')
%         xlim(frameTimes([1 end]))
%         
%         % Figure (3): ROI within mean green and red image
%         f3 = figure('Position',[3145 335 695 675]);
%         im = zeros(size(imRed{p}));
%         im(meta(p).ROI.CellMaps{c}) = 1;
%         im2 = zeros(size(imRed{p},1)+2*pixels,size(imRed{p},2)+2*pixels);
%         im2(pixels+(1:size(imRed{p},1)),pixels+(1:size(imRed{p},2))) = im;
%         map = find(im2==1);
%         stats = regionprops(im2, 'Centroid');
%         im3 = zeros(size(imRed{p},1)+2*pixels,size(imRed{p},2)+2*pixels);
%         im3(pixels+(1:size(imRed{p},1)),pixels+(1:size(imRed{p},2))) = ...
%             (imRed{p}-min(imRed{p}(:)))./(max(imRed{p}(:))-min(imRed{p}(:)));
%         img = repmat(im3,1,1,3);
%         img = reshape(img,numel(im3),3);
%         img(:,1) = img(:,1) + 0.3 * im2(:);
%         img(map,:) = img(map,:)./1.3;
%         img = reshape(img,[size(im3) 3]);
%         xPixels = round(stats.Centroid(1))+(-pixels:pixels);
%         yPixels = round(stats.Centroid(2))+(-pixels:pixels);
%         subplot(2,2,1)
%         imshow(img(yPixels,xPixels,:))
%         title('Red channel')
%         subplot(2,2,2)
%         imagesc(im3(yPixels,xPixels))
%         colormap gray
%         im3 = zeros(size(imRed{p},1)+2*pixels,size(imRed{p},2)+2*pixels);
%         im3(pixels+(1:size(imRed{p},1)),pixels+(1:size(imRed{p},2))) = ...
%             (imGreen{p}-min(imGreen{p}(:)))./(max(imGreen{p}(:))-min(imGreen{p}(:)));
%         img = repmat(im3,1,1,3);
%         img = reshape(img,numel(im3),3);
%         img(:,1) = img(:,1) + 0.3 * im2(:);
%         img(map,:) = img(map,:)./1.3;
%         img = reshape(img,[size(im3) 3]);
%         subplot(2,2,3)
%         imshow(img(yPixels,xPixels,:))
%         title('Green channel')
%         subplot(2,2,4)
%         imagesc(im3(yPixels,xPixels))
%         colormap gray
%         
%         % Figure (2): Scatter plots
%         f2 = figure('position',[1923 335 1205 675]);
%         
%         subplot(2,2,1)
%         plot(running, trDF(:,p), 'k.', 'MarkerSize',3)
%         rho = corr(running, trDF(:,p));
%         xlabel('Running')
%         ylabel('Green')
%         title(sprintf('Rho = %.3f', rho))
%         subplot(2,2,2)
%         plot(running, filteredRed, 'k.', 'MarkerSize',3)
%         rho = corr(running, filteredRed);
%         xlabel('Running')
%         ylabel('Filtered red')
%         title(sprintf('Rho = %.3f', rho))
%         
%         subplot(2,2,3)
%         plot(filteredRed, trDF(:,p), 'k.', 'MarkerSize',3)
%         hold on
%         xLim = get(gca,'XLim');
%         plot(xLim, predict(mdl, xLim'), 'Color', orange,'LineWidth',2)
%         rho = corr(filteredRed, trDF(:,p));
%         xlabel('Filtered red')
%         ylabel('Green')
%         title(sprintf('Rho = %.3f', rho))
%         subplot(2,2,4)
%         plot(running, correctedGreen, 'k.', 'MarkerSize',3)
%         rho = corr(running, correctedGreen);
%         xlabel('Running')
%         ylabel('Corrected \DeltaF/F')
%         title(sprintf('Rho = %.3f', rho))
%         
%         pause
%         close([f1 f2 f3])
%     end
end

%% Use globel red signal
planes = 2:4;
exp = 1; % 3;
folderSuite2P = 'C:\DATA\F\M150611_SS048\2015-12-02\1_2_3_4_6'; % 'C:\DATA\F\M150305_SS041\2015-04-23\1_2_3_4';
file2P = 'F_M150611_SS048_2015-12-02_plane%d_Nk650.mat'; % 'F_M150305_SS041_2015-04-23_plane%d_Nk650.mat';
folderBinFiles = 'J:\DATA';
binFile = 'SS048_2015-12-02_1-4_6_plane%d_red.bin'; % 'SS041_2015-04-23_3_plane%d_red.bin';
folderMeta = 'C:\DATA\InfoStructs\M150611_SS048\2015-12-02\1'; % 'C:\DATA\InfoStructs\M150305_SS041\2015-04-23\3';
fileMeta = '2015-12-02_1_M150611_SS048_2P_plane%03d_ROI.mat';

nFrames = zeros(max(planes),exp);
planeData = cell(1,max(planes));
numPix = 8; % average over this number of pixels
% filterCutoff = 1/240; % oscillations slower than 1 minute
% filterOrder = 3;
prctileNeuron = 8;
prctileWindow = 180; % in sec

for iPlane = planes
    data = load(fullfile(folderSuite2P, sprintf(file2P,iPlane)));
    nFrames(iPlane,1:exp) = data.ops.Nframes(1:exp);
    Lx = data.ops.Lx;
    Ly = data.ops.Ly;
    
    fid = fopen(fullfile(folderBinFiles, sprintf(binFile,iPlane)));
    fread(fid,Ly*Lx*sum(nFrames(iPlane,1:end-1)), '*int16');
    data = fread(fid,Ly*Lx*nFrames(iPlane,end),'*int16');
    fclose(fid);
    data = reshape(single(data),Ly,Lx,[]);
    newSize = Ly/numPix;
    dataSmall = mean(reshape(data, [Ly, numPix, newSize, size(data,3)]),2);
    dataSmall = mean(reshape(dataSmall, numPix,newSize,newSize,size(data,3)),1);
    dataSmall = squeeze(dataSmall);
    dataSmall(:,:,1) = dataSmall(:,:,2);
    planeData{iPlane} = reshape(dataSmall,newSize*newSize,[])';
    clear data dataSmall
end

data = load(fullfile(folderMeta, sprintf(fileMeta,planes(1))));
frameTimes = ppbox.getFrameTimes(data.meta);
samplingRate = 1/median(diff(frameTimes));
ballData = nonVis.getRunningSpeed(data.meta);
filtWindow = ceil(3 / median(diff(ballData.t)));
if mod(filtWindow,2) == 0
    filtWindow = filtWindow-1;
end
total = ballData.total ./ median(diff(ballData.t)) ./ 53;
running = sgolayfilt(total, 3, filtWindow);
running = interp1(ballData.t, running, frameTimes, 'pchip')';

% concatenate all planes
allRed = cat(2,planeData{:});

% get 1st PCA from data where moving 8th prctile is subtracted
n = round(samplingRate * prctileWindow/2);
prctileRed = zeros(size(allRed));
parfor k = 1:size(allRed,1)
    if mod(k,50)==0
        fprintf('%d ',k)
    end
    tmpTraces = allRed(max(1,k-n) : min(size(allRed,1),k+n),:);
    prctileRed(k,:) = prctile(tmpTraces, prctileNeuron);
end
fprintf('\n')
allData = cat(2,planeData{:}) - prctileRed;
allData = bsxfun(@minus, allData, mean(allData,1));
gpuDevice(1);
allData = gpuArray(allData);
[U,S,V] = svd(allData,'econ');
globalRed.percentiled = gather(U(:,1));
rhoRed.percentiled = corr(running,globalRed.percentiled);
if rhoRed.percentiled<0
    globalRed.percentiled = -globalRed.percentiled;
    rhoRed.percentiled = -rhoRed.percentiled;
end

% get 1st PCA from detrended data
allData = detrend(allRed);
allData = bsxfun(@minus, allData, mean(allData,1));
gpuDevice(1);
allData = gpuArray(allData);
[U,S,V] = svd(allData,'econ');
globalRed.detrended = gather(U(:,1));
rhoRed.detrended = corr(running,globalRed.detrended);
if rhoRed.detrended<0
    globalRed.detrended = -globalRed.detrended;
    rhoRed.detrended = -rhoRed.detrended;
end

% get 1st PCA from frame normalised data
allData = bsxfun(@rdivide, allRed, mean(allRed,2));
allData = bsxfun(@minus, allData, mean(allData,1));
gpuDevice(1);
allData = gpuArray(allData);
[U,S,V] = svd(allData,'econ');
globalRed.normalised = gather(U(:,1));
rhoRed.normalised = corr(running,globalRed.normalised);
if rhoRed.normalised<0
    globalRed.normalised = -globalRed.normalised;
    rhoRed.normalised = -rhoRed.normalised;
end

figure
subplot(4,1,1)
plot(frameTimes,running,'k')
ylabel('running speed')
set(gca,'box','off')
axis tight
ax(1)=gca;
cols = lines(3);
subplot(4,1,2)
plot(frameTimes,globalRed.percentiled,'Color',cols(1,:))
title('1st PCA of red channel with moving 8th pcrtl. subtracted')
set(gca,'box','off')
axis tight
ax(2) = gca;
subplot(4,1,3)
plot(frameTimes,globalRed.detrended,'Color',cols(2,:))
title('1st PCA of detrended red channel')
set(gca,'box','off')
axis tight
ax(3) = gca;
subplot(4,1,4)
plot(frameTimes,globalRed.normalised,'Color',cols(3,:))
title('1st PCA of frame normalised red channel')
set(gca,'box','off')
axis tight
ax(4) = gca;

% for comparison, calculate red trace based on each plane separately
% red = NaN(size(allData,1),max(planes));
% for iPlane = planes
%     data = bsxfun(@minus, planeDataFilt{iPlane}, mean(planeDataFilt{iPlane},1));
%     data = gpuArray(data);
%     [U,S,V] = svd(data,'econ');
%     red(:,iPlane) = gather(U(:,1));
%     r = corr(running',red(:,iPlane));
%     if r<0
%         red(:,iPlane) = -red(:,iPlane);
%     end
% end
% plot global trace and all plane traces
% figure, hold on
% plot(frameTimes,globalRed,'r')
% for k=planes,plot(frameTimes,red(:,k)-.12*(k-1),'k'),end
% set(gca,'YTick',(-4:0).*.12,'YTickLabel',{'Plane 5','Plane 4','Plane 3','Plane 2','Global'})
% xlim(frameTimes([1 end]))

green = [];
cellPlanes = [];
cellIDs = [];
for iPlane = planes
    data = load(fullfile(folderMeta, sprintf(fileMeta,iPlane)));
    ids = find(data.meta.ROI.isDuplicate==0);
    t = ppbox.getFrameTimes(data.meta);
%     tr = bsxfun(@rdivide, data.meta.Fcorr(:,ids) - ...
%         data.meta.F0(:,ids), bsxfun(@max, 1, mean(data.meta.F0(:,ids),1)));
%     traces = [traces,interp1(t,tr,frameTimes,'pchip')];
    green = [green,interp1(t,data.meta.Fcorr(:,ids),frameTimes,'pchip')];
    cellIDs = [cellIDs,ids'];
    cellPlanes = [cellPlanes, ones(1,length(ids)) * iPlane];
end
prctlTraces = zeros(size(green));
for k = 1:size(green,1)
    tmpTraces = green(max(1,k-n) : min(size(green,1),k+n),:);
    prctlTraces(k,:) = prctile(tmpTraces, prctileNeuron);
end
% filtered = filtfilt(b,a,prctlTraces);
% tracesFiltered = traces - filtered;
tracesFiltered = bsxfun(@rdivide, green - prctlTraces, ...
    bsxfun(@max,1,mean(prctlTraces,1)));
% tracesFiltered = traces;

[rho,p] = corr(running, tracesFiltered);
[~,order] = sort(rho,'ascend');

for k=1:length(order)
    ax = zeros(3,1);
    figure('Position',[5 42 1915 1074])
    subplot(4,1,1)
    plot(frameTimes,running,'k')
    ylabel('Running speed')
    title(sprintf('Cell %d in plane %d',cellIDs(order(k)),cellPlanes(order(k))))
    ax(1) = gca;
    subplot(4,1,2)
    plot(frameTimes,globalRed.percentiled,'r')
    ylabel('Global red component')
    ax(2) = gca;
    subplot(4,1,3)
    plot(frameTimes,tracesFiltered(:,order(k)),'Color',[0 .5 0])
    ylabel('Green F')
    ax(3) = gca;
    linkaxes(ax,'x')
    xlim(frameTimes([1 end]))
    
    subplot(4,4,13)
    plot(running,tracesFiltered(:,order(k)),'k.','MarkerSize',1)
    xlabel('Running')
    ylabel('Green F')
    title(sprintf('Rho=%.3f',rho(order(k))))
    subplot(4,4,14)
    plot(running,globalRed.percentiled,'k.','MarkerSize',1)
    xlabel('Running')
    ylabel('Red')
    title(sprintf('Rho=%.3f',rhoRed.percentiled))
    subplot(4,4,15)
    plot(globalRed.percentiled,tracesFiltered(:,order(k)),'k.','MarkerSize',1)
    xlabel('Red')
    ylabel('Green F')
    r = corr(globalRed.percentiled,tracesFiltered(:,order(k)));
    title(sprintf('Rho=%.3f',r))
    subplot(4,4,16)
    mdl = fitlm(globalRed.percentiled, tracesFiltered(:,order(k)));
    pred = feval(mdl,globalRed.percentiled);
    traceRest = tracesFiltered(:,order(k))-pred;
    plot(running,traceRest,'k.','MarkerSize',1)
    xlabel('Running')
    ylabel('Green-Red')
    r = corr(running,traceRest);
    title(sprintf('Rho=%.3f',r))
    
    pause
    close gcf
end