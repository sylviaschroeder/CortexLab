%% Pixel average: tuning and trace + running
% Example 1 (drifting gratings)
fid = fopen('J:\DATA\registered\M160706_SS066\2016-08-23\plane2.bin', 'r');
data = fread(fid,  512*512*9844, '*int16');
fclose(fid);
data = reshape(single(data),512*512,[]);

meta = ppbox.infoPopulate('M160706_SS066','2016-08-23',1);
meta.planeFrames = 2 + (0:9843)' .* meta.nPlanes;
frameTimes = ppbox.getFrameTimes(meta);

datAvg = mean(data,1)';
[F_filt,F0] = preproc.removeSlowDrift(datAvg,frameTimes,180,20);

[~, stimSequence, stimMatrix, ~, samplingRate] = ...
    ssLocal.getStimulusResponseInfo(meta);
gratings.plotTuningOfBrainRegions(permute(datAvg,[2 3 1]), samplingRate, stimMatrix, ...
    stimSequence);

smoothing = 3; %in sec
filtPoly = 3;
ballData = nonVis.getRunningSpeed(meta);
filtWindow = ceil(smoothing / median(diff(ballData.t)));
if mod(filtWindow,2) == 0
    filtWindow = filtWindow-1;
end
nonVisData = sgolayfilt(ballData.total, filtPoly, filtWindow);
nonVisTime = ballData.t;

% Example 2 (flickering checkerboard)
fid = fopen('J:\DATA\registered\M160706_SS066\2016-09-08\plane9.bin', 'r');
data = fread(fid,  128*128*4393, '*int16');
fclose(fid);
data = reshape(single(data),128*128,[]);

meta = ppbox.infoPopulate('M160706_SS066','2016-09-08',2);
meta.planeFrames = 9 + (0:4393)' .* meta.nPlanes;
frameTimes = ppbox.getFrameTimes(meta);
frameTimes(end) = [];

datAvg = mean(data,1)';
[F_filt,F0] = preproc.removeSlowDrift(datAvg,frameTimes,180,20);

[~, stimSequence, stimMatrix, ~, samplingRate] = ...
    ssLocal.getStimulusResponseInfo(meta);
gratings.plotTuningOfBrainRegions(permute(datAvg,[2 3 1]), samplingRate, stimMatrix, ...
    stimSequence);

smoothing = 3; %in sec
filtPoly = 3;
ballData = nonVis.getRunningSpeed(meta);
filtWindow = ceil(smoothing / median(diff(ballData.t)));
if mod(filtWindow,2) == 0
    filtWindow = filtWindow-1;
end
nonVisData = sgolayfilt(ballData.total, filtPoly, filtWindow);
nonVisTime = ballData.t;
nonVisInt = interp1(nonVisTime,nonVisData,frameTimes);

%% z-stack imaging
% define data
db_boutons;

rhos = zeros(1,length(db));
pVals = zeros(1,length(db));
for iExp = 1:length(db)
    % get data
    movies = cell(1, length(db(iExp).planes));
    % avgImgs = [];
    nFrames = zeros(1, length(db(iExp).planes));
    for p = 1:length(db(iExp).planes)
        exps = sprintf('%d_',db(iExp).exps);
        exps(end)= [];
        d = load(fullfile('C:\DATA\F\', db(iExp).subject, db(iExp).date, exps, ...
            sprintf('F_%s_%s_plane%d_Nk200.mat', db(iExp).subject, ...
            db(iExp).date, db(iExp).planes(p))));
        ops = d.ops;
        j = find(db(iExp).exp==db(iExp).exps);
        nFrames(p) = ops.Nframes(j);
        fid = fopen(fullfile('J:\DATA\registered', db(iExp).subject, ...
            db(iExp).date, exps, ...
            sprintf('plane%d.bin',db(iExp).planes(p))));
        fread(fid, ops.Ly*ops.Lx*sum(ops.Nframes(1:(j-1))), '*int16');
        data = fread(fid, ops.Ly*ops.Lx*nFrames(p), '*int16');
        data = reshape(data,ops.Ly,ops.Lx,nFrames(p));
        fclose(fid);
        movies{p} = data;
        %     avgImgs = cat(3, avgImgs, mean(data,3));
    end
    numFrames = min(nFrames);
    iPlane = round(length(db(iExp).planes)/2);
    meta = ppbox.infoPopulate(db(iExp).subject,db(iExp).date,db(iExp).exp);
    meta.planeFrames = db(iExp).planes(iPlane) + (0:numFrames-1)' .* meta.nPlanes;
    frameTimes = ppbox.getFrameTimes(meta);
    smoothing = 3; %in sec
    filtPoly = 3;
    ballData = nonVis.getRunningSpeed(meta);
    filtWindow = ceil(smoothing / median(diff(ballData.t)));
    if mod(filtWindow,2) == 0
        filtWindow = filtWindow-1;
    end
    running = sgolayfilt(ballData.total./median(diff(ballData.t))./ 53, ...
        filtPoly, filtWindow);
    running = interp1(ballData.t,running,frameTimes);
    
    % compare grand average of fluorescence (all pixels, all planes) with
    % running
    planeTraces = zeros(numFrames, length(db(iExp).planes));
    for iPlane = 1:length(db(iExp).planes)
        planeTraces(:,iPlane) = squeeze(mean(reshape( ...
            movies{iPlane}(:,:,1:numFrames),[],numFrames),1));
    end
    % remove first frame
    frameTimes(1) = [];
    running(1) = [];
    planeTraces(1,:) = [];
    % remove slow trends
    winSpan = 60; % in s
    fs = 1/median(diff(frameTimes));
    n = ceil(winSpan * fs / 2);
    drifts = zeros(size(planeTraces));
    for k = 1:size(planeTraces,1)
        tmpTraces = planeTraces(max(1,k-n) : min(size(planeTraces,1),k+n),:);
        drifts(k,:) = prctile(tmpTraces, 8);
    end
    planeTraces = planeTraces - drifts;
    totalTrace = mean(planeTraces,2);
    [rho,p] = corr(running',totalTrace);
    rhos(iExp) = rho;
    pVals(iExp) = p;
end

% plot rhos and mark one example experiment
example = 1;
figure
plot(rhos(example),0,'v','MarkerFaceColor','r','MarkerEdgeColor','none')
hold on
rest = setdiff(1:length(rhos),example);
plot(rhos(rest),zeros(1,length(rest)),'v','MarkerFaceColor','k','MarkerEdgeColor','none')
xlim([-.55 .7])
set(gca,'box','off')
xlabel('Corr. coeff.')
title('Corr. of bouton activity with running speed')

% smooth traces
smoothTraces = zeros(size(planeTraces));
for k = 1:size(planeTraces,2)
    smoothTraces(:,k) = smooth(planeTraces(:,k),10);
end
normTraces = bsxfun(@rdivide, bsxfun(@minus, smoothTraces, mean(smoothTraces,1)), ...
    std(smoothTraces,0,1));
% upperLim = repmat(prctile(normTraces,95,1),size(normTraces,1),1);
% normTraces(normTraces > upperLim) = upperLim(normTraces > upperLim);
% lowerLim = repmat(prctile(normTraces,5,1),size(normTraces,1),1);
% normTraces(normTraces < lowerLim) = lowerLim(normTraces < lowerLim);

figure
subplot(6,1,1)
plot(frameTimes, running, 'r')
ylabel('Running speed (cm/s)')
set(gca,'box','off')
axis tight
subplot(6,1,2)
plot(frameTimes, totalTrace, 'k')
ylabel('F (grand average)')
set(gca,'box','off')
axis tight
subplot(6,1,3:6)
imagesc(frameTimes([1 end]), db(iExp).planes([6 end]), normTraces(:,6:end)')
colormap gray
set(gca,'YTick',[6 16],'YTickLabel',[30 70],'box','off')
xlabel('Time (s)')
ylabel('Approximate depth from surface of SC')
colorbar('Ticks',[min(reshape(normTraces(:,6:end),[],1)) ...
    max(reshape(normTraces(:,6:end),[],1))], 'TickLabels', {'low','high'})
% plot(frameTimes, bsxfun(@plus, normTraces, (0:length(db.planes)-1).*(-5)))
figure
plot(running, totalTrace, 'k.','MarkerSize',1)
[rho,p] = corr(running',totalTrace);
title(sprintf('Rho = %.3f (p = %.3f)', rho, p))
axis([-1 22 -25 110])
set(gca,'box','off')
xlabel('Running speed')
ylabel('F (grand average)')

% plot average frame for each plane
for iPlane=1:length(movies)
    figure
    imagesc(mean(movies{iPlane},3))
    colormap gray
    title(sprintf('Plane %d',db.planes(iPlane)))
    axis image off
end



% determine and show aperture to analyse
x = [42 46];
y = [59 63];
gap = 10;
grandImg = zeros(size(avgImgs,1)*4+gap*3, size(avgImgs,2)*4+gap*3);
apertureImg = zeros((y(2)-y(1)+1)*4+gap*3, (x(2)-x(1)+1)*4+gap*3);
frameX = zeros(4,size(avgImgs,3));
frameY = zeros(4,size(avgImgs,3));
for k=1:4
    for j=1:4
        grandImg((k-1)*(size(avgImgs,1)+gap)+1 : k*size(avgImgs,1)+(k-1)*gap, ...
            (j-1)*(size(avgImgs,2)+gap)+1 : j*size(avgImgs,2)+(j-1)*gap) = ...
            avgImgs(:,:,(k-1)*4+j);
        apertureImg((k-1)*(y(2)-y(1)+1+gap)+1 : k*(y(2)-y(1)+1)+(k-1)*gap, ...
            (j-1)*(x(2)-x(1)+1+gap)+1 : j*(x(2)-x(1)+1)+(j-1)*gap) = ...
            avgImgs(y(1):y(2),x(1):x(2),(k-1)*4+j);
        frameY(:,(k-1)*4+j) = y([1 1 2 2]) + (k-1)*(size(avgImgs,1)+gap);
        frameX(:,(k-1)*4+j) = x([1 2 2 1]) + (j-1)*(size(avgImgs,2)+gap);
    end
end
frameY = [frameY; frameY(1,:)];
frameX = [frameX; frameX(1,:)];

figure
imshow(grandImg, [0 max(grandImg(:))])
set(gcf,'Position', [713 42 1204 1074])
hold on
plot(frameX, frameY, 'r')
title('Averages of full frames of all planes')

% figure
% imshow(apertureImg,[0 max(apertureImg(:))])
% set(gcf,'Position', [713 42 1204 1074])
% title('Averages of apertures of all planes')

% average pixels in aperature and make depth-vs-time plot
planes=1:16;
[~,~,nFr] = cellfun(@size,movies);
nFr = min(nFr);
depthPlot = zeros(length(planes), nFr);
for p=1:length(planes)
    ap = movies{planes(p)}(y(1):y(2),x(1):x(2),1:nFr);
    depthPlot(p,:) = mean(reshape(ap,[],nFr),1);
end
% dp = depthPlot-min(depthPlot(:));
% dp = dp ./ max(dp(:)) .* 254 + 1;
dp = depthPlot;
figure
imagesc(dp)
colormap hot
xlabel('Frame')
ylabel('Plane')
title('Depth plot of all experiments')

experiments = 2:3;
exp = 2;
indExp = find(experiments == exp);
p = round(nPlanes/2);
meta = ppbox.infoPopulate('M160706_SS066','2016-09-08',exp);
meta.planeFrames = p + (0:nFrames{p}(indExp)-1)' .* meta.nPlanes;
frameTimes = ppbox.getFrameTimes(meta);
smoothing = 3; %in sec
filtPoly = 3;
ballData = nonVis.getRunningSpeed(meta);
filtWindow = ceil(smoothing / median(diff(ballData.t)));
if mod(filtWindow,2) == 0
    filtWindow = filtWindow-1;
end
nonVisData = sgolayfilt(ballData.total, filtPoly, filtWindow);
nonVisTime = ballData.t;
runInt=interp1(nonVisTime,nonVisData,frameTimes);

frameInd = sum(nFrames{end}(1:indExp-1)) + (1:length(frameTimes));
frameInd(frameInd>size(dp,2)) = [];
dpDetrend = bsxfun(@plus, detrend(dp(:,frameInd)'), median(dp(:,frameInd),2)');
figure
subplot(2,1,1)
plot(nonVisTime,nonVisData,'k')
ylabel('Running speed')
title(sprintf('Experiment %d',exp))
ax1 = gca;
subplot(2,1,2)
imagesc(frameTimes([1 end]),[1 size(dp,1)], dpDetrend')
colormap hot
xlabel('Time (s)')
ylabel('Plane')
title('Depth plot')
ax2 = gca;
linkaxes([ax1 ax2],'x')
xlim(frameTimes([1 end]))

dpZTrans = bsxfun(@rdivide,bsxfun(@minus,dpDetrend,mean(dpDetrend,1)), ...
    std(dpDetrend,0,1));
figure
subplot(5,1,1)
plot(nonVisTime,nonVisData,'k')
ylabel('Running speed')
title(sprintf('Experiment %d',exp))
ax1 = gca;
subplot(5,1,2:5)
plot(frameTimes(1:size(dpZTrans,1)), bsxfun(@plus,dpZTrans,(0:length(planes)-1)*(-5)),'k')
title('Depth plot')
ax2 = gca;
ylim([-5*length(planes) 5])
linkaxes([ax1 ax2],'x')
xlim(frameTimes([1 end]))

boutonRows = 9:10;
bouton=mean(dp(boutonRows,frameInd));
figure
plot(runInt(1:length(bouton)),bouton,'k.')
title(sprintf('Experiment %d',exp))
xlabel('Running speed')
ylabel('Bouton fluorescence')

figure
subplot(2,1,1)
plot(nonVisTime, nonVisData, 'k')
ylabel('Running speed')
title(sprintf('Experiment %d',exp))
ax1 = gca;
subplot(2,1,2)
plot(frameTimes(1:length(bouton)),bouton, 'k')
xlabel('Time (s)')
title(sprintf('Bouton activity (planes %d-%d)',boutonRows(1), boutonRows(end)))
ax2 = gca;
linkaxes([ax1 ax2],'x')
xlim(frameTimes([1 end]))