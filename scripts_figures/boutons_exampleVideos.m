%% Make movie for darkness
% % SS070, 2016-10-18, plane 9, exp. 3
% movieFile = 'J:\DATA\registered\M160923_SS070\2016-10-18\1_2_3_4\interpolated\M160923_SS070_2016-10-18_1_2_3_4_plane9.bin';
% regFile = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\F\M160923_SS070\2016-10-18\1_2_3_4\regops_M160923_SS070_2016-10-18.mat';
% metaFile_dark = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\InfoStructs\M160923_SS070\2016-10-18\3\2016-10-18_3_M160923_SS070_2P_plane009_ROI.mat';
% outputFile_dark = 'C:\STORAGE\OneDrive - University College London\Presentations\boutonMedia\SS070_2016-10-18_3_darkness-running.avi';
% plane = 9;
% exp_dark = 3;
% SS076, 2017-10-02, plane 6, exp. 1
% movieFile = 'J:\DATA\registered\SS076\2017-10-02\1_2_3\interpolated\SS076_2017-10-02_1_2_3_plane6.bin';
% regFile = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\F\SS076\2017-10-02\1_2_3\regops_SS076_2017-10-02.mat';
% metaFile_dark = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\InfoStructs\SS076\2017-10-02\1\2017-10-02_1_SS076_2P_plane006_ROI.mat';
% outputFile_dark = 'C:\STORAGE\OneDrive - University College London\Presentations\boutonMedia\SS076_2017-10-02_1_darkness-running.avi';
% outputFile_dark_triggered = 'C:\STORAGE\OneDrive - University College London\Presentations\boutonMedia\SS076_2017-10-02_1_darkness-runningTriggered.avi';
% plane = 6;
% exp_dark = 1;
% M170821_SS075, 2017-09-13, plane 5, exp. 1
movieFile = 'J:\Data\registered\M170821_SS075\2017-09-13\1_2_4\interpolated\M170821_SS075_2017-09-13_1_2_4_plane5.bin';
regFile = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\F\M170821_SS075\2017-09-13\1_2_4\regops_M170821_SS075_2017-09-13.mat';
metaFile_dark = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\InfoStructs\M170821_SS075\2017-09-13\1\2017-09-13_1_M170821_SS075_2P_plane005_ROI.mat';
outputFile_dark = 'C:\STORAGE\OneDrive - University College London\Presentations\boutonMedia\M170821_SS075_2017-09-13_1_darkness-running';
outputFile_dark_triggered = 'C:\STORAGE\OneDrive - University College London\Presentations\boutonMedia\M170821_SS075_2017-09-13_1_darkness-runningTriggered';
plane = 5;
exp_dark = 1;


% for running triggered movie
thresh = 1;
movDur = [-10 26];
speedUp_dark_triggered = 2;
minSamplesRunning = 20;

speedUp_dark = 10;
smoothStd = 1;
% load movie
data = load(regFile);
ops = data.ops1{plane};
frames = ops.Nframes;
ops1{1} = ops;
ds = ops.DS_allPlanes;
fr = sum(frames(1:exp_dark-1))+1:sum(frames(1:exp_dark));
if length(size(ops.CorrFrame_allPlanes)) == 2
    invalid = isnan(ops.CorrFrame_allPlanes(:,plane));
    ds(invalid,:,:) = [];
    ds = ds(fr,:,:);
    ops1{1}.CorrFrame_allPlanes(invalid,:) = [];
    ops1{1}.CorrFrame_allPlanes = ops1{1}.CorrFrame_allPlanes(fr,:);
else
    planeID = ops.planesToInterpolate == plane;
    invalid = all(isnan(ops.CorrFrame_allPlanes(:,planeID,planeID)),2);
    ds(invalid,:,:,:,:) = [];
    ds = ds(fr,:,:,:,:);
    ops1{1}.CorrFrame_allPlanes(invalid,:,:) = [];
    ops1{1}.CorrFrame_allPlanes = ops1{1}.CorrFrame_allPlanes(fr,:,:);
end
ops1{1}.DS_allPlanes = ds;
ops1{1}.usedPlanes = ops1{1}.usedPlanes(fr,:);
ops1{1}.badframes = false(1,length(fr));
ops1 = getRangeNoOutliers(ops,ops1);
ops1{1}.yrange = ops.yrange;
ops1{1}.xrange = ops.xrange;
ops = ops1{1};
Ly = ops.Ly;
Lx = ops.Lx;
fid = fopen(movieFile, 'r');
if exp_dark > 1
    fseek(fid, Ly*Lx*sum(frames(1:exp_dark-1))*2,'bof');
end
mov = fread(fid, Ly*Lx*frames(exp_dark), '*int16');
fclose(fid);
mov = reshape(mov, Ly, Lx, []);
mov = mov(ops.yrange, ops.xrange, :);
% get frame times
data = load(metaFile_dark);
meta = data.meta;
frameTimes = ppbox.getFrameTimes(meta);
sr = 1/median(diff(frameTimes));

% try: convolve pixels with gauss window; apply SVD and only take first few
% components
stdSamples = round(smoothStd * sr);
convWindow = normpdf(-4*stdSamples:4*stdSamples, 0, stdSamples);
L = length(convWindow);
mov = cat(3, repmat(mean(mov(:,:,1:floor(L/2)),3),1,1,L), ...
    mov, repmat(mean(mov(:,:,end-floor(L/2):end),3),1,1,L));
for y = 1:size(mov,1)
    for x = 1:size(mov,2)
        mov(y,x,:) = conv(squeeze(mov(y,x,:)), convWindow, 'same');
    end
end
mov(:,:,[1:L, end-L+1:end]) = [];

% get running speed
ballData = nonVis.getRunningSpeed(meta);
stdSamples = round(smoothStd / median(diff(ballData.t)));
convWindow = normpdf(-4*stdSamples:4*stdSamples, 0, stdSamples);
running = conv(ballData.total, convWindow, 'same');
running = running / median(diff(ballData.t)) / 53;
running = interp1(ballData.t, running, frameTimes);

running_binary = running > thresh;
sta = find(diff(running_binary)==1);
sto = find(diff(running_binary)==-1);
if sta(1)>sto(1)
    sto(1) = [];
end
ind = find(sta(2:end) - sto(1:length(sta)-1) < minSamplesRunning);
sta(ind+1) = [];
sto(ind) = [];
ind = (sto - sta(1:length(sto))) < minSamplesRunning;
sta(ind) = [];
if sta(end)>length(running_binary)-minSamplesRunning
    sta(end) = [];
end

movDurFr = round(movDur .* sr);
timeInFrames = movDurFr(1):movDurFr(2);
tim = timeInFrames / sr;
framesTriggered = sta'+1 + timeInFrames;
framesTriggered(framesTriggered<1 | framesTriggered>size(mov,3)) = NaN;
ind = framesTriggered;
ind(isnan(ind)) = 1;
% runTriggered = double(running_binary(ind));
runTriggered = running(ind);
runTriggered(isnan(framesTriggered)) = NaN;
runTriggered = nanmean(runTriggered,1);
runTriggered = (runTriggered-min(runTriggered)) ./ ...
    (max(runTriggered)-min(runTriggered));
avgMov = zeros(size(mov,1),size(mov,2),size(framesTriggered,2));
for k = 1:size(framesTriggered,2)
    ind = framesTriggered(~isnan(framesTriggered(:,k)),k);
%     if timeInFrames(k) < 1
%         ind = ind(running_binary(ind));
%     else
%         ind = ind(~running_binary(ind));
%     end
    avgMov(:,:,k) = mean(mov(:,:,ind),3);
end
avgMov = avgMov - mean(mov(:,:,~running_binary),3);

% make running triggered movie
maxi = max(abs(avgMov(:))) * .6;
running_scaled = round(runTriggered.*length(ops.yrange));
xpos = [-round(length(ops.xrange)/40) 0];
% vid = VideoWriter(outputFile_dark_triggered,'Uncompressed AVI');
vid = VideoWriter(outputFile_dark_triggered,'MPEG-4');
vid.FrameRate = sr * speedUp_dark_triggered;
% vid.Quality = 100;
open(vid);
figure
weights = linspace(0,1,50)';
cm = [[0 0 1].*flip(weights)+[1 1 1].*weights; ...
    [1 1 1].*flip(weights)+[1 0 0].*weights];
cm(length(weights),:) = [];
colormap(cm)
an = annotation('textbox',[.15 .83 .1 .06],'String', ...
    sprintf('%.0f s',0),'EdgeColor','none','FontSize',14, ...
    'FontWeight','bold','HorizontalAlignment','right');
an2 = annotation('textbox',[0.73 0.2 0.12 0.05],'String','20 \mum', ...
    'EdgeColor','none','FontSize',12);
anTitle = annotation('textbox',[0.4 0.83 0.22 0.06],'String','In darkness', ...
    'EdgeColor','none','FontSize',14,'FontWeight','bold');
anScale = annotation('line',[0.7 0.868],[0.25 0.25],'LineWidth',2,'Color','k');
for k = 1:length(running_scaled)
    imagesc(avgMov(:,:,k), [-maxi maxi])
%     colormap gray
    axis image off
    hold on
    fill([xpos flip(xpos)],length(ops.yrange)+.5-[0 0 [1 1].*running_scaled(k)],'r', 'EdgeColor','none','FaceColor','r')
    hold off
    an.String = sprintf('%.0f s',tim(k));
    txt = text(6,140,'running speed');
    txt.Color='r';
    txt.Rotation=90;
    txt.FontSize=12;
    txt.FontWeight='bold';
%     clear frame
    frame = getframe;
%     pause
    writeVideo(vid,frame)
%     writeVideo(vid,frame.cdata)
end
imagesc(zeros(size(avgMov,1),size(avgMov,2)), [-maxi maxi])
axis image off
hold on
fill([xpos flip(xpos)],length(ops.yrange)+.5-[0 0 0 0],'r', 'EdgeColor','none')
hold off
an.String = sprintf('');
delete(an2)
delete(anScale)
delete(anTitle)
frame = getframe;
for k = 1:round(sr * speedUp_dark_triggered)
    writeVideo(vid,frame)
end
close(vid)
close(gcf)

% make movie
mov_short = mov;
running_scaled = (running-min(running)) ./ (max(running)-min(running));
running_scaled = round(running_scaled.*length(ops.yrange));
xpos = [-round(length(ops.xrange)/40) 0];
bad = find(ops.badframes);
mov_short(:,:,bad) = [];
running_scaled(bad) = [];
running_binary(bad) = [];
meanIm_noRunning = mean(mov_short(:,:,~running_binary),3);
mov_short = double(mov_short) - meanIm_noRunning;

figure
meanIm = mean(mov_short, 3);
imagesc(meanIm, [-1 1].*max(meanIm(:)))
colormap(cm)
axis image off
hold on
fill([xpos flip(xpos)],length(ops.yrange)+.5-[0 0 0 0],'r', 'EdgeColor','none')
annotation('rectangle',[.825 .52 .0554 .14],'Color','b','LineWidth',2)
annotation('rectangle',[.31 .671 .0893 .067],'Color','r','LineWidth',2)
% annotation('rectangle',[.707 .64 .05 .0624],'Color','r','LineWidth',2)
% annotation('rectangle',[.604 .826 .075 .0716],'Color','r','LineWidth',2)
annotation('rectangle',[.348 .578 .0733 .067],'Color','r','LineWidth',2)
% annotation('rectangle',[.211 .49 .0482 .09],'Color','r','LineWidth',2)
title('Mean image')

mov_short(:,:,1:2:end) = [];
tim = (0:size(mov_short,3)-1)./(sr/2);
running_scaled(1:2:end) = [];

% clims = prctile(mov_short(:), [0 100]);
clims = [-1 1] .* max(abs(mov_short(:))) * 0.6;
vid = VideoWriter(outputFile_dark,'MPEG-4');
vid.FrameRate = 1 / median(diff(frameTimes)) * speedUp_dark / 2;
open(vid);
% 2d gauss to convolve movie with
win = normpdf(-3:3,0,1);
win = win'*win;
figure
colormap(cm)
annotation('rectangle',[.825 .52 .0554 .14],'Color','b','LineWidth',2)
annotation('rectangle',[.31 .671 .0893 .067],'Color','r','LineWidth',2)
% annotation('rectangle',[.707 .64 .05 .0624],'Color','r','LineWidth',2)
% annotation('rectangle',[.604 .826 .075 .0716],'Color','r','LineWidth',2)
annotation('rectangle',[.348 .578 .0733 .067],'Color','r','LineWidth',2)
% annotation('rectangle',[.211 .49 .0482 .09],'Color','r','LineWidth',2)
an = annotation('textbox',[.17 .8 .1 .1],'String', ...
    sprintf('%.0f s',0),'EdgeColor','none','FontSize',14, ...
    'FontWeight','bold','HorizontalAlignment','right');
annotation('textbox',[0.73 0.2 0.12 0.05],'String','20 \mum', ...
    'EdgeColor','none','FontSize',12);
annotation('line',[0.7 0.868],[0.25 0.25],'LineWidth',2,'Color','k');
annotation('textbox',[0.4 0.83 0.22 0.06],'String','In darkness', ...
    'EdgeColor','none','FontSize',14,'FontWeight','bold');
for k = 1:length(running_scaled)
%     im = conv2(double(mov_short(:,:,k)),win,'same');
    im = mov_short(:,:,k);
    imagesc(im, clims)
%     colormap gray
    axis image off
    hold on
    fill([xpos flip(xpos)],length(ops.yrange)+.5-[0 0 [1 1].*running_scaled(k)],'r', 'EdgeColor','none','FaceColor','r')
    hold off
    an.String = sprintf('%.0f s',tim(k));
    txt = text(6,140,'running speed');
    txt.Color='r';
    txt.Rotation=90;
    txt.FontSize=12;
    txt.FontWeight='bold';
    frame = getframe;
    writeVideo(vid,frame)
end
close(vid)
close(gcf)

% plot mean image
meanIm = mean(mov_short,3);
figure
imagesc(meanIm,prctile(double(meanIm(:)), [1 99]))
colormap gray
axis image off
% plot bouton masks, colour codes for correlation with running
indNaN = any(isnan([running',meta.F_final]),2);
rhos = corr(running(~indNaN)', meta.F_final(~indNaN,:));
maxi = max(abs(rhos));
weights = linspace(0,1,11)';
cm = [[0 0 1].*flip(weights)+[1 1 1].*weights; ...
    [1 1 1].*flip(weights)+[1 0 0].*weights];
cm(length(weights),:) = [];
base = -(maxi + 2*maxi/(size(cm,1)));
masks = ones(Ly,Lx) .* base;
for j = 1:length(meta.targetFrameROI)
    masks(meta.targetFrameROI{j}) = rhos(j);
end
masks = masks(ops.yrange, ops.xrange);
figure
imagesc(masks,[base maxi])
colormap([[0 0 0]; cm])
axis image off
colorbar('Ticks',floor(maxi*20)/20 .* [-1 0 1],'Limits',[-maxi maxi])

%% Make movie for gratings
% SS070, 2016-10-18, plane 9, exp. 1
% movieFile = 'J:\DATA\registered\M160923_SS070\2016-10-18\1_2_3_4\interpolated\M160923_SS070_2016-10-18_1_2_3_4_plane9.bin';
% regFile = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\F\M160923_SS070\2016-10-18\1_2_3_4\regops_M160923_SS070_2016-10-18.mat';
% metaFile_gratings = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\InfoStructs\M160923_SS070\2016-10-18\1\2016-10-18_1_M160923_SS070_2P_plane009_ROI.mat';
% tuningFile = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\boutons\nonVisualEffects\pupil\tuning_prefDirSigmaDIFixed_isTuned.mat';
% outputFile_gratings = 'C:\STORAGE\OneDrive - University College London\Presentations\boutonMedia\SS070_2016-10-18_1_gratings.avi';
% plane = 9;
% exp_gratings = 1;
% subject = 5;
% planeID = 3;
% arrowShift = [20 45];
% SS069, 2016-10-21, plane 6, exp. 1
% movieFile = 'J:\DATA\registered\M160923_SS069\2016-10-21\1_2_3\interpolated\M160923_SS069_2016-10-21_1_2_3_plane6.bin';
% regFile = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\F\M160923_SS069\2016-10-21\1_2_3\regops_M160923_SS069_2016-10-21.mat';
% metaFile_gratings = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\InfoStructs\M160923_SS069\2016-10-21\1\2016-10-21_1_M160923_SS069_2P_plane006_ROI.mat';
% tuningFile = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\boutons\nonVisualEffects\pupil\tuning_prefDirSigmaDIFixed_isTuned.mat';
% outputFile_gratings = 'C:\STORAGE\OneDrive - University College London\Presentations\boutonMedia\SS069_2016-10-21_1_plane6_gratings.avi';
% plane = 6;
% exp_gratings = 1;
% subject = 4;
% planeID = 1;
% arrowShift = [49 20];
% movLimits = [300 2560];
% SS075, 2017-09-13, plane 5, exp. 4
movieFile = 'J:\Data\registered\M170821_SS075\2017-09-13\1_2_4\interpolated\M170821_SS075_2017-09-13_1_2_4_plane5.bin';
regFile = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\F\M170821_SS075\2017-09-13\1_2_4\regops_M170821_SS075_2017-09-13.mat';
metaFile_gratings = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\InfoStructs\M170821_SS075\2017-09-13\3\2017-09-13_3_M170821_SS075_2P_plane005_ROI.mat';
tuningFile = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\boutons\nonVisualEffects\pupil\tuning_prefDirSigmaDIFixed_isTuned.mat';
outputFile_gratings = 'C:\STORAGE\OneDrive - University College London\Presentations\boutonMedia\M170821_SS075_2017-09-13_3_plane5_gratings';
plane = 5;
exp_gratings = 3;
subject = 7;
planeID = 1;
arrowShift = [25 30];
movLimits = [400 2000];

speedUp_gratings = 5;
smoothStd = .5;

% load movie
data = load(regFile);
ops = data.ops1{plane};
frames = ops.Nframes;
ops1{1} = ops;
ds = ops.DS_allPlanes;
fr = sum(frames(1:exp_gratings-1))+1:sum(frames(1:exp_gratings));
if length(size(ops.CorrFrame_allPlanes)) == 2
    invalid = isnan(ops.CorrFrame_allPlanes(:,plane));
    ds(invalid,:,:) = [];
    ds = ds(fr,:,:);
    ops1{1}.CorrFrame_allPlanes(invalid,:) = [];
    ops1{1}.CorrFrame_allPlanes = ops1{1}.CorrFrame_allPlanes(fr,:);
else
    pl = ops.planesToInterpolate == plane;
    invalid = all(isnan(ops.CorrFrame_allPlanes(:,pl,pl)),2);
    ds(invalid,:,:,:,:) = [];
    ds = ds(fr,:,:,:,:);
    ops1{1}.CorrFrame_allPlanes(invalid,:,:) = [];
    ops1{1}.CorrFrame_allPlanes = ops1{1}.CorrFrame_allPlanes(fr,:,:);
end
ops1{1}.DS_allPlanes = ds;
ops1{1}.usedPlanes = ops1{1}.usedPlanes(fr,:);
ops1{1}.badframes = false(1,length(fr));
ops1 = getRangeNoOutliers(ops,ops1);
ops1{1}.yrange = ops.yrange;
ops1{1}.xrange = ops.xrange;
ops = ops1{1};
Ly = ops.Ly;
Lx = ops.Lx;
fid = fopen(movieFile, 'r');
if exp_gratings > 1
    fseek(fid, Ly*Lx*sum(frames(1:exp_gratings-1))*2,'bof');
end
mov = fread(fid, Ly*Lx*frames(exp_gratings), '*int16');
mov = reshape(mov, Ly, Lx, []);
mov = mov(ops.yrange, ops.xrange, :);
% get frame times and stimulus info
data = load(metaFile_gratings);
meta = data.meta;
[~, stimSeq, stimMatrix, frameTimes] = ...
    ssLocal.getStimulusResponseInfo(meta);
[directions, blanks] = gratings.getOrientations(stimSeq);
sr = 1/median(diff(frameTimes));

stdSamples = round(smoothStd * sr);
convWindow = normpdf(-4*stdSamples:4*stdSamples, 0, stdSamples);
L = length(convWindow);
mov = cat(3, repmat(mean(mov(:,:,1:floor(L/2)),3),1,1,L), ...
    mov, repmat(mean(mov(:,:,end-floor(L/2):end),3),1,1,L));
for y = 1:size(mov,1)
    for x = 1:size(mov,2)
        mov(y,x,:) = conv(squeeze(mov(y,x,:)), convWindow, 'same');
    end
end
mov(:,:,[1:L, end-L+1:end]) = [];
% get tuning info
data = load(tuningFile);
tuning = data.tuning;
inds = tuning(subject).plane(planeID).isTuned==1;
tunedUnits = tuning(subject).plane(planeID).cellIDs(inds);
prefDirs = cat(2,tuning(subject).plane(planeID).cond(1).cell(inds).parameters);
prefDirs = round(prefDirs(1,:));

% make movie
clims = prctile(double(mov(:)), [1 99]);
bad = find(ops.badframes);
mov(:,:,bad) = [];
stimMatrix(:,bad) = [];
vid = VideoWriter(outputFile_gratings,'MPEG-4');
vid.FrameRate = 1 / median(diff(frameTimes)) * speedUp_gratings;
open(vid);
figure
imagesc(mov(:,:,1))
axis image off
ax = gca;
ax.Units = 'points';
axPos = ax.Position;
arrowPos = axPos([1 2]) + [arrowShift(1) axPos(4)-arrowShift(2)];
an2 = annotation('textbox',[0.73 0.16 0.12 0.05],'String','20 \mum', ...
    'EdgeColor','none','FontSize',12,'Color','w');
anTitle = annotation('textbox',[0.15 0.17 0.22 0.06],'String','Gratings', ...
    'EdgeColor','none','FontSize',14,'FontWeight','bold','Color','w');
anScale = annotation('line',[0.7 0.868],[0.21 0.21],'LineWidth',2,'Color','w');
for k = movLimits(1):min(movLimits(2),size(mov,3))
    imagesc(mov(:,:,k), clims)
    colormap gray
    axis image off
    stim = find(stimMatrix(:,k));
    if isempty(stim) || stim > size(directions,1)
%         a = annotation('rectangle');
%         a.Units = 'points';
%         a.Position = [arrowPos(1)-15 arrowPos(2)-15 30 30];
%         a.FaceColor = 'w';
%         a.Color = 'none';
    else
        angle = directions(directions(:,2)==stim,1)/180*pi;
        a = annotation('arrow');
        a.Units = 'points';
        a.Position = [arrowPos(1)-cos(angle)*15, arrowPos(2)+sin(angle)*15, ...
            cos(angle)*30, -sin(angle)*30];
        a.LineWidth = 4;
        a.Color = 'w';
        a.HeadWidth = 20;
    end
    frame = getframe;
    writeVideo(vid,frame)
    if exist('a','var')
        delete(a)
    end
end
close(vid)
close(gcf)

% plot mean image
meanIm = mean(mov_short,3);
figure
imagesc(meanIm,[min(meanIm(:)), 0.8*max(meanIm(:))])
colormap gray
axis image off
% plot bouton masks, colour codes for pref. direction
masks = ones(Ly,Lx).*-1;
for j = 1:length(tunedUnits)
    masks(meta.targetFrameROI{tunedUnits(j)}) = prefDirs(j);
end
masks = masks(ops.yrange, ops.xrange);
figure
imagesc(masks, [-1 360])
colormap([0 0 0; hsv(361)])
axis image off
colorbar('Ticks',0:90:360,'Limits',[0 360])