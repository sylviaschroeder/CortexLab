%%

%Written by Mario Dipoppa during Autumn 2017

%% clear workspace
clearvars;
close all;
clc;


%% load session metadata
load('metaData')
fprintf('There are %d datasets\n',length(metaData.fileName));
%% choose session
idataset = 1; %example good one: 10
%% load raw data
load(['rawData_' metaData.fileName{idataset}])
%% load eye position and time frames

time_eye = rawData.time_eye;
nt = length(time_eye); %length of eye position vector

eyepos = rawData.eyepos; %horizontal (1) and vertical (2) eye position

timeFrames = rawData.timeFrames; %this can be arbitrary and outside the function
nf = length(timeFrames); %length of calcium frame vector 




%% interpolate nan values in eye position

nanVals = zeros(1,nt);
nanFrames = find(isnan(mean(eyepos))); %find nan values of eye position
realFrames = find(~isnan(mean(eyepos))); %find correct values of eye position
nanVals(nanFrames) = ones(1,length(nanFrames));

if sum(isnan(time_eye))==length(time_eye)
    error('Bug in eye time frames\n')
end
nanSync = interp1(time_eye, nanVals, timeFrames);

%% interpolate eyeposition using time of calcium frames
eyepos_sync = nan(2,nf);
for ix = 1:2
    eyepos_real = eyepos(ix,:);
    eyepos_real(nanFrames) = interp1(realFrames, eyepos(ix,realFrames), nanFrames, 'pchip');
    eyepos_sync(ix,:) = interp1(time_eye, eyepos_real, timeFrames);
end



%% find saccades
StepThre=20;
saccades = findSaccades(eyepos_sync,nanSync,StepThre);


%% average eye position between saccades
[fix_mean,fix_pos] = modelFixation(saccades,eyepos_sync);
time_saccades = saccades;



%% create blue/white/red colormap
BlueWhiteRed = flipud(RedWhiteBlue);

%% plot eye movements
max_eye = max(eyepos,[],2);
min_eye = min(eyepos,[],2);
dt = 1000;
ncols = 4;
nrows = 3;

figure;
for it = 1:min(ncols*nrows,floor(nt/dt))
    subplot(nrows,ncols,it)
    ic = mod(it-1,ncols)+1;
    ir = floor((it-1)/ncols)+1;
    time_current = (it-1)*dt + (1:500);
    time_old = (1:time_current(1)-1);
    hold on
    plot(eyepos(1,time_old),eyepos(2,time_old),'Color',0.7.*ones(1,3))
    plot(eyepos(1,time_current),eyepos(2,time_current),'k')
    xlim([min_eye(1) max_eye(1)])
    ylim([min_eye(2) max_eye(2)])
    if ic ~= 1
        set(gca,'YColor','w','YTick',[])
    end
    if ir ~= nrows
        set(gca,'XColor','w','XTick',[])
    end
    if (ic == 1) && (ir == nrows)
        xlabel('x eye pos. (Deg)')
        ylabel('y eye pos. (Deg)')
    end
    axis image
    set(gca,'YDir','reverse')
end

%% plot distribution of eye position
axisNm = {'x','y'}; %axis name

ctrs = cell(1,2);
for ix = 1:2
    ctrs{ix} = (min_eye(ix):max_eye(ix));
end
[N,C] = hist3(eyepos',ctrs);
maxN = prctile(N(:),99);
figure;
hold on
imagesc(ctrs{1},ctrs{2},N',[-maxN maxN])
plot(fix_mean(1,:),fix_mean(2,:),'ko')
colormap(BlueWhiteRed);
axis image
set(gca,'YDir','reverse')
xlabel('x eye pos. (Deg)')
        ylabel('y eye pos. (Deg)')
%% plot time course of eye position
figure;
for ix = 1:2
    subplot(3,1,ix)
    hold on
    
    plot(time_eye,eyepos(ix,:),'r')
    plot(timeFrames,eyepos_sync(ix,:),'c')
    plot(timeFrames,fix_pos(ix,:),'k')
    xlim([-inf inf])
    ylim([-inf inf])
    ylabel([axisNm{ix} ' pos. (Deg)'])

end
subplot(3,1,3)
hold on
plot(time_eye,nanVals,'r')
plot(timeFrames,nanSync,'k')
xlim([-inf inf])
ylim([-0.1 1.1])
ylabel('Interpolated times')
xlabel('Time (s)')















