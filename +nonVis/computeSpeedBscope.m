function [mouseSpeed,alignedTime,mouseSpeedRot]= computeSpeedBscope(Timeline,options)

if nargin<2 || ~isfield(options, 'nDirs')
    options.nDirs=4;
end
if nargin<2 || ~isfield(options, 'forDir')
    options.forDir=1;
end
if nargin<2 || ~isfield(options, 'latDir')
    options.latDir=3;
end
if nargin<2 || ~isfield(options, 'rotDir1')
    options.rotDir1=2;
end
if nargin<2 || ~isfield(options, 'rotDir2')
    options.rotDir2=4;
end


if nargin<2 || ~isfield(options, 'ballToDegree')
    % divide ay/by by 8.3 to translate pixels to degrees of ball rotation (ball diameter is 20 cm)
    options.ballToDegree=1/8.3;
end
if nargin<2 || ~isfield(options, 'ballToRoom')
    % divide ax/bx by 53 to get from pixels to cm
    options.ballToRoom=1/53;
end
if nargin<2 || ~isfield(options, 'diffTimes')
    % time difference between acquisition
    options.diffTimes=15.6/1000;
end

nDirs=options.nDirs;
forDir=options.forDir;
latDir=options.latDir;
rotDir1=options.rotDir1;
rotDir2=options.rotDir2;
ballToDegree=options.ballToDegree;
ballToRoom=options.ballToRoom;
diffTimes=options.diffTimes;

nT=Timeline.ballUDPCount;
%nEvents=Timeline.mpepUDPCount;
events=Timeline.ballUDPEvents;
%ballTimes=Timeline.ballUDPTimes(1:nT);
balVec=nan(nT,nDirs);
realTime=nan(1,nT);
%realTime2=nan(1,nT);
%ballTime2=nan(1,nT);
%jumps=[];
%kExample=0;
for t=1:nT    
    temp = strsplit(events{t});
    balVec(t,:) = cellfun(@str2num, temp(2:end));
    realTime(t)=str2num(temp{1});
end

%% 

delayT=mean(Timeline.ballUDPTimes(round(nT/4):nT))-mean(realTime(round(nT/4):nT)/1000.);
alignedTime=realTime/1000.+delayT;

maxt=max(max((realTime-realTime(1))/1000.),max(Timeline.ballUDPTimes(1:nT)));
mint=min(min((realTime-realTime(1))/1000.),min(Timeline.ballUDPTimes(1:nT)));
figure;
hold on
plot(alignedTime,Timeline.ballUDPTimes(1:nT),'o')
plot(alignedTime,alignedTime,'k','linewidth',2.)


hold off
%% 


figure;
%subplot(3,1,3)
hist(diff(Timeline.ballUDPTimes(1:nT))*1000,[0:2:63])
xlim([0,15.6*4])
%% 
ballInterp=balVec;
for t=2:nT-1
    if balVec(t,forDir) ==0 && balVec(t,rotDir1)==0
        if (abs(balVec(t-1,forDir))>0  || abs(balVec(t-1,rotDir1))>0) && (abs(balVec(t+1,forDir))>0  || abs(balVec(t+1,rotDir1))>0) 

        ballInterp(t,forDir)=(balVec(t+1,forDir)+balVec(t-1,forDir))/2;
        ballInterp(t,rotDir1)=(balVec(t+1,rotDir1)+balVec(t-1,rotDir1))/2;
%         else
%             t
        end
    end
if balVec(t,latDir) ==0 && balVec(t,rotDir2)==0
        if (abs(balVec(t-1,latDir))>0  || abs(balVec(t-1,rotDir2))>0) && (abs(balVec(t+1,latDir))>0  || abs(balVec(t+1,rotDir2))>0)

        ballInterp(t,latDir)=(balVec(t+1,latDir)+balVec(t-1,latDir))/2;
        ballInterp(t,rotDir2)=(balVec(t+1,rotDir2)+balVec(t-1,rotDir2))/2;

        end
    end
end
%% 


figure;
hold on
plot(alignedTime,ballInterp(1:nT,forDir).*(ballToRoom/diffTimes),'b:')
plot(alignedTime,ballInterp(1:nT,rotDir1).*(ballToRoom/diffTimes),'r:')
plot(alignedTime,ballInterp(1:nT,latDir).*(ballToRoom/diffTimes),'g:')
plot(alignedTime,ballInterp(1:nT,rotDir2).*(ballToRoom/diffTimes),'k:')
plot(alignedTime,balVec(1:nT,forDir).*(ballToRoom/diffTimes),'b')
plot(alignedTime,balVec(1:nT,rotDir1).*(ballToRoom/diffTimes),'r')
plot(alignedTime,balVec(1:nT,latDir).*(ballToRoom/diffTimes),'g')
plot(alignedTime,balVec(1:nT,rotDir2).*(ballToRoom/diffTimes),'k')
xlim([alignedTime(1) alignedTime(nT)])
ylabel('velocity (cm/s)','fontsize',16)
xlabel('time (s)','fontsize',16)
legend({'forward (mouse #1)', 'rotation (mouse #1)', 'lateral (mouse #2)', 'rotation (mouse #2)'} ,'location','SouthWest')
set(gca,'fontsize',16)

%% 

figure;
subplot(2,1,1)
hold on
plot(balVec(1:nT,1))
plot(balVec(1:nT,2),'r')
plot(ballInterp(1:nT,forDir),'bo')
plot(ballInterp(1:nT,rotDir1),'mo')
subplot(2,1,2)
plot(balVec(1:nT,3),'g')
plot(balVec(1:nT,4),'k')
plot(ballInterp(1:nT,latDir),'g')
plot(ballInterp(1:nT,rotDir2),'k')

%% 

figure;
hold on
%plot(balVec(1:nT,1))
%plot(balVec(1:nT,2),'r')
plot(ballInterp(1:nT,forDir),'b')
plot(ballInterp(1:nT,rotDir1),'r')
%plot(balVec(1:nT,3),'g')
%plot(balVec(1:nT,4),'k')
plot(ballInterp(1:nT,latDir),'g')
plot(ballInterp(1:nT,rotDir2),'k')
%% 

mouseMov=sqrt(ballInterp(1:nT,forDir).^2+ballInterp(1:nT,latDir).^2);
mouseRot=(ballInterp(1:nT,rotDir1)+ballInterp(1:nT,rotDir2))./2;
mouseSpeed=mouseMov.*(ballToRoom/diffTimes);
mouseSpeedRot=mouseRot.*(ballToDegree/diffTimes);

%fprintf('shouldnt ball to degree and ball to room be dependent?\n')
%% 

smoothedVel=smooth(mouseSpeed,10);

figure;
hold on
plot(alignedTime,mouseSpeed)
plot(alignedTime,smoothedVel,'r')
%% 

figure;
hist(smoothedVel,[0:1:20])

sum(smoothedVel>1)/length(smoothedVel)
sum(smoothedVel>3)/length(smoothedVel)
%% 

figure;
plot(alignedTime,mouseSpeedRot)


%% 

