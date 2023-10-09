function [frameTimes tl_flag] = getFrameTimes(videoPath, timeLinePath)

% This function will return the frametimes of the eye-tracking movie
% The frame times will be aligned to the Timeline time axis
% If corresponding Timeline file does not exist, frame times will be
% relative only.

% May 2014 - MK Created

warning off
data = load([videoPath '.mat']);
warning on
eyeLog = data.eyeLog;

data = load(timeLinePath);
Timeline = data.Timeline;
    
vReader = VideoReader([videoPath '.mj2']);
nFrames = vReader.NumFrames;
fprintf('There are %d frames in the video file\n', nFrames);
fprintf('There are %d timestamps in the log file\n', length(eyeLog.TriggerData));

% find the last ExpStart-ExpEnd pair
endInd = [];
startInd = [];
for iEvent = length(eyeLog.udpEvents):-1:1
    if strfind(eyeLog.udpEvents{iEvent}, 'ExpEnd')
    %if ~isempty(strfind(eyeLog.udpEvents{iEvent}, 'ExpEnd')) ...
    %        || ~isempty(strfind(eyeLog.udpEvents{iEvent}, 'ExpInterrupt'))
        endInd = iEvent;
        break;
    end
end
for iEvent = (endInd-1):-1:1
    if strfind(eyeLog.udpEvents{iEvent}, 'ExpStart')
        startInd = iEvent;
        break;
    end
end

udpEvents = eyeLog.udpEvents(startInd:endInd);
udpEventTimes = eyeLog.udpEventTimes(startInd:endInd);

nEvents = length(udpEventTimes);

if exist('Timeline','var')
    tl_flag = 1;
else
    tl_flag = 0;
    warning('Timeline was not found!!!');
end

if  tl_flag && ~isequal(nEvents, Timeline.mpepUDPCount)
    warning('Number of UDP events logged by Timeline and by EyeCamera is different. Something must be wrong!');
    nEvents = min(nEvents, Timeline.mpepUDPCount);
end

% converting absolute times to times in seconds
eyeTimes = nan(nEvents, 1);
for iEvent=1:nEvents
    eyeTimes(iEvent) = datenum(udpEventTimes{iEvent})*(24*60*60);
end
frameTimes = nan(nFrames, 1);
for iFrame = 1:nFrames
    if iFrame<=length(eyeLog.TriggerData)
        frameTimes(iFrame) = datenum(eyeLog.TriggerData(iFrame).AbsTime)*(24*60*60);
    else
        frameTimes(iFrame)=NaN;
    end
end

if tl_flag
    tlTimes = Timeline.mpepUDPTimes(1:nEvents);
    nEvents2Discard = 1;
    idx = nEvents2Discard+1:nEvents;
    timeDiff = median(eyeTimes(idx) - tlTimes(idx));
    frameTimes = frameTimes - timeDiff;    
end

return;

%% =============some plotting for debugging purposes=========
eyeTimes = eyeTimes - timeDiff;

figure
stem(eyeTimes, ones(nEvents, 1), 'b');
hold on;
stem(tlTimes, ones(nEvents, 1), 'r:');
legend('eyeUDPs', 'tlUDPs');
xlabel('time [seconds]');
title('UDP messges Timing (aligned)');

figure
dd = diff(eyeTimes-tlTimes);
plot(dd(nEvents2Discard+1:end))
title('UDP timing jitter (eyeCamera - Timeline)');
xlabel('UDP message number');
ylabel('Time difference [sec]');

figure
hist(dd(nEvents2Discard+1:end), 20);
title('Time jitter histogram');
xlabel('Time difference [sec]');


