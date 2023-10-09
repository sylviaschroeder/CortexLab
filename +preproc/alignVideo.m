
function alignVideo(mouseName, thisDate, expNum, movieName, varargin)
% function alignVideo(mouseName, thisDate, expNum, movieName[, params])

folderVideo = 'Z:\UCLData\Ephys_Task\Subjects';

reposName = 'master';
timelineExpNums = expNum;
tlSyncName = 'camSync';
silentMode = true;
recompute = false;
nFramesToLoad = 3000;

switch movieName
    case 'eye'
        strobeName = 'eyeCameraStrobe';
    case 'face'
        strobeName = 'faceCamStrobe';
    otherwise
        strobeName = 'unknown';
end

if ~isempty(varargin)
    params = varargin{1};
    if isfield(params, 'reposName')
        reposName = params.reposName;
    end
    if isfield(params, 'timelineExpNums')
        timelineExpNums = params.timelineExpNums;
    end
    if isfield(params, 'strobeName')
        strobeName = params.strobeName;
    end
    if isfield(params, 'tlSyncName')
        tlSyncName = params.tlSyncName;
    end
    if isfield(params, 'silentMode')
        silentMode = params.silentMode;
    end
    if isfield(params, 'recompute')
        recompute = params.recompute;
    end
    if isfield(params, 'nFramesToLoad')
        nFramesToLoad = params.nFramesToLoad;
    end
end

assert(~strcmp(strobeName, 'unknown'), 'please provide params.strobeName');

%%

% mouseName = 'Dam';
% thisDate = '2016-08-26';
% expNum = 1;
% movieName = 'eye';

movieDir = fullfile(folderVideo, mouseName, thisDate, num2str(expNum));
intensFile = fullfile(movieDir, [movieName '_avgIntensity.mat']);

if recompute && exist(intensFile, 'file')
    delete(intensFile);
end

if ~exist(intensFile, 'file')            
    fprintf(1, 'computing average intensity of first/last frames...\n');
    if silentMode
        ROI = avgMovieIntensity(movieDir, movieName, [], true, [], [], nFramesToLoad);
    else
        ROI = avgMovieIntensity(movieDir, movieName, [], true, 'ask', [], nFramesToLoad);
    end
end

fprintf(1, 'loading avg intensity\n');
load(intensFile)

%% first detect the pulses in the avgIntensity trace

expectedNumSyncs = numel(timelineExpNums)*2; % one at the beginning and end of each timeline file

vidIntensThresh = [15 20];
[intensTimes, intensUp, intensDown] = schmittTimes(1:numel(avgIntensity), avgIntensity, vidIntensThresh);
attemptNum = 1; loadAttemptNum = 1;
while(numel(intensDown)~=expectedNumSyncs)
    % try some different approaches to get the right threshold
    % automatically...
    switch attemptNum
        case 1
            vidIntensThresh = min(avgIntensity)*[1.2 1.4];
        case 2
            intensMed = median(avgIntensity);
            intensMin = min(avgIntensity);
            vidIntensThresh = intensMin+(intensMed-intensMin)*[0.4 0.6];
        case 3
            vidIntensThresh = intensMin+(intensMed-intensMin)*[0.15 0.25];
        otherwise
            switch loadAttemptNum
                case 1
                    fprintf(1, 'trying to load more frames...\n')
                    avgMovieIntensity(movieDir, movieName, [], true, ROI, [], 10000);
                    load(intensFile)
                    attemptNum = 0;
                case 2
                    fprintf(1, 'trying to load all frames...\n')
                    avgMovieIntensity(movieDir, movieName, [], true, ROI);
                    load(intensFile)
                    attemptNum = 0;
                otherwise
                    fprintf(1, 'cannot find a threshold that works. You tell me...\n');
                    figure; plot(avgIntensity);
                    keyboard
            end
            loadAttemptNum = loadAttemptNum+1;
    end
    
    [intensTimes, intensUp, intensDown] = schmittTimes(1:numel(avgIntensity), avgIntensity, vidIntensThresh);
    attemptNum = attemptNum +1;
end
assert(numel(intensDown)==expectedNumSyncs, 'could not find correct number of syncs');
fprintf(1, 'found the sync pulses in the video\n');

%%
% usually these are TTL and pretty much anything between 0 adn 5 will work
% but here I use a small value because for the whisker camera I didn't have
% TTL so this works for that too.
tlStrobeThresh = [0.08 0.15]; 


tlSyncThresh = [2 3];


tVid = NaN(size(avgIntensity));

for tInd = 1:numel(timelineExpNums)
    
    fprintf(1, 'loading timeline\n');
    load(fullfile(movieDir, sprintf('%s_%d_%s_Timeline.mat', thisDate, ...
        timelineExpNums(tInd), mouseName)));

    %check strobe counts against the number of frames between sync events that
    %were actually in the video

    % find the timeline samples where cam sync pulses started (went
    % from 0 to 5V)
    syncIndex = find(strcmp({Timeline.hw.inputs.name}, tlSyncName));
    tlSync = Timeline.rawDAQData(:,syncIndex);
    [~, tlSyncOnSamps, ~] = schmittTimes(1:numel(tlSync), tlSync, tlSyncThresh);

    % find the strobe times for the camera
    strobeIndex = find(strcmp({Timeline.hw.inputs.name}, strobeName));
    tlStrobe = Timeline.rawDAQData(:,strobeIndex);
    [~,strobeSamps,~] = schmittTimes(1:numel(tlStrobe), tlStrobe, tlStrobeThresh);
    
    vidSyncOnFrames = intensDown([1 2]+(tInd-1)*2);

    numStrobesFoundBetweenSyncs = sum(strobeSamps>=tlSyncOnSamps(1) & strobeSamps<tlSyncOnSamps(2));
    numFramesFoundBetweenSyncs = diff(vidSyncOnFrames);

    framesMissed = numStrobesFoundBetweenSyncs-numFramesFoundBetweenSyncs;

    %unique(diff(strobeSamps))/Timeline.hw.daqSampleRate

    framesMissedPerSec = framesMissed/(diff(tlSyncOnSamps)/Timeline.hw.daqSampleRate);
    
    fprintf(1, 'missed %d frames, for %.2f frames missed/sec\n', framesMissed, framesMissedPerSec);
    if abs(framesMissed)<=2 && abs(framesMissed)>0
        fprintf(1, 'values of +/-2 are normal, can happen when you get exposures across part of the cam sync onsets\n');
    elseif abs(framesMissed)>2
        fprintf(1, 'too many missed frames! Figure out how to reduce it. Will interpolate frame times linearly.\n');
    end
    
    % get the actual timestamps for the video in question. 
    % need to make this a function

    tt = Timeline.rawDAQTimestamps;
    
    if framesMissed==0
        % best case, use the strobe times exactly
        tVid(vidSyncOnFrames(1)+1:vidSyncOnFrames(2)) = ...
            tt(strobeSamps(strobeSamps>=tlSyncOnSamps(1) & strobeSamps<tlSyncOnSamps(2)));
    else
        % otherwise interpolate from first to last
        tVid(vidSyncOnFrames(1)+1:vidSyncOnFrames(2)) = ...
            linspace(tt(strobeSamps(find(strobeSamps>=tlSyncOnSamps(1),1))), ...
            tt(strobeSamps(find(strobeSamps<tlSyncOnSamps(2),1, 'last'))), ...
            numel(vidSyncOnFrames(1)+1:vidSyncOnFrames(2)));
    end
    
    % fill in the timestamps before and after sync pulses
    vidFs = 1/(mean(diff(strobeSamps))/Timeline.hw.daqSampleRate);
    
    if tInd==1
        % only fill in before sync pulses if it's the first timeline file,
        % otherwise we'll overwrite the previous timestamps.
        tVid(1:vidSyncOnFrames(1)) = (-vidSyncOnFrames(1):-1)/vidFs+tVid(vidSyncOnFrames(1)+1);
    end
    tVid(vidSyncOnFrames(2)+1:end) = (1:numel(tVid(vidSyncOnFrames(2)+1:end)))/vidFs+tVid(vidSyncOnFrames(2));
    
end

saveName = fullfile(movieDir, [movieName '_timeStamps.mat']);
fprintf(1, 'saving to %s\n', saveName)

save(saveName, 'tVid', 'vidFs');