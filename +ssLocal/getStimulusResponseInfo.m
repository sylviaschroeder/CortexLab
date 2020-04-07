function [stimTimes, stimSequence, stimMatrix, frameTimes, ...
    samplingRate] = getStimulusResponseInfo(meta)

% Get stimulus information
stimTimes = ppbox.getStimTimes(meta);
d = sscanf(meta.expDate, '%d-%d-%d');
if d(1)>2017 || (d(1)==2017 && d(2)>=9)
    stimSequence = ppbox.getStimSequence(meta.subject, meta.expDate, meta.exp);
else
    stimSequence = ppbox.getStimSequence(meta);
end

% Get recorded frames and times
% frames = loadArr(fullfile(meta.folderProcessed, meta.basenameRegistered));
frameTimes = ppbox.getFrameTimes(meta);
samplingRate = 1 / median(diff(frameTimes));

% Match stimulus and frame times
stimMatrix = ppbox.buildStimMatrix(stimSequence, stimTimes, frameTimes);