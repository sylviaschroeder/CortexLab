function [stimTraces, time] = getStimTraces(trace, traceTime, ...
    stimTimes, stimFrameTimes, timeBefore, timeAfter)

% create matrix of relevant time points [stimFrames x time]
frameDur = median(diff(stimFrameTimes));
time = (-round(timeBefore/frameDur) : round(timeAfter/frameDur)) .* frameDur;
timesAllFrames = reshape(([stimTimes.onset] + stimFrameTimes)', [], 1) + time;

% get calcium responses at each time point by interpolating calcium trace
stimTraces = interp1(traceTime, trace, timesAllFrames, 'pchip');
stimTraces(timesAllFrames < traceTime(1) | timesAllFrames > traceTime(end)) = NaN;