function [stimOnsets, stimOffsets, stimIDs] = mpepPhotodiode(diode, Fs, ...
    protocol, timeline)
% Tries to figure out when stimuli started and ended in the photodiode
% trace
%
% - diode is the photodiode trace
% - datFs is the sampling frequency of the photodiode signal

% assumption: syncSquareType is "flickergray"

pdThreshUp = 0.7;
pdThreshDown = 0.1;


diode = medfilt1(diode,10);
    
flipTimes = ephys.detectPDiodeUpDown(diode, Fs, pdThreshUp, pdThreshDown);
disp(['   found ' num2str(length(flipTimes)) ' photodiode flips']);

diffFlips = diff([0; flipTimes; flipTimes(end)+1]); % will be 15 or 18 ms for an up or down phase of the flicker
st = find(diffFlips>0.02);
allStimOn = flipTimes(st(1:end-1));
allStimOff = flipTimes(st(2:end)-1);
% actualDurs = allStimOff-allStimOn; actualDurs = actualDurs(:);

% analyze mpep udp events
tlBlockStarts = [];
tlBlockEnds = [];
tlStimStarts = [];
tlStimEnds = [];
for mp = 1:timeline.mpepUDPCount
    if strcmp(timeline.mpepUDPEvents{mp}(1:10), 'BlockStart')
        tlBlockStarts(end+1) = timeline.mpepUDPTimes(mp);
    elseif strcmp(timeline.mpepUDPEvents{mp}(1:8), 'BlockEnd')
        tlBlockEnds(end+1) = timeline.mpepUDPTimes(mp);
    elseif strcmp(timeline.mpepUDPEvents{mp}(1:9), 'StimStart')
        tlStimStarts(end+1) = timeline.mpepUDPTimes(mp);
    elseif strcmp(timeline.mpepUDPEvents{mp}(1:7), 'StimEnd')
        tlStimEnds(end+1) = timeline.mpepUDPTimes(mp);
    end
end

% The first tlStimStarts often shows up well before it should. Correct
% this here, if necessary.
%     tlStimStarts(1) = tlStimStarts(1)+0.25;

% try to match the photodiode times to the timeline times
% assume that the first allStimOn is the same as the first tlStimStarts
matchedStimOn = zeros(size(tlStimStarts));
matchedStimOff = zeros(size(tlStimStarts));
targetDiffOn = zeros(size(tlStimStarts));
targetDiffOff = zeros(size(tlStimStarts));
matchedStimOn(1) = allStimOn(1);
for ss = 2:length(tlStimStarts)
    targetTime = tlStimStarts(ss)-tlStimStarts(1)+matchedStimOn(1);
    matchedStimOn(ss) = allStimOn(abs(allStimOn-targetTime)==min(abs(allStimOn-targetTime)));
    targetDiffOn(ss) = matchedStimOn(ss)-targetTime;
end
for se = 1:length(tlStimEnds)
    targetTime = tlStimEnds(se)-tlStimStarts(1)+matchedStimOn(1);
    matchedStimOff(se) = allStimOff(abs(allStimOff-targetTime)==min(abs(allStimOff-targetTime)));
    targetDiffOff(se) = matchedStimOff(se)-targetTime;
end
stimOnsets = matchedStimOn;
stimOffsets = matchedStimOff;
if max(abs(targetDiffOn))>0.2 || max(abs(targetDiffOff))>0.2
    disp(' warning, photodiode signals might be wrong.');
    keyboard;
    
    figure
    hold on
    plot((1:length(diode))/Fs, diode)
    plot(allStimOn, ones(size(allStimOn))+500, 'g.')
    plot(allStimOff, ones(size(allStimOff)), 'r.')
    plot(matchedStimOn, ones(size(matchedStimOn))-500, 'm.')
    plot(matchedStimOff, ones(size(matchedStimOff))-1000, 'k.')
    plot(tlStimStarts-tlStimStarts(1)+allStimOn(1), ones(size(tlStimStarts))-1500, '.', 'Color', [1 0.5 0])
    plot(tlStimEnds-tlStimStarts(1)+allStimOn(1), ones(size(tlStimEnds))-2000, '.', 'Color', [0 0.5 1])
    legend({'photodiode signal', 'all stimulus onsets', 'all stimulus offsets', 'onsets matched to TL', 'offsets matched to TL', 'TL start times', 'TL end times'})
else
    disp(' photodiode seems to match TL, at least to some approximation.');
end

seqnums = protocol.seqnums;
% numStims = protocol.npfilestimuli;
sInd = 1;
stimIDs = zeros(1, numel(seqnums));
for r = 1:size(seqnums, 2)
    for s = 1:size(seqnums,1)
        stimIDs(sInd) = find((seqnums(:,r)-(r-1)*size(seqnums,1))==s);
        sInd = sInd+1;
    end
end