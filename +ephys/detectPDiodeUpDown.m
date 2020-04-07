


function flipTimes = detectPDiodeUpDown(sig, Fs, threshUp, threshDown)
% version of photodiode detection that works for the photodiode setting
% (called "flickergrey" in the config files) in which the syncsquare goes
% from white to black or black to white each frame, but in between stimuli
% it goes to a middle grey level. The point is that if your stimulus ended
% on a black frame, you wouldn't know when it turned off. So a "flip time"
% is defined as when a null frame crosses threshUp upwards, an up frame crosses
% threshUp downwards, or a down frame crosses threshDown upwards.
% Edit 2014-12-15: for reasons unclear, sometimes it goes from grey
% straight to black. So adding threshDownDownwards.

threshUpUpwards = find(sig(1:end-1)<threshUp & sig(2:end)>=threshUp);
threshUpDownwards = find(sig(1:end-1)>=threshUp & sig(2:end)<threshUp);
threshDownUpwards = find(sig(1:end-1)<threshDown & sig(2:end)>=threshDown);
threshDownDownwards = find(sig(1:end-1)>=threshDown & sig(2:end)<threshDown);

% threshUpUpwards commonly follows threshDownUpwards by a small or
% nonexistent lag, here assumed to be <4ms. So all of these threshold
% crossings together constitute the flipTimes, but many of them are doubled
% up and we can safely do away with these doublings. 

flipSamps = sort([threshUpUpwards; threshUpDownwards; threshDownDownwards; threshDownUpwards]);
flipTimes = flipSamps/Fs;

doubleThresh = 0.006; % seconds
dFlipTimes = diff([0; flipTimes]);
flipTimes = flipTimes(dFlipTimes > doubleThresh);

