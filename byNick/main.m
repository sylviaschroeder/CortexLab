lickInd = find(strcmp({Timeline.hw.inputs.name},'piezoLickDetector'));
sig = Timeline.rawDAQData(:,lickInd);
Fs = Timeline.hw.daqSampleRate;
times = Timeline.rawDAQTimestamps;

[lickTimes, filtSig] = piezoLickDetector(times, sig, Fs);

figure; 

plot(times, sig);
hold on; 
plot(times, filtSig+max(sig)-min(filtSig));
[xx,yy] = rasterize(lickTimes);
plot(xx,yy-min(sig),'k')