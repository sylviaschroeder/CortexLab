

function [lickTimes, filtSig] = piezoLickDetector(times, sig, Fs)

smSize = round(Fs/20);
lickSeparation = round(Fs/100);

thresh = 0.01;

filtSig = conv(abs(sig-mean(sig)), gausswin(smSize)./sum(gausswin(smSize)), 'same');

[PKS,LOCS]= findpeaks(filtSig, 'MinPeakDistance', lickSeparation, 'MinPeakProminence',thresh);

% lickTimes = times(LOCS(PKS>thresh));
lickTimes = times(LOCS);