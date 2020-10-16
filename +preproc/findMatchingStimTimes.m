function [stimOnTimes, stimOffTimes] = findMatchingStimTimes(blStimStarts, ...
    blStimEnds, blPhotodiode, tlPhotodiode, blockToTL)

lag = -max(blPhotodiode * blockToTL(1) - tlPhotodiode);
toPDTimeFrameLag = @(t)t*blockToTL(1) + lag;

stimOnTimes = follows(...
    toPDTimeFrameLag(blStimStarts), tlPhotodiode);

stimOffTimes = follows(...
    toPDTimeFrameLag(blStimEnds), tlPhotodiode);
end

function t = follows(a, b)
    n = numel(a);
    t = zeros(size(a));
    ti = t;
    for ii = 1:n
        ti(ii) = find(b > a(ii), 1);
        t(ii) = b(ti(ii));
    end

    d = t - a;
    range = max(d) - min(d);
    assert((range/mean(d)) < 4.2, 'delta range is much larger than the mean');
end