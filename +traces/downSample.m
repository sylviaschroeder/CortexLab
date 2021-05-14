function convTrace = downSample(trace, time, sigma, newTime, keepNaNs)

if nargin < 5
    keepNaNs = true;
end

dt = median(diff(time));
sigSamples = round(sigma / dt);
win = normpdf(-4*sigSamples : 4*sigSamples, 0, sigSamples);
invalid = isnan(trace);
invalidBins = histcounts(time(invalid), newTime) > 0;
halfWin = (length(win)-1)/2;
padded = [ones(ceil(halfWin),1) .* nanmean(trace(1 : ceil(halfWin))); ...
    trace; ...
    ones(floor(halfWin),1) .* nanmean(trace(end-floor(halfWin) : end))];
convTrace = conv(padded, win, 'valid');
invalid = isnan(convTrace);
convTrace = interp1(time(~invalid), convTrace(~invalid), newTime, 'pchip');

if keepNaNs
    convTrace(invalidBins) = NaN;
end