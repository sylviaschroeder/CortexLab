function [pupilDiam, pupilDevFromStandard] = getPupilDiam(pupilData, smoothSpan)

if nargin < 2
    smoothSpan = 5;
end

% areaSmoothed = smooth(pupilData.area, smoothingSpan, 'rlowess');
pupilDiam = sqrt(4 * pupilData.area / pi);
pupilDiam(isnan(pupilData.area) | pupilData.blink | ~pupilData.goodFit) = NaN;
if isfield(pupilData, 'blinkSoft')
    pupilDiam(pupilData.blinkSoft | pupilData.blinkManual) = NaN;
end
pupilDiam = medfilt1(pupilDiam, smoothSpan);

if nargout > 1
    xMode = mode(round(pupilData.x/2)*2);
    yMode = mode(round(pupilData.y/2)*2);
    x = pupilData.x - xMode;
    y = pupilData.y - yMode;
    pupilDevFromStandard = sqrt(x.^2+y.^2);
    pupilDevFromStandard(isnan(pupilData.area) | pupilData.blink | ...
        ~pupilData.goodFit) = NaN;
end