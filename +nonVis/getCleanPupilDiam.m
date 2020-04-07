function [pupilDiameter, pupilTime] = getCleanPupilDiam(pupilData, pupilTime)

pupilTime(length(pupilData.x)+1:end) = [];
pupilDiam = sqrt(4 * pupilData.area / pi)';
badTimes = [false diff(pupilDiam) < -10 * nanstd(diff(pupilDiam))];
badTimes = repmat(find(badTimes), 3, 1);
badTimes = bsxfun(@plus,badTimes,(-1:1)');
badTimes = unique(badTimes(:));
tmp = false(size(pupilDiam));
tmp(badTimes) = true;
badTimes = tmp;
ind = ~badTimes' & ~pupilData.blink & pupilData.goodFit;
pupilDiameter = interp1(pupilTime(ind), pupilDiam(ind), pupilTime, 'pchip');