function wheelData = getRunningSpeed_wheel(wheelMeasures, time, smoothingStd)

if nargin < 3
    smoothingStd = .5; % in sec
end

maxValue = max(wheelMeasures);
diffs = abs(diff(wheelMeasures));
if max(diffs) > 0.5 * maxValue
    ind = wheelMeasures > 0.5 * maxValue;
    wheelMeasures(ind) = wheelMeasures(ind) - maxValue;
    diffs = diff(wheelMeasures);
end
if sum(diffs<0) > sum(diffs>0)
    diffs = -diffs;
end
timeStep = median(diff(time));

if smoothingStd > 0
    smoothingSamples = round(smoothingStd / timeStep);
    win = normpdf(-3*smoothingSamples:3*smoothingSamples, 0, ...
        smoothingSamples);
    wheelData.total = conv(diffs ./ timeStep, win, 'same');
else
    wheelData.total = diffs ./ timeStep;
end

wheelData.t = time(1:end-1);