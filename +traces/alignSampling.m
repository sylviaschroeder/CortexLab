function newTraces = alignSampling(traces, time, delayIDs, delays)

newTraces = traces;
for k = length(delays) %1:length(delays)
    if delays(k) == 0
        continue
    end
    indUnits = find(delayIDs == k);
    
    trP = interp1(time + delays(k), traces(:, indUnits), time, 'pchip');
    trM = interp1(time + delays(k), traces(:, indUnits), time, 'makima');
    trS = interp1(time + delays(k), traces(:, indUnits), time, 'spline');
    for n = 1:10
        figure
        hold on
        plot(time + delays(k), traces(:, indUnits(n)), 'k')
        plot(time, trP(:, n))
        plot(time, trM(:, n))
        plot(time, trS(:, n))
    end
end