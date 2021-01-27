function data = getTimelineData(folder, vars)

data.wheelMoves = readNPY(fullfile(folder, '_ss_wheelMoves.intervals.npy'));
data.wheelDisplacements = readNPY(fullfile(folder, '_ss_wheelMoves.displacement.npy'));
data.wheelPeakTimes = readNPY(fullfile(folder, '_ss_wheelMoves.peak_times.npy'));
data.wheelAmplitudes = readNPY(fullfile(folder, '_ss_wheelMoves.amplitude.npy'));

if nargin < 2
    return
end

if any(strcmp({'lickSignal','wheelPos', 'wheelVel'}, vars))
    data.t_TL = readNPY(fullfile(folder, '_ss_signals.timestamps.npy'));
end
if strcmp('lickSignal', vars)
    data.lickSignal = readNPY(fullfile(folder, '_ss_lickPiezo.raw.npy'));
end
if strcmp('wheelPos', vars)
    data.wheelPos = readNPY(fullfile(folder, '_ss_wheel.position.npy'));
end
if strcmp('wheelVel', vars)
    data.wheelVel = readNPY(fullfile(folder, '_ss_wheel.velocity.npy'));
end
