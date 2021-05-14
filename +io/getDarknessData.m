function data = getDarknessData(folder)

data.interval = [];
if isfile(fullfile(folder, '_ss_recordings.darkness_intervals.npy'))
    data.interval = readNPY(fullfile(folder, ...
        '_ss_recordings.darkness_intervals.npy'));
end