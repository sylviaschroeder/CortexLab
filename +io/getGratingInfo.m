function data = getGratingInfo(folder)

if isfile(fullfile(folder, '_ss_grating._ss_gratingID.npy'))
    data.interval = readNPY(fullfile(folder, '_ss_grating._ss_gratingID.npy'));
else
    data = [];
    return
end
data.IDs = readNPY(fullfile(folder, '_ss_grating._ss_gratingID.npy'));
data.times = readNPY(fullfile(folder, '_ss_grating.intervals.npy'));
data.directions = readNPY(fullfile(folder, '_ss_gratingID.directions.npy'));