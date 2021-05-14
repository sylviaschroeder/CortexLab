function data = getRunningData(folder)

data.running = readNPY(fullfile(folder, '_ss_running.speed.npy'));
data.time = readNPY(fullfile(folder, '_ss_running.timestamps.npy'));