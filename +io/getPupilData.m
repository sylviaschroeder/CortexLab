function data = getPupilData(folder)

data.pupilSize = readNPY(folder, 'eye.diameter.npy');
data.time = readNPY(folder, 'eye.timestamps.npy');