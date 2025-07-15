function [piezoPerVolume, planes_volume] = getPiezoMovement(ops)

piezoMicronsPerVolt = 400 / 5;

% get directory with piezo recording
metaDataPath = ops{'data_path'};
metaDataPath = string(metaDataPath{1});
metaDataPath = char(metaDataPath(1));
metaDataPath = ['Z' metaDataPath(2:end)];

% get number of planes
numPlanes = double(ops{'nplanes'});

% load data recorded with NiDAQ
f = dir(fullfile(metaDataPath, 'nidaqChannels*.csv'));
if isempty(f)
    return
end
nidaqChannels = readlines(fullfile(metaDataPath, f(1).name));
nidaqChannels = split(nidaqChannels(1), ',');
f = dir(fullfile(metaDataPath, 'NiDaqInput*.bin'));
if isempty(f)
    return
end
fid = fopen(fullfile(metaDataPath, f(1).name), 'r');
nidaq = fread(fid, 'float64');
fclose(fid);
  % discard last datapoints if not all channels have been recorded
nidaq = nidaq(1 : floor(length(nidaq) / length(nidaqChannels)) * ...
    length(nidaqChannels));
nidaq = reshape(nidaq, length(nidaqChannels), [])';
% get piezo data
piezo = nidaq(:, strcmp(nidaqChannels, 'piezo'));
t_nidaq = (1:size(nidaq,1))' ./ 1000;
% get times of 2p images
frameClock = nidaq(:, strcmp(nidaqChannels, 'frameclock'));
frameInds = [false; diff(frameClock > 4.0) > 0.5];
frameTimes = t_nidaq(frameInds);

% align piezo trace to frame times of each imaging plane
% start time of each imaged volume
volTimes = frameTimes(1 : numPlanes : end);
dt = mean(diff(volTimes));
[piezoPerVolume, t_volume] = traces.getAlignedTraces(piezo, t_nidaq, ...
    volTimes, [0 1.1*dt]);
piezoPerVolume = mean(piezoPerVolume(:,2:end), 2, 'omitnan');
planes_volume = t_volume ./ dt .* numPlanes;

% subtract highest point of piezo trace
piezoPerVolume = piezoPerVolume - min(piezoPerVolume);
% convert voltage of piezo into microns
piezoPerVolume = piezoPerVolume * piezoMicronsPerVolt;