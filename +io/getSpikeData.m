function data = getSpikeData(folder, vars)

data.times = readNPY(fullfile(folder, 'spikes.times.npy'));
data.clusters = readNPY(fullfile(folder, 'spikes.clusters.npy'));
data.depths = readNPY(fullfile(folder, 'spikes.depths.npy'));
data.amps = readNPY(fullfile(folder, 'spikes.amps.npy'));
data.templates = readNPY(fullfile(folder, 'spikes.templates.npy'));
data.IDs = readNPY(fullfile(folder, 'clusters.ids.npy'));
data.groups = readNPY(fullfile(folder, 'clusters.groups.npy'));

if nargin < 2
    return
end

if strcmp('waveforms', vars)
    data.waveforms = readNPY(fullfile(folder, 'templates.waveforms.npy'));
    data.waveformChans = readNPY(fullfile(folder, 'templates.waveformsChannels.npy'));
end
if strcmp('coordinates', vars)
    data.chanCoordinates = readNPY(fullfile(folder, 'channels.localCoordinates.npy'));
end