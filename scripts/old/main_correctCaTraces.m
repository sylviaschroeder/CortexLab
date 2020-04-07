%% Define data set
subject = 'M141007_SS029';
expDate = '2014-11-12';
exp = 2;
plane = 2;
info = ppbox.infoPopulate(subject, expDate, exp);

%% Define parameters
pixelSize = 1; % 1 for zoom 3, 0.68 for zoom 2
neuropilInnerRing = 3; % in um; distance between neuron ROI and inner ring 
                       % of  neuropil ROI
neuropilOuterRing = 20; % in um; distance betwen neuron ROI and outer 
                        % ring of neuropil ROI
percentileBaseline = 8; % as in Dombeck and Tank (2014)
movingWindowBaseline = 300; % in sec; size of moving window, in which 
                            % baseline is determined

%% Get necessary data
filePath = fullfile(info.folderProcessed, sprintf('%s_plane%03d_registered', ...
    info.basename2p, plane));
frames = loadArr(filePath);
% load meta
load(fullfile(info.folderProcessed, ...
    sprintf('%s_plane%03d_ROI', info.basename2p, plane)), 'meta')

frameTimes = ppbox.getFrameTimes(meta);
samplingRate = 1 / median(diff(frameTimes));

%% Determine correction factor for neuropil (usually 0.7)
neuropilRatios = ssLocal.extractNeuropilRatio(frames, pixelSize, ...
    neuropilInnerRing, neuropilOuterRing);
display(['Ratios of F_vessel and F_surround: ' num2str(mean(neuropilRatios, 1))])

neuropilFactor = 0.7;

%% Calculate F_corrected
neuropilTraces = ssLocal.extractNeuropil(frames, meta.ROI.CellMaps, ...
    pixelSize, neuropilInnerRing, neuropilOuterRing);
deltaF = meta.F - neuropilFactor * neuropilTraces;
baselines = ssLocal.extractBaselines(deltaF, samplingRate, ...
    percentileBaseline, movingWindowBaseline);
deltaF = (deltaF - baselines) ./ baselines;

%% make new meta structure and save
meta.deltaF = deltaF;
meta.pixelSize = pixelSize;
meta.neuropilFactor = neuropilFactor;
meta.neuropilInnerRing = 3;
meta.neuropilOuterRing = 20;
meta.percentileBaseline = 8;
meta.movingWindowBaseline = 300;

save(fullfile(info.folderProcessed, ...
    sprintf('%s_plane%03d_F_corrected', info.basename2p, plane)), 'meta')