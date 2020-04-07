%% Define data set
subject = 'M150410_SS045';
expDate = '2015-05-04';
exp = 1;
plane = 1;
channel = '_channel001';
rect = '_rect01_01';

info = ppbox.infoPopulate(subject, expDate, exp);

%% Define parameters
neuropilInnerRing = 3; % in um; distance between neuron ROI and inner ring 
                       % of  neuropil ROI
neuropilOuterRing = 20; % in um; distance betwen neuron ROI and outer 
                        % ring of neuropil ROI

%% Get necessary data
filePath = fullfile(info.folderProcessed, sprintf('%s_plane%03d%s%s_registered', ...
    info.basename2p, plane, channel, rect));
frames = loadArr(filePath);
% load meta
[sz, ~, infoROI] = loadArrInfo(strrep(filePath, 'registered', 'ROI'));

[sizeHoriz, sizeVert] = ssLocal.getFieldOfViewSize(info.zoomFactor);
pxPerMicronHoriz = size(infoROI.targetFrame,2) / sizeHoriz;
pxPerMicronVert = size(infoROI.targetFrame,1) / sizeVert;

%% Get neuropil masks
addpath('\\ZSERVER\Code\Neuropil Correction');

neurons = strcmp(infoROI.ROI.CellClasses, 's');
cellMasks = false(sum(neurons), sz(1), sz(2));
allMasks = zeros(sz([1 2]));
j = 1;
for k = 1:length(infoROI.ROI.CellMaps)
    if neurons(k)
        mask = false(sz([1 2]));
        mask(infoROI.ROI.CellMaps{k}) = true;
        cellMasks(j,:,:) = mask;
        j = j + 1;
    end
    allMasks(infoROI.ROI.CellMaps{k}) = 1;
end
options.inNeurop = neuropilInnerRing;
options.outNeurop = neuropilOuterRing;
neuropilMasks = createNeuropilMasks(cellMasks, allMasks, pxPerMicronHoriz, ...
    pxPerMicronVert, options);

%% Estiamted corrected neural signals
neuropilTraces = reshape(neuropilMasks, size(neuropilMasks,1), ...
    prod(sz([1 2]))) * double(reshape(frames, prod(sz([1 2])), sz(3)));
neuropilTraces = bsxfun(@rdivide, neuropilTraces, sum(sum(neuropilMasks, 2), 3));
[F_corrected, parameters] = estimateNeuropil(infoROI.F(:,neurons)', neuropilTraces);