subject='M150410_SS044';
expDate='2015-04-28';
exp=4;
plane = 4;
neuron = 35;

% for removal of slow drift in neural responses (high-pass filter)
window = 400; %in sec
percentile = 5; 

info=ppbox.infoPopulate(subject, expDate, exp);
basenameROI = sprintf('%s_plane%03d_ROI', info.basename2p, plane);
load(fullfile(info.folderProcessed, basenameROI), 'meta');

frameTimes = ppbox.getFrameTimes(meta);
stimTimes = ppbox.getStimTimes(meta);
[stimFrames, stimFrameTimes, stimPosition] = whiteNoise.getStimulusFrames(meta);

RFtypes = {'Absolute', 'ON', 'OFF'};
colors = 'rb';

% Preprocess calcium trace
trace = meta.F(:,neuron);
[~,F_0] = ssLocal.removeSlowDrift(trace, frameTimes, window, percentile);
trace = (trace - F_0) ./ F_0;

%% Use Pearson correlation or kernel method to calculate spatiotemporal RF
% use 'corrs' or 'kernels' as last argument
[receptiveFields, RFtimes] = whiteNoise.getReceptiveField( ...
    trace, frameTimes, stimFrames, stimFrameTimes, ...
    stimTimes, RFtypes, 1, 'corrs');

%% Plot  spatiotemporal RF
smoothedRFs = whiteNoise.plotReceptiveField(receptiveFields, RFtimes, ...
    stimPosition, RFtypes);

figure
RFinds = [];
for type = 1:length(RFtypes)
    if max(abs(smoothedRFs{type}(:))) / std(smoothedRFs{type}(:)) > 5
        infRF(type) = whiteNoise.fitGaussianToRF(smoothedRFs{type}, ...
            RFtimes, stimPosition, 0, RFtypes{type}, 0);
        RFinds = [RFinds, type];
    end
end
infRF = infRF(RFinds);
whiteNoise.plotFitReceptiveField([size(smoothedRFs{1},1) ...
    size(smoothedRFs{1},2)], ...
    {infRF.gaussianFit}, stimPosition, RFtypes(RFinds), [infRF.RFsign]);

% figure
% whiteNoise.plotFitReceptiveField([size(receptiveFields{1},1), ...
%     size(receptiveFields{1},1)], {RFinfo.gaussianFit}, stimPosition, ...
%     RFtypes, [RFinfo.RFsign]);