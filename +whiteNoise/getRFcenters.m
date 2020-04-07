function RFcenters = getRFcenters(meta)

RFtype = 'Absolute';

frameTimes = ppbox.getFrameTimes(meta);
stimTimes = ppbox.getStimTimes(meta);
[stimFrames, stimFrameTimes, stimPosition] = whiteNoise.getStimulusFrames(meta);

RFcenters = NaN(length(meta.ROI.CellMaps), 2);

ROImaps = meta.ROI.CellMaps;
for neuron = 1:length(ROImaps)
    [receptiveField, RFtimes] = whiteNoise.getReceptiveField( ...
        meta.F(:, neuron), frameTimes, stimFrames, stimFrameTimes, ...
        stimTimes, RFtype, 0, 'corrs');
    smoothedRF = smooth3(receptiveField{1}, 'gaussian');
    if max(abs(smoothedRF(:))) / std(smoothedRF(:)) < 5
        continue
    end
    infoRF = whiteNoise.fitGaussianToRF(smoothedRF, ...
        RFtimes, stimPosition, 0, RFtype, 0);
    RFcenters(neuron,1) = stimPosition(1) + (infoRF.gaussianFit(2)-0.5) / ...
        size(smoothedRF,2) * diff(stimPosition(1:2));
    RFcenters(neuron,2) = stimPosition(3) + (infoRF.gaussianFit(4)-0.5) / ...
        size(smoothedRF,1) * diff(stimPosition(3:4));
end