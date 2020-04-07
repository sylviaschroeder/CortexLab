subject = 'M141002_SS026';
expDate = '2014-10-29';
exp = 4;
info=ppbox.infoPopulate(subject, expDate, exp);
iPlane = 2;
neuron = 7;

filePath = fullfile(info.folderProcessed, sprintf('%s_plane%03d_ROI', ...
    info.basename2p, iPlane));
load(filePath, 'meta')

[~, ~, stimulusMatrix, frameTimes] = ssLocal.getStimulusResponseInfo(meta);
timeFrames = frameTimes;
calciumTrace = meta.F(:,neuron)';

ballData = nonVis.getRunningSpeed(meta);
velocityFrames = interp1(ballData.t, ballData.total, frameTimes, 'pchip');