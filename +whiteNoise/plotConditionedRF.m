function plotConditionedRF(meta, neuron, condition, plotConditionData, plane)

RFtypes = {'ON', 'OFF', 'Absolute'};
method = 'corrs';

frameTimes = ppbox.getFrameTimes(meta);
stimTimes = ppbox.getStimTimes(meta);
[stimFrames, stimFrameTimes, stimPosition] = whiteNoise.getStimulusFrames(meta);

[negativeCond, positiveCond, conditionTime, labels] = ...
    ssLocal.classifyData(meta, condition, plotConditionData);

titleLine = sprintf('Responses to white noise of neuron %d in plane %d', ...
    neuron, plane);
ssLocal.compareRFbetweenConditions(meta.F(:, neuron), frameTimes, ...
    [negativeCond, positiveCond], conditionTime, labels, stimFrames, ...
    stimFrameTimes, stimTimes, stimPosition, RFtypes, method, titleLine)

end