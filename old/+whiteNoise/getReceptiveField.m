function [receptiveField, frameTimes] = getReceptiveField(trace, ...
    traceTimes, stimFrames, stimFrameTimes, repetitionTimes, RFtype, ...
    plotResponseTraces)
%GETRECEPTIVEFIELD Returns spatiotemporal receptive field.
%   [receptiveField, frameTimes] = GETRECEPTIVEFIELD(trace, ...
%    traceTimes, stimFrames, stimFrameTimes, repetitionTimes, RFtype, ...
%    plotResponseTraces) correlates the neural response with the noise
%    stimulus.
%
%   receptiveField      {1 x RFtypes}; in each entry: [rows x cols x time]
%                       containing Pearson's correlation between the pixel
%                       of the stimulus with the response at the specified 
%                       time lag 
%   frameTimes          in sec
%
%   trace               [n x 1]; calcium trace of neuron
%   traceTimes          [n x 1]; sample times of calcium trace
%   stimFrames          [rows x cols x m]; noise stimulus
%   stimFrameTimes      [1 x m]; times of stimulus frames
%   repetitionTimes     struct; contains stimulus on- and offsets of all
%                       stimulus repetitions
%   RFtype              'ON', 'OFF', 'Linear', or 'Absolute'
%   plotResponseTraces  0 or 1

RFlimits = [0 1.0]; % in sec

stimTimeStep = median(diff(stimFrameTimes));
RFtimesInFrames = floor(RFlimits(1) / stimTimeStep) : ...
    ceil(RFlimits(2) / stimTimeStep);
lagInFrames = max(abs(RFtimesInFrames));
corrTimesInFrames = -lagInFrames : lagInFrames;
[~, RFFrameInds] = intersect(corrTimesInFrames, RFtimesInFrames);

if plotResponseTraces == 1
    screenSize = get(0, 'ScreenSize');
    figure('Position', [10 695 screenSize(3)-20 420])
    hold on
end

stim = reshape(stimFrames, [], size(stimFrames, 3));

trialTraces = NaN(length(repetitionTimes.onset), length(stimFrameTimes));
for rep = 1:length(repetitionTimes.onset)
    % interpolate trace to times of stimulus frames
    medianTrace = interp1(traceTimes, trace, ...
        repetitionTimes.onset(rep) + stimFrameTimes);
    if isnan(medianTrace(1)) && rep == 1 % record of neural response started later than stimulus
        medianTrace(1) = trace(1);
    end
    if isnan(medianTrace(end)) && rep == length(repetitionTimes.onset) % record of neural response finished earlier than stimulus
        medianTrace(end) = trace(end);
    end
    trialTraces(rep, :) = medianTrace;
    
    if plotResponseTraces == 1
        plot(stimFrameTimes, medianTrace, 'Color', ...
            ones(1, 3) * (rep - 1) / length(repetitionTimes.onset))
    end
end
medianTrace = nanmedian(trialTraces, 1);
if plotResponseTraces == 1
    plot(stimFrameTimes, medianTrace, 'r', 'LineWidth', 2)
    axis tight
    xlabel('Time (in s)')
    ylabel('Raw Ca-response')
end

% z-score neural response
medianTrace = (medianTrace - mean(medianTrace)) ./ std(medianTrace);

% calculate RFs of all specified types
if ischar(RFtype)
    RFtype = {RFtype};
end
receptiveField = cell(1, length(RFtype));
for type = 1:length(RFtype)
    st = stim;
    % adapt stimulus to match RF type
    switch RFtype{type}
        case 'ON'
            st(st < 0) = 0;
        case 'OFF'
            st(st > 0) = 0;
        case 'Linear' % don't do anything
        case 'Absolute'
            st = abs(st);
        otherwise
            disp('RFtype must be ''ON'', ''OFF'', ''Linear'', or ''Absolute''!')
            frameTimes = [];
            return
    end
    st = (st - mean(st(:))) ./ std(st(:));
    
    RF = zeros(size(stimFrames, 1) * size(stimFrames, 2), ...
        length(RFtimesInFrames));
    for pix = 1:size(stim, 1)
        c = xcorr(medianTrace, st(pix,:), lagInFrames, 'unbiased');
        RF(pix, :) = c(RFFrameInds);
    end
    receptiveField{type} = reshape(RF, size(stimFrames, 1), ...
        size(stimFrames, 2), []);
end
frameTimes = corrTimesInFrames(RFFrameInds) * stimTimeStep;