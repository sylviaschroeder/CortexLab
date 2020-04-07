function [detectedFrames, errorMedian, error] = ...
    detectMovingFrames(movie, sampleRate, options)

if nargin<3 || ~isfield(options, 'nParThreads')
    options.nParThreads = 4;
end
if nargin<3 || ~isfield(options, 'nFramesPerChunk')
    nFrames = size(movie, 3);
    options.nFramesPerChunk = nFrames;
end
if nargin<3 || ~isfield(options, 'upperThreshold') % to detect deviations
    options.upperThreshold = 2; % in STDs
end
if nargin<3 || ~isfield(options, 'lowerThreshold') % to detect deviations
    options.lowerThreshold = 1; % in STDs
end
if nargin<3 || ~isfield(options, 'window')
    options.window = 2; % in sec
end
if nargin<3 || ~isfield(options, 'doPlot')
    options.doPlot = 1;
end

options.spFilter = 1;
options.resolution = 0;
medianFilterOrder = round(100*sampleRate);
medianWindow = round(5*sampleRate);
% smoothPar = 10^(-11);
% filter parameters (exact choice doesn't seem to matter much)
filterWindow = 3; % in s
filterOrder = 3;

[~,~,~, error] = img.regTranslationsNew(single(movie), ...
    single(mean(movie,3)), options);

% subtract straight line fit
errorDetrend = detrend(error);

% fit very smooth spline to data and subtract it
% splineFit = fit((1:length(error))', error', ...
%     'smoothingspline', 'SmoothingParam', smoothPar);
% errorSpline = error - splineFit(1:length(error))';

% fit median filter and subtract
ind = round(medianFilterOrder/2);
filtered = medfilt1([ones(1,ind)*median(errorDetrend(1:medianWindow)), errorDetrend, ...
    ones(1,ind)*median(errorDetrend(end-medianWindow+1:end))], ...
    medianFilterOrder);
errorMedian = errorDetrend - filtered(ind+1:end-ind);

% apply smoothing filter to data (Savitzky-Golay FIR filter)
errorFilt = sgolayfilt(errorMedian, filterOrder, ...
    round((filterWindow * sampleRate) / 2) * 2 + 1);

% detect large deviations using a Schmitt-Trigger
errorStd = std(errorFilt);
detectedFrames = general.SchmittTrigger(abs(errorFilt), ...
    options.upperThreshold * errorStd, options.lowerThreshold * errorStd);

% add time before and after detected frames (as detecting may be late)
detectedFrames = conv(double(detectedFrames), ...
    ones(1, ceil(options.window * sampleRate)), 'same');
detectedFrames = detectedFrames >= 1;

if options.doPlot == 1
    figure('Position', [180 680 1720 420])
    subplot(1,3,1)
    plot(1:length(error), error)
    xlim([1 length(error)])
    title('Difference to average frame')
    xlabel('# Frame')
    
    subplot(1,3,2)
    plot(1:length(errorDetrend), errorDetrend, 'c')
    hold on
    plot(1:length(errorMedian), errorMedian, 'k')
    xlim([1 length(error)])
    xlabel('# Frame')
    title('Detrended error and median filter')
    
    subplot(1,3,3)
    h = [0 0];
    plot(1:length(error), errorDetrend-errorMedian, 'k');
    h(1) = plot(1:length(error), errorFilt, 'b');
    hold on
    plot([1 length(error)], [1 1]*options.upperThreshold * errorStd, 'b--')
    plot([1 length(error)], [1 1]*options.lowerThreshold * errorStd, 'b:')
    h(2) = plot(1:length(error), double(detectedFrames) * ...
        std(errorFilt), 'r');
    xlim([1 length(error)])
    legend(h, 'error', 'detected')
    xlabel('# Frame')
    title('Preprocessed error and detected frames')
end