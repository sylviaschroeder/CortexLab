function [dx, dy] = alignRectangles( ...
    targetAverage, targetPos, targetValidX, targetValidY, ...
    otherAverage, otherPos, otherValidX, otherValidY, options)

if nargin<9 || ~isfield(options, 'resolution')
    options.resolution = 2; % [1/pix],
    % e.g. if options.resolution = 10 the registration resolution will be
    % down to 1/10 pixels
end

regOptions.nParThreads = 1;
regOptions.resolution = options.resolution;
regOptions.spFilter = 1;

targetIm = zeros(targetPos.ymax-targetPos.ymin+1, ...
    targetPos.xmax-targetPos.xmin+1);
targetIm(targetValidY-targetPos.ymin+1, targetValidX-targetPos.xmin+1) = ...
    targetAverage;
otherIm = zeros(otherPos.ymax-otherPos.ymin+1, ...
    otherPos.xmax-otherPos.xmin+1);
otherIm(otherValidY-otherPos.ymin+1, otherValidX-otherPos.xmin+1) = ...
    otherAverage;

minx = max(targetValidX(1) - targetPos.xmin + 1, ...
    otherValidX(1) - otherPos.xmin + 1);
maxx = min(targetValidX(end) - targetPos.xmin + 1, ...
    otherValidX(end) - otherPos.xmin + 1);
miny = max(targetValidY(1) - targetPos.ymin + 1, ...
    otherValidY(1) - otherPos.ymin + 1);
maxy = min(targetValidY(end) - targetPos.ymin + 1, ...
    otherValidY(end) - otherPos.ymin + 1);

targetIm = targetIm(miny : maxy, minx : maxx);
otherIm = otherIm(miny : maxy, minx : maxx);

[dx, dy] = img.regTranslationsNew(single(otherIm), single(targetIm), regOptions);
fprintf('\n')