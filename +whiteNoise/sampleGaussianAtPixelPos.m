function fitRF = sampleGaussianAtPixelPos(rf, stimPos, fitPars)

[x0, y0] = meshgrid(1:size(rf,2), 1:size(rf,1));
% vectors x1 and y1 specify gridpoints with distance of 1
% degree (diff(stimPos(...))); values match position of
% gridpoints in pixels of stimulus row/column; all
% locations within edges of stimulus are included (no cropping)
x1 = linspace(0.5, size(rf,2)+0.5, diff(stimPos(1:2)));
y1 = linspace(0.5, size(rf,1)+0.5, -diff(stimPos(3:4)));
[x1, y1] = meshgrid(x1, y1);
% fitted Gaussian was based on cropped, upsampled RF (using
% size of x2 (y2)); need to add points beyond edges of x2 (y2)
% to match grid size and relative location of x1 (y1)
[x3, y3] = meshgrid((1:size(x1,2)) - sum(x1(1,:)<1), ...
    (1:size(y1,1)) - sum(y1(:,1)<1));
fitRF = whiteNoise.D2GaussFunctionRot(fitPars, cat(3, x3, y3));
% sample fitted RF at positions measured for original RF map
fitRF = interp2(x1, y1, fitRF, x0, y0);