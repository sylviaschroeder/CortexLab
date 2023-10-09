function rf_visDeg = interpolateRFtoVisualDegrees(rf, stimPos)

[x0, y0] = meshgrid(1:size(rf,2), 1:size(rf,1));
% vectors x1 and y1 specify gridpoints with distance of 1
% degree (diff(stimPos(...))); values match position of
% gridpoints in pixels of stimulus row/column
x1 = linspace(0.5, size(rf,2)+0.5, diff(stimPos(1:2)));
y1 = linspace(0.5, size(rf,1)+0.5, -diff(stimPos(3:4)));
% need to delete pixel values outside  given stimulus pixels,
% so we can use interpolation (rather than extrapolation) when
% mapping the RF from stimulus pixels to visual degrees
x2 = x1;
x2(x1<1 | x1>size(rf,2)) = [];
y2 = y1;
y2(y1<1 | y1>size(rf,1)) = [];
[x2, y2] = meshgrid(x2, y2);
rf_visDeg = interp2(x0, y0, rf, x2, y2);