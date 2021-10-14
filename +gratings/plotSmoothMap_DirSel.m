function plotSmoothMap_DirSel(X, Y, Pref, Sel, N, minN)
% Pref represents the preferred direction of each patch in
% the map, i.e. values are expected to by between 0 and 360 degrees
% Sel represents direction selectivity of each path in the
% map, i.e. values are expected to be between 0 and 1

% colormap for preferred directions
cm = hsv(360);
% transform preferred directions to indices into colors
im = round(mod(Pref(:)-1, 360)) + 1;
im(isnan(im)) = 1;
im = cm(im, :);
% represent direction selectivity by "mixing" colour for preferred
% direction with gray value
im = Sel(:) .* im + (1 - Sel(:)) .* [0.5 0.5 0.5];
im = reshape(im, size(Pref,1), size(Pref,2), 3);
% set pixels with too few data points (N < minN) to transparent
aVal = double(N >= minN);

figure
image(X([1 end]), Y([1 end]), im, 'AlphaData', aVal)
xlabel('Azimuth')
ylabel('Elevation')

% colormatrix
figure
c = repmat(permute(cm, [3 1 2]), 20, 1, 1);
gray = linspace(0, 1, 20)';
binSz = diff(gray(1:2));
c = gray .* c + (1 - gray) .* 0.5;
image([0.5 359.5], [0+0.5*binSz, 1-0.5*binSz], c)
set(gca, 'XTick', 0:90:360, 'YTick', 0:0.25:1)
xlabel('Direction (visual degrees)')
ylabel('Selectivity')
title('Colormap')