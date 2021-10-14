function plotSmoothMap_DirSel(X, Y, Pref, Sel, N, minN)
% Pref represents the preferred direction of each patch in
% the map, i.e. values are expected to by between 0 and 360 degrees
% Sel represents direction selectivity of each path in the
% map, i.e. values are expected to be between 0 and 1

% colormap for preferred directions
cm = hsv(n_cMap);
% transform preferred directions to indices into colors
im = cm(round(mod(Pref(:), 360)), :);
% represent direction selectivity by "mixing" colour for preferred
% direction with gray value
im = Sel(:) .* im + (1 - Sel(:)) .* [0.5 0.5 0.5];
im = reshape(im, size(Pref,1), size(Pref,2), 3);
% set pixels with too few data points (N < minN) to transparent
aVal = double(N >= minN);

figure
image(X([1 end]), Y([1 end]), im, 'AlphaData', aVal)
