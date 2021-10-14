function plotSmoothMap(X, Y, Z, N, type, maxN)
% if type == 'pref', Z represents the preferred direction of each patch in
% the map, i.e. values are expected to by between 0 and 360 degrees
% if type ~= 'pref', Z represents direction selectivity of each path in the
% map, i.e. values are expected to be between 0 and 1

% Parameters
n_cMap = 100;

if nargin < 6
    maxN = max(N(:));
end

if strcmp(type, 'pref')
    cm = hsv(n_cMap);
    zNorm = round(mod(Z, 360) ./ 360 .* n_cMap);
else
    cm = winter(n_cMap);
    zNorm = round(Z .* n_cMap);
end


figure
image(X([1 end]), Y([1 end]), zNorm, 'AlphaData', N, 'AlphaDataMapping', 'scaled')
colormap(cm)
ax = gca;
ax.ALim = [0 maxN];