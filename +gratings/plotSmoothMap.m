function plotSmoothMap(X, Y, Z, N, type, maxN)
% if type == 'pref', Z represents the preferred direction of each patch in
% the map, i.e. values are expected to by between 0 and 360 degrees
% if type ~= 'pref', Z represents direction selectivity of each path in the
% map, i.e. values are expected to be between 0 and 1

% Parameters
n_cMap = 360;

if nargin < 6
    maxN = max(N(:));
end

if strcmp(type, 'pref')
    cm = hsv(n_cMap);
    zNorm = round(mod(Z, 360) ./ 360 .* n_cMap);
else
    cm = parula(n_cMap);
    zNorm = round(Z .* n_cMap);
end

figure
image(X([1 end]), Y([1 end]), zNorm, 'AlphaData', N, 'AlphaDataMapping', 'scaled')
colormap(cm)
ax = gca;
ax.ALim = [0 maxN];
if strcmp(type, 'pref')
    c = colorbar('Ticks', 0:90:360);
    c.Label.String = 'Preferred direction';
else
    c = colorbar('Ticks', 0:90:360, 'TickLabels', 0:0.25:1);
    c.Label.String = 'Direction selectivity';
end
xlabel('Azimuth')
ylabel('Elevation')