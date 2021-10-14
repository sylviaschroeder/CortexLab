function [X, Y, Pref, Var, N] = getSmoothMap(rf_x, rf_y, prefDir, dirSelectivity)

% Parameters
gridBin = 2; % visual degrees, distance between values to be calculated for map
radius = 10; % visual degrees, radius for including nearest neighbours for 
             % each grid point used to calculate map values

% Get grid (x and y)
x_range = [floor(min(rf_x) / gridBin) * gridBin, ceil(max(rf_x) / gridBin) * gridBin];
y_range = [floor(min(rf_y) / gridBin) * gridBin, ceil(max(rf_y) / gridBin) * gridBin];
[X, Y] = meshgrid(x_range(1):gridBin:x_range(2), y_range(1):gridBin:y_range(2));

% Transform preferred directions and direction selectivity into vector
% format
[dir_x, dir_y] = pol2cart(deg2rad(prefDir), dirSelectivity);

% Calculate average preference, variance and number of data points for each
% grid point
Pref = NaN(size(X));
Var = NaN(size(X));
N = NaN(size(X));
for k = 1:numel(X)
    % distance of grid point to RF of each neuron
    d = sqrt((X(k) - rf_x).^2 + (Y(k) - rf_y).^2);
    % find neighbours within radius
    neighbours = d <= radius;
    % polar coordainates of vector average of direction preferences
    [theta, rho] = cart2pol(mean(dir_x(neighbours)), mean(dir_y(neighbours)));
    Pref(k) = theta;
    Var(k) = rho;
    N(k) = sum(neighbours);
end