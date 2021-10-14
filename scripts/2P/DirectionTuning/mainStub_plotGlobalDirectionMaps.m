
% load data: plotData.m

% Determine average preferred direction and variance across nearby neurons
valid = ~isnan(data.RFX) & ~isnan(data.RFY) & ~isnan(data.direction) & ...
    data.sigdir==1;
[X, Y, Pref, Var, N] = gratings.getSmoothMap(data.RFX(valid), ...
    data.RFY(valid), data.direction(valid));

% Plot map of direction preferences (colour); transparency reflects number
% of neurons per spot
gratings.plotSmoothMap(X, Y, Pref, N, 'pref', 20);

% Plot map of direction selectivity (colour); transparency reflects number
% of neurons per spot
gratings.plotSmoothMap(X, Y, Var, N, 'var', 20);

% Plot map of direction preferences (color hue), gray scale reflects
% selectivity; patches with fewer than 5 neurons are completely transparent
gratings.plotSmoothMap_DirSel(X, Y, Pref, Var, N, 5);