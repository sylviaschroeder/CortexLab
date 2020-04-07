function [fitParameters, adjustedRsquared] = fitTuningCurves(tuningCurves, ...
    orientations)

% tuningCurves      [stimulus x neuron x condition]; contains
%                   response of each neuron to each stimulus in each
%                   condition; last condition: all data not divided into
%                   conditions
% orientations      [stimulus x 1]; contains orientation/direction of each
%                   stimulus

% fitParameters     [neuron x condition x parameter]; contains parameters
%                   resulting from fit for each neuron and each condition;
%                   parameters: preferred orientation (Dp), ampl. at pref.
%                   dir. (Rp), ampl. at opp. dir. (Rn), offset of tuning 
%                   curve from zero (Ro), tuning width (sigma)
% adjustedRsquared  [neuron x condition]; contains adjusted R_squared of
%                   fit for each neuron and each condition (plus all data 
%                   not divided into conditions)

numFitIterations = 10;

fitParameters = NaN(size(tuningCurves,2), size(tuningCurves,3), 5);

adjustedRsquared = NaN(size(tuningCurves,2), size(tuningCurves,3));

for iCell = 1:size(tuningCurves,2)
    % first fit data, which was divided into conditions, to fix preferred
    % orientation and tuning width
    tun = tuningCurves(:,iCell,end);
    ind = ~isnan(tun);
    pars = fitori(orientations(ind), tun(ind), [], [], '', numFitIterations);
    fitParameters(iCell,:,1) = pars(1); % preferred orientation
    fitParameters(iCell,:,5) = pars(5); % tuning width (sigma)
    fitParameters(iCell,3,2:4) = pars(2:4);
    modelFit = oritune(pars, orientations);
    n = sum(ind);
    Rsqu = 1 - sum((modelFit(ind) - tun(ind)).^2) / (n-6) / ...
        (sum((tun(ind) - mean(tun(ind))).^2) / (n-1));
    adjustedRsquared(iCell,end) = Rsqu;
    
    % fit trials of each condition, fixing preferred orientation and tuning
    % width
    for cond = 1:size(tuningCurves,3)-1
        tun = tuningCurves(:,iCell,cond);
        ind = ~isnan(tun);
        pars = fitori(orientations(ind), tun(ind), [], ...
            [fitParameters(iCell,cond,1), NaN, NaN, NaN, ...
            fitParameters(iCell,cond,5)], '', numFitIterations);
        fitParameters(iCell,cond,2:4) = pars(2:4);
        modelFit = oritune(pars, orientations);
        n = sum(ind);
        Rsqu = 1 - sum((modelFit(ind) - tun(ind)).^2) / (n-4) / ...
            (sum((tun(ind) - mean(tun(ind))).^2) / (n-1));
        adjustedRsquared(iCell,cond) = Rsqu;
    end
end