function parameters = getModelParameters(model)

parameters = [];

predNames = model.PredictorNames;
coefNames = model.CoefficientNames;
coefs = model.Coefficients.Estimate;

% get stimulus parameters (intercepts)
intercept = strcmp('(Intercept)', coefNames);
stimInds = strncmp('stim_', coefNames, 5);
if sum([intercept stimInds]) == 0
    stimP = [];
elseif sum(stimInds) == 0
    stimP = coefs(intercept);
else
    stimP = ones(1, sum(stimInds)+1) * coefs(intercept) + ...
        [0, coefs(stimInds)'];
end
% get nonvis parameters
first = strcmp(predNames{1}, coefNames);
otherInds = strncmp([predNames{1} ':stim_'], coefNames, 12);
% first = strcmp('nonvis', coefNames);
% otherInds = strncmp('nonvis:stim_', coefNames, 12);
if sum([first otherInds]) == 0
    nonvisTotal = [];
    nonvisEachStim = [];
elseif sum(otherInds) == 0
    nonvisTotal = coefs(first);
    nonvisEachStim = [];
else
    tmp = ones(1, sum(otherInds)+1) * coefs(first) + ...
        [0, coefs(otherInds)'];
    nonvisTotal = mean(tmp);
    nonvisEachStim = tmp - mean(tmp);
end
% get baseline parameters
first = strcmp(predNames{2}, coefNames);
otherInds = strncmp([predNames{2} ':stim_'], coefNames, 10);
% first = strcmp('base', coefNames);
% otherInds = strncmp('base:stim_', coefNames, 10);
if sum([first otherInds]) == 0
    baseTotal = [];
    baseEachStim = [];
elseif sum(otherInds) == 0
    baseTotal = coefs(first);
    baseEachStim = [];
else
    tmp = ones(1, sum(otherInds)+1) * coefs(first) + ...
        [0, coefs(otherInds)'];
    baseTotal = mean(tmp);
    baseEachStim = tmp - mean(tmp);
end
% get interaction parameter (between baseline and nonvisual)
% ind = strcmp('nonvis:base', coefNames);
% if sum(ind) == 0
%     interactP = [];
% else
%     interactP = coefs(ind);
% end

parameters.stimuli = stimP;
parameters.([predNames{1} 'Total']) = nonvisTotal;
parameters.([predNames{1} 'EachStim']) = nonvisEachStim;
parameters.([predNames{2} 'Total']) = baseTotal;
parameters.([predNames{2} 'EachStim']) = baseEachStim;
% parameters.interactPars = interactP;