function error = threeVarsBaselineConstrained(respTraining, varsTraining, ...
    respTest, varsTest, bestModel)

error = NaN;

% in case a whole stimulus (e.g. blank) was excluded by setting values to
% NaN, don't detect those as the left out stimulus in crossvalidation
resp = respTraining;
nanStim = all(isnan(resp),1);
resp(:,nanStim) = 0;
[trial,stim] = find(isnan(resp));
rTrain = num2cell(respTraining, 1);
rTrain{stim}(trial) = [];
if sum(nanStim) > 0
    rTrain{nanStim} = [];
end
nvTrain = num2cell(varsTraining{1}, 1);
nvTrain{stim}(trial) = [];
blTrain = num2cell(varsTraining{2}, 1);
blTrain{stim}(trial) = [];

coefNames = bestModel.CoefficientNames;
coefs = bestModel.Coefficients.Estimate;
numStim = size(respTraining,2);
terms = [0 0]; % [nonvis, interact]
% inital parameters
x0 = cell(1, 4);

% get stimulus parameters (intercepts)
intercept = strcmp('(Intercept)', coefNames);
stimInds = strncmp('stim_', coefNames, 5);
if sum(stimInds) == 0
    return
else
    x0{1} = ones(numStim,1) * coefs(intercept) + [0; coefs(stimInds)];
end
% get nonvis parameters
first = strcmp('nonvis', coefNames);
otherInds = strncmp('nonvis:stim_', coefNames, 12);
if sum(otherInds) > 0
    terms(1) = 1;
    x0{2} = ones(numStim,1) * coefs(first) + [0; coefs(otherInds)];
end
% get baseline parameters
first = strcmp('base', coefNames);
otherInds = strncmp('base:stim_', coefNames, 10);
if sum(otherInds) == 0
    return
else
    x0{3} = ones(numStim,1) * coefs(first) + [0; coefs(otherInds)];
end
% get interaction parameter (between baseline and nonvisual)
ind = strcmp('nonvis:base', coefNames);
if sum(ind) > 0
    terms(2) = 1;
    x0{4} = coefs(ind);
end

% initial estimtate of linear dependence of baseline on intercepts
b =  [ones(length(x0{1}),1), x0{1}] \ x0{3};
x0{3} = b;
x0 = cat(1, x0{:});
options = optimoptions('fminunc');
options.MaxFunEvals = 500 * length(x0);
options.Display = 'off';
options.Algorithm = 'quasi-newton';
x = fminunc(@(x)myminfun(rTrain, nvTrain, blTrain, x, terms), x0, options);
if all(terms == 1) % full model
    x1 = x([stim, numStim+stim, end-2, end-1, end]);
elseif terms(1) == 1 % no interaction
    x1 = [x([stim, numStim+stim, end-1, end]); 0];
else % no nonvis
    x1 = [x(stim); 0;x([end-1, end]); 0];
end
error = mymodel(varsTest{1}(trial,stim), varsTest{2}(trial,stim), x1) - ...
    respTest(trial,stim);


function error = myminfun(resp, nonvis, baseline, parameters, terms)

errors = cell(size(resp));
numStim = length(resp);
for iStim = 1:numStim
    if isempty(resp{iStim})
        errors{iStim} = 0;
        continue
    end
    if all(terms == 1) % full model
        pars = parameters([iStim, numStim+iStim, end-2, end-1, end]);
    elseif terms(1) == 1 % no interaction
        pars = [parameters([iStim, numStim+iStim, end-1, end]); 0];
    else % no nonvis
        pars = [parameters(iStim); 0; parameters([end-1, end]); 0];
    end
    errors{iStim} = mymodel(nonvis{iStim}, baseline{iStim}, pars) - ...
        resp{iStim};
end
errors = cat(1, errors{:});
error = sqrt(sum(errors .^ 2));

function estimate = mymodel(nonvis, baseline, parameters)

estimate = parameters(1) + parameters(2) * nonvis + ...
    (parameters(3) + parameters(4) * parameters(1)) * baseline + ...
    parameters(5) * nonvis .* baseline;