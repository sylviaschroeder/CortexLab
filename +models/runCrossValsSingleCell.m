function crossVals = runCrossValsSingleCell(responses, nonVisual, ...
    baselines)
% responses     [trial x stimulus]
% nonVisual     [trial x stimulus]
% baselines     [trial x stimulus]

% Models including 3 variables: stimulus, nonvisual signal and baseline
formulas = {'resp ~ 1 + stim*(nonvis+base)', ...
    'resp ~ 1 + stim*nonvis + base', ...
    'resp ~ 1 + stim*base + nonvis', ...
    'resp ~ 1 + stim + nonvis + base', ...
    'resp ~ 1 + stim*base', ...
    'resp ~ 1 + stim*nonvis', ...
    'resp ~ 1 + stim + base', ...
    'resp ~ 1 + stim + nonvis', ...
    'resp ~ 1 + nonvis + base', ...
    'resp ~ 1 + nonvis', ...
    'resp ~ 1 + base', ...
    'resp ~ 1 + stim'};
crossVals.formulas = formulas;
explVars = NaN(1, length(formulas));
for m = 1:length(formulas)
    errors = models.crossvalidate(@models.testFitlmModel, {responses}, ...
        {nonVisual, baselines}, formulas{m});
    explVars(m) = 1 - sum(errors{1}(:).^2) / ...
        sum((responses(:)-mean(responses(:))).^2);
end
crossVals.explVars = explVars;