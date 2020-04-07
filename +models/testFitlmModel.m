function error = testFitlmModel(respTraining, varsTraining, respTest, ...
    varsTest, model)
% respTraining  [trial x stimulus]; contains 1 NaN
% varsTraining  {1 x 2}, each entry: [trial x stimulus]; contains 1 NaN
% respTest      [trial x stimulus]; all NaN, except one entry
% varsTest      {1 x 2}, each entry: [trial x stimulus]; all NaN, except one entry

[trial, stimulus] = find(isnan(respTraining));
ind = find(isnan(respTraining));
resp = respTraining(:);
resp(ind) = [];
nonvis = varsTraining{1}(:);
nonvis(ind) = [];
base = varsTraining{2}(:);
base(ind) = [];
stim = reshape(repmat(1:size(respTraining,2), size(respTraining,1), 1), [], 1);
stim(ind) = [];

tbl = table(nonvis, base, stim, resp);
mdl = fitlm(tbl, model, 'CategoricalVars', 3);
testData = [varsTest{1}(trial,stimulus), varsTest{2}(trial,stimulus), stimulus];
error = predict(mdl, testData) - respTest(trial,stimulus);