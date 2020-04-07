function error = stimOnlySingleCell(respTraining, varsTraining, ...
    respTest, varsTest)
% respTraining  [trial x stimulus]; contains 1 NaN
% varsTraining  [trial x stimulus]; contains 1 NaN
% respTest      [trial x stimulus]; all NaN, except one entry
% varsTest      [trial x stimulus]; all NaN, except one entry

[trial,stim] = find(~isnan(respTest));
prediction = nanmean(respTraining(:,stim));

error = respTest(trial,stim) - prediction;