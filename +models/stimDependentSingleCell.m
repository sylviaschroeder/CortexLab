function error = stimDependentSingleCell(respTraining, varsTraining, ...
    respTest, varsTest)
% respTraining  [trial x stimulus]; contains 1 NaN
% varsTraining  [trial x stimulus]; contains 1 NaN
% respTest      [trial x stimulus]; all NaN, except one entry
% varsTest      [trial x stimulus]; all NaN, except one entry

[trial,stim] = find(isnan(respTraining));
rTrain = mat2cell(respTraining, size(respTraining,1), ...
    ones(1, size(respTraining,2)));
rTrain{stim}(trial) = [];
vTrain = mat2cell(varsTraining, size(varsTraining,1), ...
    ones(1, size(varsTraining,2)));
vTrain{stim}(trial) = [];

coeffs = general.bilinfit(vTrain, rTrain);

[trial,stim] = find(~isnan(respTest));
prediction = coeffs(stim,1) * varsTest(trial, stim) + coeffs(stim,2);
error = respTest(trial,stim) - prediction;