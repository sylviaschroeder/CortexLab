function error = stimSeparateSingleCell(respTraining, varsTraining, ...
    respTest, varsTest)
% respTraining  [trial x stimulus]; contains 1 NaN
% varsTraining  [trial x stimulus]; contains 1 NaN
% respTest      [trial x stimulus]; all NaN, except one entry
% varsTest      [trial x stimulus]; all NaN, except one entry

[trial,stim] = find(isnan(respTraining));
rTrain = respTraining(:,stim);
rTrain(trial) = [];
vTrain = varsTraining(:,stim);
vTrain(trial) = [];

x = [vTrain, ones(size(vTrain))] \ rTrain;

[trial,stim] = find(~isnan(respTest));
prediction = x(1) * varsTest(trial, stim) + x(2);
error = respTest(trial,stim) - prediction;