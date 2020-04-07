function error = nonvisOnlySingleCell(respTraining, varsTraining, ...
    respTest, varsTest)
% respTraining  [trial x stimulus]; contains 1 NaN
% varsTraining  [trial x stimulus]; contains 1 NaN
% respTest      [trial x stimulus]; all NaN, except one entry
% varsTest      [trial x stimulus]; all NaN, except one entry

rTrain = respTraining(:);
rTrain(isnan(rTrain)) = [];
vTrain = varsTraining(:);
vTrain(isnan(vTrain)) = [];

x = [vTrain, ones(size(vTrain))] \ rTrain;

[trial,stim] = find(~isnan(respTest));
prediction = x(1) * varsTest(trial,stim) + x(2);
error = respTest(trial,stim) - prediction;