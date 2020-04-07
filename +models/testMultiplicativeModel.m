function error = testMultiplicativeModel(respTraining, varsTraining, ...
    respTest, varsTest)
% respTraining  [trial x stimulus]; contains 1 NaN
% varsTraining  [trial x stimulus]; contains 1 NaN
% respTest      [trial x stimulus]; all NaN, except one entry
% varsTest      [trial x stimulus]; all NaN, except one entry

xx = mat2cell(varsTraining,size(respTraining,1),ones(1,size(respTraining,2)));
yy = mat2cell(respTraining,size(respTraining,1),ones(1,size(respTraining,2)));

[trial, stim] = find(~isnan(respTest));
xx{stim}(trial) = [];
yy{stim}(trial) = [];

coeffs = models.bilinfit(xx,yy);

prediction = coeffs(stim,1) * varsTest(trial,stim) + coeffs(stim,2);
error = prediction - respTest(trial,stim);