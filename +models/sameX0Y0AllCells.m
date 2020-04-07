function error = sameX0Y0AllCells(respTraining, varsTraining, ...
    respTest, varsTest)
% respTraining  {neuron x 1}; each entry: [trial x stimulus]; contains 1
%               NaN for 1 neuron
% varsTraining  {neuron x 1}; each entry: [trial x stimulus]; contains 1
%               NaN for 1 neuron
% respTest      {neuron x 1}; 1 non-empty entry: [trial x stimulus]; all
%               NaNs except 1 entry
% varsTest      {neuron x 1}; 1 non-empty entry: [trial x stimulus]; all
%               NaNs except 1 entry

neuron = find(~cellfun(@isempty, respTest));
[trial,stim] = find(~isnan(respTest{neuron}));
rTrain = cellfun(@num2cell, respTraining, ...
    num2cell(ones(length(respTraining),1)), 'UniformOutput', false);
rTrain{neuron}{stim}(trial) = [];
rTrain = cat(2, rTrain{:});
vTrain = cellfun(@num2cell, varsTraining, ...
    num2cell(ones(length(varsTraining),1)), 'UniformOutput', false);
vTrain{neuron}{stim}(trial) = [];
vTrain = cat(2, vTrain{:});

coeffs = general.bilinfit(vTrain, rTrain);

ind = (neuron-1)*size(respTraining{1},2) + stim;
prediction = coeffs(ind,1) * varsTest{neuron}(trial, stim) + coeffs(ind,2);
error = respTest{neuron}(trial,stim) - prediction;