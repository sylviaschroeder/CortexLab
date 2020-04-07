function error = fitTuningCurveConditions_forCrossVal(respTraining, ...
    respTest, variables, testIndices)
% respTraining  [trial x stimulus]; contains 1 NaN
% respTest      [trial x stimulus]; all NaN, except one entry
% varsTest      {1 x 4}, {1}: [stimulus x 2], 1st col: direction, 2nd col:
%                             stimID
%                        {2}: [n x 1], stim IDs of blank stimuli, n is
%                        number of conditions
%                        {3}: [trial x stimulus], condition of each trial
%                        {4}: [1 x p], parameters that are fixed across
%                        conditions

% ind = find(isnan(respTraining));
% [~,stim] = find(isnan(respTraining));
directions = variables{1};
blanks = variables{2};
conditions = variables{3};
fixedPars = variables{4};

parameters = gratings.fitTuningCurveConditions(respTraining', directions, blanks, ...
    conditions', fixedPars, 1);

if nargin < 4
    testIndices = find(isnan(respTraining));
end
[~,stim] = ind2sub(size(respTraining), testIndices);

conds = conditions(testIndices);
dirs = NaN(length(testIndices),1);
for k = 1:length(testIndices)
    j = find(stim(k) == directions(:,2));
    if ~isempty(j)
        dirs(k) = directions(j,1);
    end
end
invalid = isnan(dirs);
conds(invalid) = [];
dirs(invalid) = [];


% if nargin < 4
%     testIndices = isnan(respTraining);
% else
%     tmp = false(size(respTraining));
%     tmp(testIndices) = true;
%     testIndices = tmp;
% end
% 
% resp = respTraining(:,directions(:,2));
% conds = conditions(:,directions(:,2));
% ind = testIndices(:,directions(:,2));
% 
% [trials,stim] = find(ind);
% ind = find(ind);
cs = unique(conds);
prediction = NaN(length(dirs),1);
for c = 1:length(cs)
    j = conds==cs(c);
%     prediction(j) = orituneWrappedConditions(parameters(:,c), directions(stim(j),1), ...
%         conds(ind(j)));
    prediction(j) = gratings.orituneWrappedConditions(parameters(:,c), dirs(j), conds(j));
end
% test = respTest(:,directions(:,2));
% err = prediction - test(ind);
% error = NaN(size(respTraining,2),1);
% error(directions(:,2),:) = err;
test = respTest(testIndices);
test(invalid) = [];
error = NaN(length(testIndices), 1);
error(~invalid) = prediction - test;