function error = fitTuningCurveConditionsAndLaser_forCrossVal(respTraining, ...
    respTest, variables, testIndices)
% respTraining  [trial x stimulus]; contains 1 NaN
% respTest      [trial x stimulus]; all NaN, except one entry
% varsTest      {1 x 4}, {1}: [stimulus x 2], 1st col: direction, 2nd col:
%                             stimID
%                        {2}: [trial x stimulus], condition of each trial
%                        {3}: [1 x stimulus], laser on or not
%                        {4}: [1 x p], parameters that are fixed across
%                        conditions

% ind = find(isnan(respTraining));
% [~,stim] = find(isnan(respTraining));
directions = variables{1};
conditions = variables{2};
laserOn = variables{3};
fixedPars = variables{4};

[parsLaserOff, parsLaserOn] = gratings.fitTuningCurveConditionsAndLaser( ...
    respTraining', directions, conditions', laserOn', fixedPars, 1);

if nargin < 4
    testIndices = find(isnan(respTraining));
end
[~,stim] = ind2sub(size(respTraining), testIndices);

conds = conditions(testIndices);
condsUni = unique(conds(~isnan(conds)));
prediction = NaN(length(testIndices),1);

j = find(~laserOn(stim));
c = conds(j);
cUni = unique(c(~isnan(c)));
if length(cUni) < length(condsUni)
    if cUni(1) == condsUni(1)
        parsLaserOff(6:10) = [];
    else
        parsLaserOff = parsLaserOff(1:5) + parsLaserOff(6:10);
    end
end
k = isnan(c);
c(k) = [];
j(k) = [];
prediction(j) = gratings.orituneWrappedConditions(parsLaserOff, ...
    directions(stim(j)), c);

j = find(laserOn(stim));
c = conds(j);
cUni = unique(c(~isnan(c)));
if length(cUni) < length(condsUni)
    if cUni(1) == condsUni(1)
        parsLaserOn(6:10) = [];
    else
        parsLaserOn = parsLaserOn(1:5) + parsLaserOn(6:10);
    end
end
k = isnan(c);
c(k) = [];
j(k) = [];
prediction(j) = gratings.orituneWrappedConditions(parsLaserOn, ...
    directions(stim(j)), c);
test = respTest(testIndices);
error = prediction - test;




% resp = respTraining(:,directions(:,2));
% conds = conditions(:,directions(:,2));
% ind = find(isnan(resp));
% [~,stim] = find(isnan(resp));
% condsUni = unique(conds(~isnan(conds)));
% prediction = NaN(length(ind),1);
% j = find(~laserOn(stim));
% c = conds(ind(j));
% cUni = unique(c(~isnan(c)));
% if length(cUni) < length(condsUni)
%     if cUni(1) == condsUni(1)
%         parsLaserOff(6:10) = [];
%     else
%         parsLaserOff = parsLaserOff(1:5) + parsLaserOff(6:10);
%     end
% end
% k = isnan(c);
% c(k) = [];
% j(k) = [];
% prediction(j) = orituneWrappedConditions(parsLaserOff, ...
%     directions(stim(j),1), c);
% j = find(laserOn(stim));
% c = conds(ind(j));
% cUni = unique(c(~isnan(c)));
% if length(cUni) < length(condsUni)
%     if cUni(1) == condsUni(1)
%         parsLaserOn(6:10) = [];
%     else
%         parsLaserOn = parsLaserOn(1:5) + parsLaserOn(6:10);
%     end
% end
% k = isnan(c);
% c(k) = [];
% j(k) = [];
% prediction(j) = orituneWrappedConditions(parsLaserOn, ...
%     directions(stim(j),1), c);
% % for c = 1:length(cs)
% %     j = conds(ind)==cs(c) & ~laserOn(stim);
% %     prediction(j) = orituneWrappedConditions(parsLaserOff(:,c), directions(stim(j),1), ...
% %         conds(ind(j)));
% % end
% test = respTest(:,directions(:,2));
% err = prediction - test(ind);
% error = NaN(size(respTraining,2),1);
% error(directions(:,2),:) = err;