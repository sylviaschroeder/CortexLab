function [parsLaserOff, parsLaserOn, predLaserOff, predLaserOn, ...
    blanksLaserOff, blanksLaserOn] = ...
    fitTuningCurveConditionsAndLaser(responses, stimDirections, ...
    conditions, laserOn, fixedPars, runFitTwice)

if nargin < 6
    runFitTwice = 0;
end

% resp = responses(stimDirections(:,2),:);
% notBlanks = ~ismember(1:size(responses,1), blanks)';
notBlanks = ~isnan(stimDirections);
resp = responses(notBlanks, :);

% 1. find common preferred orientation
if ismember(1, fixedPars)
    ind = ~isnan(resp);
    dirs = repmat(stimDirections(notBlanks), 1, size(resp,2));
    parsInit = fitoriWrapped(dirs(ind), resp(ind));
    prefDir = mod(parsInit(1),360);
else
    prefDir = NaN;
end

% 2. fit (a) laser-off data and (b) laser-on data for two behavioural
% conditions
condInds = {find(~laserOn & notBlanks), find(laserOn & notBlanks)};
predictions = cell(1,2);
parameters = cell(1,2);
for k = 1:2
    ind = condInds{k};
    dirs = repmat(stimDirections, 1, size(responses,2));
    dirs = dirs(ind);
    r = responses(ind);
    c = conditions(ind);
    ind = ~isnan(r) & ~isnan(c);
    % find starting values for parameters fixed across conditions
    pars = [prefDir, NaN, NaN, NaN, NaN];
    parsInit = fitoriWrapped(dirs(ind), r(ind), [], pars);
    % find starting values for free parameters in condition 1
    conds = unique(c);
    pars(fixedPars) = parsInit(fixedPars);
    ind1 = ind & c == conds(1);
    parsInit1 = fitoriWrapped(dirs(ind1), r(ind1), [], pars);
    parsInit1(3) = (parsInit1(2)-parsInit1(3))/(parsInit1(2)+parsInit1(3)); % DI
    % find starting values for free parameters in condition 2
    ind2 = ind & c == conds(2);
    parsInit2 = fitoriWrapped(dirs(ind2), r(ind2), [], pars);
    parsInit2(3) = (parsInit2(2)-parsInit2(3))/(parsInit2(2)+parsInit2(3)); % DI
    if ismember(3, fixedPars)
        DI = mean([parsInit1(3) parsInit2(3)]);
        parsInit1(3) = DI;
        parsInit2(3) = DI;
    end
    parsInit = [parsInit1, parsInit2 - parsInit1];
    paramDeltas = NaN(1,5);
    paramDeltas(fixedPars) = 0;
    paramLimits = repmat([-Inf;Inf],1,5);
    paramLimits(:,1) = prefDir; % set pref. dir.
    paramLimits(1,3) = 0; % limit min. DI to 0 (not -1) so that pref. dir. does not change
    paramLimits(1,5) = median(diff(unique(dirs))); % limit min. sigma to the sample interval between tested orientations/directions
    [parameters{k}, err] = gratings.fitoriConditions(dirs(ind), r(ind), c(ind), ...
        paramDeltas, paramLimits, 10, parsInit);
    % run fit 2nd time
    if runFitTwice == 1
        parsInit(fixedPars) = parameters{k}(fixedPars);
        paramLimits(:, fixedPars) = repmat(parameters{k}(fixedPars), 2, 1);
        [p2, err2] = gratings.fitoriConditions(dirs(ind), r(ind), c(ind), ...
            paramDeltas, paramLimits, 10, parsInit);
        if err > err2
            parameters{k} = p2;
        end
    end
    predictions{k} = NaN(size(r));
    predictions{k}(ind) = gratings.orituneWrappedConditions(parameters{k}, dirs(ind), c(ind));
end
parsLaserOff = parameters{1};
parsLaserOn = parameters{2};
predLaserOff = predictions{1};
predLaserOn = predictions{2};

blanksLaserOff = cell(1,2);
blanksLaserOff{1,1} = responses(~laserOn & repmat(~notBlanks,1, ...
    size(responses,2)) & conditions==1);
blanksLaserOff{1,2} = responses(~laserOn & repmat(~notBlanks,1, ...
    size(responses,2)) & conditions==2);
blanksLaserOn = cell(1,2);
blanksLaserOn{1,1} = responses(laserOn & repmat(~notBlanks,1, ...
    size(responses,2)) & conditions==1);
blanksLaserOn{1,2} = responses(laserOn & repmat(~notBlanks,1, ...
    size(responses,2)) & conditions==2);



% % 2. fit laser-off data for two behavioural conditions
% ind = find(~laserOn & notBlanks);
% dirs = repmat(stimDirections(ind), 1, size(responses,2));
% r = responses(ind,:);
% c1 = conditions(ind,:);
% ind = ~isnan(r) & ~isnan(c1);
% % find starting values for parameters fixed across conditions
% pars = [prefDir, NaN, NaN, NaN, NaN];
% parsInit = fitoriWrapped(dirs(ind), r(ind), [], pars);
% % find starting values for free parameters in condition 1
% conds = unique(c1);
% pars(fixedPars) = parsInit(fixedPars);
% ind1 = ind & c1 == conds(1);
% parsInit1 = fitoriWrapped(dirs(ind1), r(ind1), [], pars);
% parsInit1(3) = (parsInit1(2)-parsInit1(3))/(parsInit1(2)+parsInit1(3)); % DI
% % find starting values for free parameters in condition 2
% ind2 = ind & c1 == conds(2);
% parsInit2 = fitoriWrapped(dirs(ind2), r(ind2), [], pars);
% parsInit2(3) = (parsInit2(2)-parsInit2(3))/(parsInit2(2)+parsInit2(3)); % DI
% if ismember(3, fixedPars)
%     DI = mean([parsInit1(3) parsInit2(3)]);
%     parsInit1(3) = DI;
%     parsInit2(3) = DI;
% end
% parsInit = [parsInit1, parsInit2 - parsInit1];
% paramDeltas = NaN(1,5);
% paramDeltas(fixedPars) = 0;
% paramLimits = repmat([-Inf;Inf],1,5);
% paramLimits(:,1) = prefDir; % set pref. dir.
% paramLimits(1,3) = 0; % limit min. DI to 0 (not -1) so that pref. dir. does not change
% paramLimits(1,5) = median(diff(unique(dirs))); % limit min. sigma to the sample interval between tested orientations/directions
% [parsLaserOff, err] = fitoriConditions(dirs(ind), r(ind), c1(ind), ...
%     paramDeltas, paramLimits, 10, parsInit);
% % run fit 2nd time
% if runFitTwice == 1
%     parsInit(fixedPars) = parsLaserOff(fixedPars);
%     paramLimits(:, fixedPars) = repmat(parsLaserOff(fixedPars), 2, 1);
%     [parsLaserOff2, err2] = fitoriConditions(dirs(ind), r(ind), c1(ind), ...
%         paramDeltas, paramLimits, 10, parsInit);
%     if err > err2
%         parsLaserOff = parsLaserOff2;
%     end
% end
% predLaserOff = NaN(size(r));
% predLaserOff(ind) = orituneWrappedConditions(parsLaserOff, dirs(ind), c1(ind));
% 
% % 3. fit laser-on data for two behavioural conditions
% ind = find(laserOn & notBlanks);
% dirs = repmat(stimDirections(ind), 1, size(responses,2));
% r = responses(ind,:);
% c2 = conditions(ind,:);
% ind = ~isnan(r) & ~isnan(c2);
% % find starting values for parameters fixed across conditions
% pars = [prefDir, NaN, NaN, NaN, NaN];
% parsInit = fitoriWrapped(dirs(ind), r(ind), [], pars);
% % find starting values for free parameters in condition 1
% conds = unique(c2);
% pars(fixedPars) = parsInit(fixedPars);
% ind1 = ind & c2 == conds(1);
% parsInit1 = fitoriWrapped(dirs(ind1), r(ind1), [], pars);
% parsInit1(3) = (parsInit1(2)-parsInit1(3))/(parsInit1(2)+parsInit1(3)); % DI
% % find starting values for free parameters in condition 2
% ind2 = ind & c2 == conds(2);
% parsInit2 = fitoriWrapped(dirs(ind2), r(ind2), [], pars);
% parsInit2(3) = (parsInit2(2)-parsInit2(3))/(parsInit2(2)+parsInit2(3)); % DI
% if ismember(3, fixedPars)
%     DI = mean([parsInit1(3) parsInit2(3)]);
%     parsInit1(3) = DI;
%     parsInit2(3) = DI;
% end
% parsInit = [parsInit1, parsInit2 - parsInit1];
% paramDeltas = NaN(1,5);
% paramDeltas(fixedPars) = 0;
% paramLimits = repmat([-Inf;Inf],1,5);
% paramLimits(:,1) = prefDir; % set pref. dir.
% paramLimits(1,3) = 0; % limit min. DI to 0 (not -1) so that pref. dir. does not change
% paramLimits(1,5) = median(diff(unique(dirs))); % limit min. sigma to the sample interval between tested orientations/directions
% [parsLaserOn, err] = fitoriConditions(dirs(ind), r(ind), c2(ind), ...
%     paramDeltas, paramLimits, 10, parsInit);
% % run fit 2nd time
% if runFitTwice == 1
%     parsInit(fixedPars) = parsLaserOn(fixedPars);
%     paramLimits(:, fixedPars) = repmat(parsLaserOn(fixedPars), 2, 1);
%     [parsLaserOn2, err2] = fitoriConditions(dirs(ind), r(ind), c2(ind), ...
%         paramDeltas, paramLimits, 10, parsInit);
%     if err > err2
%         parsLaserOn = parsLaserOn2;
%     end
% end
% predLaserOn = NaN(size(r));
% predLaserOn(ind) = orituneWrappedConditions(parsLaserOn, dirs(ind), c2(ind));