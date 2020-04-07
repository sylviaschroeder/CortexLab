function [parameters, blankResponses, predictions, R2] = ...
    fitTuningCurve(responses, stimDirections, blanks, conditions, fixedPars)

% responses         [stimuli x repetitions]; response strength for each
%                   trial
% stimDirections    [stimuli x 2]; 1st col: direction, 2nd col: stim. ID
% blanks            [n x 1]; stim. IDs of blank stimuli
% conditions        [stimuli x repetitions]; integer for each trial
%                   differentiation between different conditions, such as 
%                   running and not running
% fixedPars         [1 x p]; parameters that will be fixed across
%                   conditions, default: [1 5] -> pref. dir, width

% parameters        [parameters x conditions]; different orientation tuning
%                   parameters for each condition; parameters: (1)
%                   preferred direction, (2) ampl. at pref. direction, (3)
%                   ampl. at opp. direction, (4) offset of tuning curve
%                   from zero, (5) tuning width
% blankResponses    [1 x conditions]; 

numFitIterations = 10;

if nargin < 5
    fixedPars = [1 5];
end

cUnique = unique(conditions(:));
parameters = NaN(5, length(cUnique));
blankResponses = NaN(1, length(cUnique));

dirs = repmat(stimDirections(:,1),1,size(responses,2));
resp = responses(stimDirections(:,2),:);
conds = conditions(stimDirections(:,2),:);

parInput = NaN(1,5);
if ~isempty(fixedPars)
    ind = ~isnan(resp);
    pars = fitoriWrapped(dirs(ind), resp(ind), [], [], '', numFitIterations);
    parInput(fixedPars) = pars(fixedPars);
end

pred = NaN(size(resp));

for c = 1:length(cUnique)
    ind = ~isnan(resp) & conds == cUnique(c);
    pars = fitoriWrapped(dirs(ind), resp(ind), [], parInput, '', ...
        numFitIterations);
    parameters(:,c) = pars;
    
    pred(ind) = orituneWrapped(pars, dirs(ind));
    
    r = responses(blanks,:);
    r(conditions(blanks,:)~=cUnique(c)) = NaN;
    blankResponses(c) = nanmean(r(:),1);
end

predictions = NaN(size(responses));
predictions(stimDirections(:,2),:) = pred;
R2 = 1 - nansum((resp(:)-pred(:)).^2) / nansum((resp(:)-nanmean(resp(:))).^2);