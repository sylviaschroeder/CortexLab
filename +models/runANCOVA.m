function [ancovas, ancova3] = runANCOVA(nonVisual, responses, baselines, ...
    blankStim)

nonvis = nonVisual(:);
resp = responses(:);
base = baselines(:);
stim = reshape(repmat(1:size(responses,2), size(responses, 1), 1), [], 1);

ancovas = [];
% (1) Run ANCOVA just with variables stimuli and nonvisual signal
% % all stimuli, without baseline
% [~,atab] = aoctool(nonvis, resp, stim, 0.05, '', '', '', 'off', ...
%     'separate lines');
% ancovas.allStimNoBase.stim = atab{2,6};
% ancovas.allStimNoBase.nonvis = atab{3,6};
% ancovas.allStimNoBase.interact = atab{4,6};
% % all stimuli, with baseline
% [~,atab] = aoctool(nonvis, resp + base, stim, 0.05, '', '', '', 'off', ...
%     'separate lines');
% ancovas.allStimWithBase.stim = atab{2,6};
% ancovas.allStimWithBase.nonvis = atab{3,6};
% ancovas.allStimWithBase.interact = atab{4,6};
% 
% % no blank, without baseline
% noBlanks = stim ~= blankStim;
% [~,atab] = aoctool(nonvis(noBlanks), resp(noBlanks), stim(noBlanks), 0.05, ...
%     '', '', '', 'off', 'separate lines');
% ancovas.noBlankNoBase.stim = atab{2,6};
% ancovas.noBlankNoBase.nonvis = atab{3,6};
% ancovas.noBlankNoBase.interact = atab{4,6};
% % no blank, with baseline
% [~,atab] = aoctool(nonvis(noBlanks), resp(noBlanks) + base(noBlanks), ...
%     stim(noBlanks), 0.05, '', '', '', 'off', 'separate lines');
% ancovas.noBlankWithBase.stim = atab{2,6};
% ancovas.noBlankWithBase.nonvis = atab{3,6};
% ancovas.noBlankWithBase.interact = atab{4,6};

% (2) Run ANCOVA with 3 variables: nonvisual signal, baseline, stimulus
tbl = table(nonvis, base, stim, resp);
mdl = cell(1,7);
mdl{1} = fitlm(tbl, 'interactions', 'CategoricalVars', 3);
mdl{2} = fitlm(tbl,'resp ~ 1 + base * stim','CategoricalVars',3);
mdl{3} = fitlm(tbl,'resp ~ 1 + nonvis * stim','CategoricalVars',3);
mdl{4} = fitlm(tbl,'resp ~ 1 + nonvis * base','CategoricalVars',3);
mdl{5} = fitlm(tbl,'resp ~ 1 + nonvis','CategoricalVars',3);
mdl{6} = fitlm(tbl,'resp ~ 1 + base','CategoricalVars',3);
mdl{7} = fitlm(tbl,'resp ~ 1 + stim','CategoricalVars',3);
mdl{8} = fitlm(tbl,'resp ~ 1 + stim*(nonvis + base)','CategoricalVars',3);

Rsquares = [mdl{1}.Rsquared.Adjusted, mdl{2}.Rsquared.Adjusted, ...
    mdl{3}.Rsquared.Adjusted, mdl{4}.Rsquared.Adjusted, ...
    mdl{5}.Rsquared.Adjusted, mdl{6}.Rsquared.Adjusted, ...
    mdl{7}.Rsquared.Adjusted, mdl{8}.Rsquared.Adjusted];
ancova3.adjRsquares = Rsquares;
ancova3.models = {'full', 'no nonvis', 'no basel', 'no stim', ...
    'only nonvis', 'only basel', 'only stim', 'full w/o interact'};

atab = mdl{1}.anova;
if all(atab.pValue > 0.05)
    ancova3.best =[];
else
    [~,maxInd] = max(Rsquares);
    ancova3.best = mdl{maxInd};
end