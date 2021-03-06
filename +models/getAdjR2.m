function adjR2 = getAdjR2(cellResponse, modelPars)
% getAdjR2(cellResponse, nonVisual, modelPars)
% modelPars     struct; output from nonVis.modelResponsePerTrial

% check how well alpha*kernel+baseline fits the data
predKernel = bsxfun(@plus, bsxfun(@times, modelPars.kernel, ...
    permute(modelPars.alphaEachTrial, [3 1 2])), ...
    permute(modelPars.baselineEachTrial, [3 1 2]));
adjR2.kernel = 1 - (sum((cellResponse(:) - predKernel(:)).^2) / ...
    (numel(cellResponse)-length(modelPars.kernel)-...
    numel(modelPars.alphaEachTrial)-...
    numel(modelPars.baselineEachTrial)-1)) / ...
    (sum((cellResponse(:) - mean(cellResponse(:))).^2) / ...
    (numel(cellResponse)-1));
% check how well alpha+(alpha*beta+gamma)*nonVisual fits the
% stimulus responses in each trial (alphaEachTrial)
% predStimResp = bsxfun(@plus, modelPars.alphas, ...
%     bsxfun(@times, modelPars.alphas.*modelPars.beta + ...
%     modelPars.gamma, nonVisual));
% adjR2.stimResp = 1 - (sum((modelPars.alphaEachTrial(:) - ...
%     predStimResp(:)).^2) / ...
%     (numel(modelPars.alphaEachTrial)-length(modelPars.alphas)-3)) / ...
%     (sum((modelPars.alphaEachTrial(:) - ...
%     mean(modelPars.alphaEachTrial(:))).^2) / ...
%     (numel(modelPars.alphaEachTrial)-1));
% % check how well alphaUncond+gammaUncond*nonVisual fits the
% % stimulus responses in each trial (alphaEachTrial)
% predStimRespUncon = bsxfun(@plus, modelPars.alphasUncond, ...
%     bsxfun(@times, modelPars.gammasUncond, nonVisual));
% adjR2.stimRespUncon = 1 - (sum((modelPars.alphaEachTrial(:) - ...
%     predStimRespUncon(:)).^2) / ...
%     (numel(modelPars.alphaEachTrial)- ...
%     2*length(modelPars.alphasUncond)-1)) / ...
%     (sum((modelPars.alphaEachTrial(:) - ...
%     mean(modelPars.alphaEachTrial(:))).^2) / ...
%     (numel(modelPars.alphaEachTrial)-1));
% % check how well delta*nonVisual+const fits baselines
% predBase = modelPars.delta.*nonVisual + modelPars.const;
% adjR2.base = 1 - (sum((modelPars.baselineEachTrial(:) - ...
%     predBase(:)).^2) / ...
%     (numel(modelPars.baselineEachTrial)-3)) / ...
%     (sum((modelPars.baselineEachTrial(:) - ...
%     mean(modelPars.baselineEachTrial(:))).^2) / ...
%     (numel(modelPars.baselineEachTrial)-1));
% % check how well the model fits the data (without taking
% % into account the temporal dynamics):
% % alpha+(alpha*beta+gamma)*nonVisual + delta*nonVisual+const =
% % alphaEachTrial + baselineEachTrial ???
% predStimAndBase = predStimResp + predBase;
% adjR2.stimAndBase = 1 - (sum((modelPars.alphaEachTrial(:) + ...
%     modelPars.baselineEachTrial(:) - predStimAndBase(:)).^2) / ...
%     (numel(modelPars.alphaEachTrial)-length(modelPars.alphas)-5)) / ...
%     (sum((modelPars.alphaEachTrial(:)+modelPars.baselineEachTrial(:) - ...
%     mean(modelPars.alphaEachTrial(:)+modelPars.baselineEachTrial(:))).^2) / ...
%     (numel(modelPars.alphaEachTrial)-1));
% % check how well the complete model fits the data (with
% % temporal dynamics):
% % kernel*(alpha+(alpha*beta+gamma)*nonVisual) + delta*nonVisual+const =
% % responses ???
% predTotal = bsxfun(@plus, modelPars.alphas, bsxfun(@times, modelPars.alphas, ...
%     modelPars.beta .* nonVisual)) + modelPars.gamma .* nonVisual;
% % (2) get temporal dynamics by multiplying with kernel
% predTotal = bsxfun(@times, modelPars.kernel, permute(predTotal, [3 1 2]));
% % (3) add baseline for each trial
% predTotal = bsxfun(@plus, predTotal, ...
%     permute(modelPars.delta .* nonVisual, [3 1 2])) + modelPars.const;
% adjR2.total = 1 - (sum((cellResponse(:)-predTotal(:)).^2) / ...
%     (numel(cellResponse)-length(modelPars.kernel)- ...
%     length(modelPars.alphas)-5)) / ...
%     (sum((cellResponse(:)-mean(cellResponse(:))).^2) / ...
%     (numel(cellResponse)-1));