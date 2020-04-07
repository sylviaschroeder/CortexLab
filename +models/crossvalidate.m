function errors = crossvalidate(model, neuralResponses, variables, mode)
% model             function that fits model and returns prediction error for
%                   test set (e.g. models.testFitlmModel)
% neuralResponses   {neurons x 1}; each entry: [trials x stimuli]
% variables         {neurons x numVars}; each entry: [trials x stimuli];
%                   if size(variables,1)==1, same data are used for each
%                   neurons
% (bestModel)       LinearModel from fitlm (in models.runANCOVA)

errors = cell(size(neuralResponses));
% n = nargin;
for iNeuron = 1:size(neuralResponses,1)
    resp = neuralResponses{iNeuron};
    [numTrials, numStim] = size(resp);
    switch mode
        case 'leaveOneOut'
            indsPerSet = num2cell((1:numel(resp))');
        case 'stimSets'
            trials = zeros(size(resp));
            for k = 1:numStim
                trials(:,k) = randperm(numTrials);
            end
            indsPerSet = cell(numTrials,1);
            for k = 1:numTrials
                indsPerSet{k} = find(trials == k);
            end
    end
    err = cell(size(indsPerSet));
    parfor k = 1:length(indsPerSet)
%     for k = 1:length(indsPerSet)
        respTraining = resp;
        respTraining(indsPerSet{k}) = NaN;
        respTest = NaN(size(resp));
        respTest(indsPerSet{k}) = resp(indsPerSet{k});
        err{k} = model(respTraining, respTest, variables, indsPerSet{k});
    end
    errs = NaN(size(resp));
    for k = 1:length(indsPerSet)
        errs(indsPerSet{k}) = err{k};
    end
    errors{iNeuron} = errs;
end