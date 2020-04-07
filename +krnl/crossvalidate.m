function residuals = crossvalidate(signal, t, partition, eventTimes, ...
    eventWindows, vectors, vectorWindows, lambda)

A = krnl.getToeplitz(t, eventTimes, ...
    eventWindows, vectors, vectorWindows);

parts = ones(1, partition) .* floor(length(t) / partition);
exc = mod(length(t), partition);
parts(1:exc) = parts(1:exc) + 1;

% residuals = NaN(size(signal));
residuals = cell(1, partition);
parfor p = 1:partition
    testInd = sum(parts(1:p-1)) + (1:parts(p));
    trainInd = setdiff(1:length(t), testInd);
    if isempty(A)
        prediction = ones(length(testInd),1) .* nanmean(signal(trainInd));
    else
        [B, fitInfo] = lasso(A(trainInd,:), signal(trainInd), 'Lambda', lambda);
        prediction = [ones(length(testInd),1), A(testInd,:)] * ...
            [fitInfo.Intercept; B];
    end
%     residuals(testInd) = prediction - signal(testInd);
    residuals{p} = prediction - signal(testInd);
end
residuals = cat(1, residuals{:});