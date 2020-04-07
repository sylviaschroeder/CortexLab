function lambda = searchOptimalLambda(X, Y, doPlot)
% Determines optimal lambda in ridge regression: problem Y = X * A; or
% rather Y = a0 + X * A
% Solution: A = (X'*X + lambda*I)^(-1) * X'*Y

max_lambda = 5;
stepSize = 0.2;
crossValSets = 10; % data will be divided in that many parts, each of 
                    % which will be used for validation after fitting the rest

% don't make plots if not stated
if nargin < 3
    doPlot = 0;
end
                    
% center inputs
% X_ = bsxfun(@minus, X, mean(X, 1));

% estimate intercept a0 and subtract
% a0 = mean(Y);
% Y_ = bsxfun(@minus, Y, a0);

X_ = X;
Y_ = Y;

% perform cross-validation for each possible lambda value
lambdaList = stepSize : stepSize : max_lambda;
residuals = zeros(size(lambdaList));
setSize = round(size(X,1) / crossValSets);
for lam = 1:length(lambdaList)
    fprintf('Lambda: %.1f, set (of %d)', lambdaList(lam), crossValSets)
    res = zeros(1, crossValSets);
    for set = 1:crossValSets
        fprintf(' %d', set)
        % consider only data outside the validation set
        ind = ones(1,size(X,1));
        ind((set-1)*setSize+1 : min(size(X,1), set*setSize)) = 0;
        ind = logical(ind);
        X_set = X_(ind, :);
        Y_set = Y_(ind, :);
        % perform ridge regression on subset
        X_set = vertcat(X_set, lambdaList(lam) * eye(size(X_set,2)));
        Y_set = vertcat(Y_set, zeros(size(X_set,2), size(Y_set,2)));
        A = X_set \ Y_set;
%         A = pinv(X_set' * X_set + lambdaList(lam) * eye(size(X_set,2))) * ...
%             (X_set' * Y_set);
        % predict Y on validation set
        Y_pred = X_(~ind,:) * A;
%         Y_pred = X_ * A;
        % average squared residuals across rows (e.g. time) for each column
        % (e.g. neuron), then average across columns
        resAv = sum((Y_pred - Y_(~ind,:)) .^ 2, 1) / sum(ind);
%         resAv = sum((Y_pred - Y_) .^ 2, 1) / sum(ind);
        res(set) = mean(resAv);
    end
    fprintf('\n')
    % average across all validation sets
    residuals(lam) = mean(res);
end

% determine minimal residuals and thus optimal lambda
[~, ind] = min(residuals);
lambda = lambdaList(ind);

if doPlot == 1
    figure
    plot(lambdaList, residuals, 'k')
    xlabel('lambda')
    ylabel('residuals (average across validation sets, columns and rows)')
end