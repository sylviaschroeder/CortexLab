function [thresholds, r_squared] = ...
    findBestThreshold(neuralData, nonVisData, nThresholds)
% Finds the best threshold of non-visual data (e.g. running speed) to fit
% population responses to their means (one set above threshold, the other
% below threshold)

% nonVisData    [1 x time]
% nThresholds   number of thresholds to partition data optimally

% thresholds    [nThresholds x 1]; optimal thresholds for non-visual data
% r_squared     double; sum of squared residuals

smoothHalf = 3;

if nargin < 3 || isempty(nThresholds)
    nThresholds = 1;
end

% calculate residuals for all possible combinations of thresholds
values = unique(nonVisData);
values = prctile(values, 2:2:98);
residuals = NaN(length(values)^nThresholds, 1);
threshs = [ones(1, nThresholds-1) 0]; % first set of thresholds
k = 0;
lastValid = getResiduals(neuralData, nonVisData, values(1:nThresholds));
while threshs(1) <= length(values)-nThresholds+1
    k = k+1;
    if k > numel(residuals)
        break
    end
    next = true;
    j = 0;
    while next
        if j>0
            threshs(end-j+1) = 1;
        end
        threshs(end-j) = threshs(end-j)+1;
        if threshs(end-j) <= length(values)
            next = false;
        else
            j = j+1;
        end
    end
    if length(unique(threshs)) < nThresholds || any((sort(threshs)-threshs) ~= 0)
        residuals(k) = lastValid;
        continue
    end
    residuals(k) = getResiduals(neuralData, nonVisData, values(threshs));
    lastValid = residuals(k);
end
residuals(k+1:end) = lastValid;
if nThresholds > 1
    residuals = reshape(residuals, ones(1,nThresholds).*length(values));
end
% smooth lanscape of residuals
gauss = normpdf(-smoothHalf:smoothHalf,0,1);
wind = gauss';
cutInd0 = [false(smoothHalf,1); true(length(values),1); false(smoothHalf,1)];
cutInd = cutInd0;
if nThresholds >= 2
    wind = wind * wind';
    cutInd = bsxfun(@times, cutInd, cutInd');
end
for k = 3:nThresholds
    wind = bsxfun(@times, wind, reshape(gauss,ones(1, ndims(wind)),[]));
    cutInd = bsxfun(@times, cutInd, reshape(cutInd0,ones(1,ndims(cutInd)),[]));
end
res = residuals;
if nThresholds == 1
    res = [ones(smoothHalf,1)*residuals(1); residuals; ones(smoothHalf,1)*residuals(end)];
else
    siz = size(residuals);
    dims = 1:ndims(residuals);
    for d = 1:ndims(residuals)
        swap = [d dims(2:end)];
        swap(d) = 1;
        swappedRes = permute(res, swap);
        front = reshape(swappedRes(1,:), [1 siz(swap(2:end))]);
        back = reshape(swappedRes(end,:), [1 siz(swap(2:end))]);
        swappedRes = [repmat(front, [smoothHalf ones(1,length(dims)-1)]); ...
            swappedRes; repmat(back, [smoothHalf ones(1,length(dims)-1)])];
        res = permute(swappedRes, swap);
        siz(d) = size(swappedRes,1);
    end
end
residualsSm = convn(res, wind, 'same');
residualsSm = residualsSm(cutInd==1);
% if nThresholds == 1
%     residualsSm = smooth(residuals, 0.05*length(values), 'lowess');
% else
%     residualsSm = residuals;
% end
[~, indMin] = min(residualsSm(:));
indMin2 = indMin;
thresholds = NaN(1, nThresholds);
for k = nThresholds:-1:1
    ind = ceil(indMin2 / length(values)^(k-1));
    thresholds(k) = values(ind);
    indMin2 = indMin2 - (ind-1)*length(values)^(k-1);
end
thresholds = flip(thresholds);

r_squared = 1 - residuals(indMin) / sum((neuralData-mean(neuralData)).^2);

function residual = getResiduals(popResp, nonVisData, threshs)
y = zeros(size(popResp));
threshs = [min(nonVisData)-1, threshs, max(nonVisData)+1];
for k = 1:length(threshs)-1
    ind = nonVisData > threshs(k) & nonVisData <= threshs(k+1);
    y(ind) = mean(popResp(ind));
end
residual = sum((popResp - y) .^ 2);