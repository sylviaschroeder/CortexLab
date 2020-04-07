function [eta, triggeredTraces] = getETA(signal, toepl, numSamples)

% check necessary memory
if numel(signal)*size(toepl,2)*8/1e9 > 16 % sig will be larger than 16GB (toepl will also grow to 16 GB)
    % calculate averages cell by cell
    sig = zeros(size(signal,2), size(toepl,2));
    for iCell = 1:size(signal,2)
        s = repmat(signal(:,iCell), 1, size(toepl, 2));
        s(toepl==0) = NaN;
        sig(iCell,:) = nanmean(s,1);
    end
else
    sig = repmat(signal, 1, 1, size(toepl,2));
    sig(repmat(permute(toepl, [1 3 2]), 1, size(signal,2), 1)==0) = NaN;
    sig = permute(nanmean(sig,1), [2 3 1]); % [neurons x eventSamples]
end
eta = cell(1, length(numSamples));
for ev = 1:length(numSamples)
    eta{ev} = sig(:,sum(numSamples(1:ev-1))+(1:numSamples(ev)))';
end

if nargout > 1
    tt = cell(1, size(toepl,2));
    triggeredTraces = cell(1, length(numSamples));
    for k = 1:size(toepl,2)
        tt{k} = signal(toepl(:,k)==1,:);
    end
    [siz,~] = cellfun(@size, tt);
    for ev = 1:length(numSamples)
        ind = sum(numSamples(1:ev-1)) + (1:numSamples(ev));
        m = min(siz(ind));
        for k = find(siz(ind) > m)
            tt{ind(k)} = tt{ind(k)}(1:m,:);
        end
        triggeredTraces{ev} = permute(cat(3, tt{ind}), [3 2 1]); % [time x neurons x event occurrences]
    end
end