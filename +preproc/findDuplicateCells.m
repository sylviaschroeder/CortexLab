function [isDuplicate, duplicates] = findDuplicateCells(ROIpositions, F, ops)

% ROIpositions  {1 x planes}, each entry: [ROIs x 3], each row contains
%               positions of ROI (x,y,z) (in meta.ROI.CellXYZMicrons)
% F             {exp x planes}, each entry: [ROIs x time]; each column 
%               contains calcium trace of ROI, best to concatinate 
%               experiments (because some cells don't respond to some 
%               stimuli)
% ops           .minCorr: sets minimum correlation between ROI traces to
%                         classify ROIs as duplicates
%               .maxDistXY: sets maximum Euclidean distance between ROIs in 
%                         horizontal direction to classify as duplicates
%               .maxDistZ: sets maximum distance between ROIs in vertical
%                         direction to classify as duplicates
%
% isDuplicate   {1 x planes}, each entry: [ROIs x 1], each entry: is 1 if
%               cell appears again in other plane AND does NOT have
%               strongest signal, is 0 otherwise
% duplicates    struct (1 x planes)
%               .ROI: struct (ROI x 1)
%                     .plane: [duplicates x 1], plane of duplicate
%                     .ID: [duplicates x 1], ROI ID of duplicate
%                     .distance: [duplicates x 1], distance to duplicate
%                     .corrCoeff: [duplicates x 1], correlation coefficient
%                                 with duplicate
%               also cells that have duplicates AND have the strongest
%               signal have duplicates listed here

% Idea of algorithm: for each ROI find the closest ROI (in terms of
% Euclidean distance) in adjacent plane, check whether distance is smaller
% than ops.maxDist and whether correlation coefficient between calcium
% traces is at least ops.minCorr

if nargin < 3 || ~isfield(ops, 'minCorr')
    ops.minCorr = 0.4;
end
if nargin < 3 || ~isfield(ops, 'maxDistXY')
    ops.maxDistXY = 5;
end
if nargin < 3 || ~isfield(ops, 'maxDistZ')
    ops.maxDistZ = 30;
end
if nargin < 3 || ~isfield(ops, 'filtWindow')
    ops.filtWindow = 5;
end

isDuplicate = cell(1, length(ROIpositions));
duplicates = struct('ROI', cell(1, length(ROIpositions)));

dupl = cell(1, length(ROIpositions)); % each entry: [cells x 2], each row for one duplicate 
                   % cell: (plane of cell, ROI of cell)
for iPlane = 1:length(ROIpositions)
    isDuplicate{iPlane} = zeros(size(ROIpositions{iPlane},1),1);
    duplicates(iPlane).ROI = struct('plane',cell(size(ROIpositions{iPlane},1),1), ...
        'ID', cell(size(ROIpositions{iPlane},1),1), ...
        'distance', cell(size(ROIpositions{iPlane},1),1), ...
        'corrCoeff', cell(size(ROIpositions{iPlane},1),1));
    dupl{iPlane} = cell(size(ROIpositions{iPlane},1),1);
end

for iPlane = 1:length(ROIpositions)-1
    if isempty(ROIpositions{iPlane})
        continue
    end
    for jPlane = iPlane+1:length(ROIpositions)
        if isempty(ROIpositions{jPlane})
            continue
        end
        % calculate pairwise distances between all ROIs in plane i and
        % all ROIs in plane j
        distancesXY = sqrt(sum(bsxfun(@minus, ...
            permute(ROIpositions{iPlane}(:,1:2), [1 3 2]), ...
            permute(ROIpositions{jPlane}(:,1:2), [3 1 2])) .^ 2, 3));
        distancesZ = abs(bsxfun(@minus, ROIpositions{iPlane}(:,3), ...
            ROIpositions{jPlane}(:,3)'));
        % determine for each ROI in plane i, which cell in plane j is
        % the closest (in X-Y) and what the distance is (in X-Y and Z)
        [minDistsXY, inds] = min(distancesXY, [], 2);
        minDistsZ = distancesZ(sub2ind(size(distancesZ), (1:size(distancesZ,1))', inds));
        % check whether ROIs in plane i have several "neighbours" in 
        % plane j with the same minimal distance (in X-Y, Z should be very
        % similar anyway)
        minis = bsxfun(@eq, minDistsXY, distancesXY);
        % disregard any ROI in plane i that was identified as duplicate
        % already
        % disregard any ROI pairs with a distance larger than the set
        % thresholds (in X-Y and Z)
        minis(minDistsXY > ops.maxDistXY,:) = false;
        minis(distancesZ > ops.maxDistZ) = false;
        inds(minDistsXY > ops.maxDistXY) = NaN;
        inds(minDistsZ > ops.maxDistZ) = NaN;
        % generate a matrix in which row k contains all closest neighbours (in
        % plane j) of ROI k in plane i
        numMinis = sum(minis,2);
        if any(numMinis > 1)
            inds2 = NaN(size(inds,1), max(numMinis));
            inds2(:,1) = inds;
            for k = find(numMinis > 1)
                inds2(k,1:numMinis(k)) = find(minis(k,:));
            end
            inds = inds2;
        end
        % if ROI of plane j was detected as nearest neighbour several times,
        % associate it only with closest ROI in plane i
        multDetected = bsxfun(@eq, reshape(inds',[],1), reshape(inds',1,[]));
        multDetected(eye(size(multDetected,1))==1) = 0;
        [a, b] = find(multDetected);
        % all ROIs in plane i for which ROI in plane j is the closest
        neighb1 = ceil(a / size(inds,2));
        neighb2 = ceil(b / size(inds,2));
        % ROIs in plane j that have multiple neighbours in plane i
        ROIind = mod(a, size(inds,2));
        ROIind(ROIind==0) = size(inds,2);
        ROIid = inds(sub2ind(size(inds), neighb1, ROIind));
        % disregard duplicate detections
        d = neighb1 > neighb2;
        neighb1(d) = [];
        neighb2(d) = [];
        ROIid(d) = [];
        done2 = zeros(length(ROIid),1);
        % find all neighbors in current plane close to the same ROI in next
        % plane, only keep the pair that's closest to each other
        for k = 1:length(neighb1)
            if done2(k) == 1
                continue
            end
            instances = find(ROIid == ROIid(k));
            done2(instances) = 1;
            neighbs = unique([neighb1(instances); neighb2(instances)]);
            dists = minDistsXY(neighbs);
            [~,closest] = min(dists);
            neighbs(closest) = [];
            for j = 1:length(neighbs)
                inds(neighbs(j),inds(neighbs(j),:)==ROIid(k)) = NaN;
            end
        end
        % calculate pairwise correlations between traces of all pairs
        corrs = NaN(size(inds));
        for k = find(any(~isnan(inds),2))'
            for m = 1:numMinis(k)
                c = zeros(1, size(F,1));
                sgnlToNs = zeros(2, size(F,1));
                for exp = 1:size(F,1)
                    if isempty(F{exp,iPlane})
                        continue
                    end
                    c(exp) = corr(medfilt1(F{exp,iPlane}(:,k),ops.filtWindow), ...
                        medfilt1(F{exp,jPlane}(:,inds(k,m)), ops.filtWindow));
                    sgnlToNs(1,exp) = diff(prctile(medfilt1( ...
                        F{exp,iPlane}(:,k),ops.filtWindow), [50 98])) / ...
                        mad(F{exp,iPlane}(:,k),1);
                    sgnlToNs(2,exp) = diff(prctile(medfilt1( ...
                        F{exp,jPlane}(:,inds(k,m)),ops.filtWindow), [50 98])) / ...
                        mad(F{exp,jPlane}(:,inds(k,m)),1);
                end
                weights = min(sgnlToNs .^ 2, [], 1);
                weights = weights ./ sum(weights);
                corrs(k,m) = sum(c .* weights);
            end
        end
        % disregard pairs with correlations smaller than the set threshold
        corrs(corrs < ops.minCorr) = NaN;
        % collect all pairs
        for k = find(any(~isnan(corrs),2))'
            [~,indC] = max(corrs(k,:));
            dupl{iPlane}{k} = [dupl{iPlane}{k}; [jPlane inds(k,indC)]];
        end
    end
end

F_ = cell(1, size(F,2));
for iPlane = 1:length(dupl)
    F_{iPlane} = cat(1, F{:,iPlane});
end
for iPlane = 1:length(dupl)
    cells = find(~cellfun(@isempty, dupl{iPlane}));
    for iCell = 1:length(cells)
        sameCell = NaN(1, length(dupl));
        sameCell(iPlane) = cells(iCell);
        for p = iPlane:length(dupl)
            id = sameCell(p);
            if isnan(id)
                continue
            end
            d = dupl{p}{id};
            dupl{p}{id} = [];
            for k = 1:size(d,1)
                sameCell(d(k,1)) = d(k,2);
            end
        end
        ind = find(~isnan(sameCell));
        d = [ind', sameCell(ind)'];
        signalToNoise = NaN(size(d,1),1);
        for k = 1:size(d,1)
            trace = F_{d(k,1)}(:,d(k,2));
            signal = diff(prctile(medfilt1(trace, ops.filtWindow),[50,98]));
            noise = mad(trace, 1);
            signalToNoise(k) = signal / noise;
%             signalToNoise(iCell) = median(trace);
        end
        pl = zeros(1, size(d,1));
        ID = zeros(1, size(d,1));
        dist = NaN(size(d,1), size(d,1));
        correlations = NaN(size(d,1), size(d,1));
        for k = 1:size(d,1)
            pl(k) = d(k,1);
            ID(k) = d(k,2);
            for jCell = 1:size(d,1)
                if jCell == k
                    continue
                end
                dist(jCell,k) = sqrt(sum((ROIpositions{pl(k)}(ID(k),:) - ...
                    ROIpositions{d(jCell,1)}(d(jCell,2),:)) .^2));
                
                c = zeros(1, size(F,1));
                sgnlToNs = zeros(2, size(F,1));
                for exp = 1:size(F,1)
                    if isempty(F{exp,pl(k)})
                        continue
                    end
                    c(exp) = corr(medfilt1(F{exp,pl(k)}(:,ID(k)),ops.filtWindow), ...
                        medfilt1(F{exp,d(jCell,1)}(:,d(jCell,2)), ops.filtWindow));
                    sgnlToNs(1,exp) = diff(prctile(medfilt1( ...
                        F{exp,pl(k)}(:,ID(k)),ops.filtWindow), [50 98])) / ...
                        mad(F{exp,pl(k)}(:,ID(k)),1);
                    sgnlToNs(2,exp) = diff(prctile(medfilt1( ...
                        F{exp,d(jCell,1)}(:,d(jCell,2)), ...
                        ops.filtWindow), [50 98])) / ...
                        mad(F{exp,d(jCell,1)}(:,d(jCell,2)),1);
                end
                weights = min(sgnlToNs .^ 2, [], 1);
                weights = weights ./ sum(weights);
                correlations(jCell,k) = sum(c .* weights);
            end
        end
        dupls = find(nanmean(correlations,1) >= ops.minCorr);
        d = d(dupls,:);
        dist = dist(dupls, dupls);
        correlations = correlations(dupls, dupls);
        [~,bestCell] = max(signalToNoise(dupls));
        for k = 1:size(d,1)
            if k ~= bestCell
                isDuplicate{d(k,1)}(d(k,2)) = 1;
            end
            otherCells = d(setdiff(1:size(d,1),k),:);
            duplicates(d(k,1)).ROI(d(k,2)).plane = otherCells(:,1);
            duplicates(d(k,1)).ROI(d(k,2)).ID = otherCells(:,2);
            ind = ~isnan(dist(:,k));
            duplicates(d(k,1)).ROI(d(k,2)).distance = dist(ind,k);
            duplicates(d(k,1)).ROI(d(k,2)).corrCoeff = correlations(ind,k);
        end
    end
end