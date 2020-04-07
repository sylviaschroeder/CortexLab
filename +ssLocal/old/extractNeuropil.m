function neuropilTraces = extractNeuropil(frames, cellMaps, pixelSize, ...
    neuropilInnerRing, neuropilOuterRing)

bufferFromROIs = neuropilInnerRing / pixelSize;
neuropilOutEdgeDistance = neuropilOuterRing / pixelSize;
neuropilTraces = zeros(size(frames, 3), length(cellMaps));

maskImage = zeros(size(frames, 1), size(frames, 2));
frames = reshape(frames, [], size(frames, 3));

for k = 1:length(cellMaps)
    maskImage(cellMaps{k}) = k;
end

props = regionprops(maskImage, 'Centroid', 'ConvexHull');
neuropilIndices = cell(1, length(cellMaps));
ROIbufferIndices = cell(1, length(cellMaps));
for roi = 1:length(props)
    radials = props(roi).ConvexHull - ...
        repmat(props(roi).Centroid, size(props(roi).ConvexHull, 1), 1);
    radLengths = sqrt(sum(radials.^2, 2));
    neuropilHull = repmat(props(roi).Centroid, size(props(roi).ConvexHull, 1), 1) + ...
        radials .* repmat((radLengths + neuropilOutEdgeDistance) ./ ...
        radLengths, 1, 2);
    mask = roipoly(maskImage, neuropilHull(:,1), neuropilHull(:,2));
    neuropilIndices{roi} = find(mask);
    
    bufferHull = repmat(props(roi).Centroid, size(props(roi).ConvexHull, 1), 1) + ...
        radials .* repmat((radLengths + bufferFromROIs) ./ ...
        radLengths, 1, 2);
    mask = roipoly(maskImage, bufferHull(:,1), bufferHull(:,2));
    ROIbufferIndices{roi} = find(mask);
end
ROIbufferIndices = unique(cat(1, ROIbufferIndices{:}));

for roi = 1:length(props)
    neuropil = setdiff(neuropilIndices{roi}, ROIbufferIndices);
    neuropilTraces(:, roi) = mean(frames(neuropil, :), 1);
end