function distances = getNeuronDistances(ROImaps, mapSizePixels, mapSizeMM)


if nargin < 3
    mapSizeMM = mapSizePixels;
end

ROIcentres = NaN(length(ROImaps), 2);
for roi = 1:length(ROImaps)
    im = false(mapSizePixels);
    im(ROImaps{roi}) = true;
    stats = regionprops(im, 'Centroid');
    ROIcentres(roi,:) = stats.Centroid;
end
% convert pixels to microns
factors = mapSizeMM ./ mapSizePixels;
ROIcentres = bsxfun(@times, ROIcentres, factors);

distances = sqrt(bsxfun(@minus, ROIcentres(:,1), ROIcentres(:,1)') .^ 2 + ...
    bsxfun(@minus, ROIcentres(:,2), ROIcentres(:,2)') .^ 2);