function mergedImage = mergeImageWithMask(image, mask, colors)
% image     mean image with size of x- and y-range, all entries in the
%           range [0 1]
% mask      output from spatial.plotCellMasks
% colors    output from spatial.plotCellMasks
maskRGB = zeros([size(mask), 3]);
maskRGB = reshape(maskRGB, numel(mask), 3);
for col = 1:max(mask(:))
    ind = find(mask == col);
    maskRGB(ind,:) = repmat(colors(col+1,:), length(ind), 1);
end

image = image-min(image(:));
image = image./max(image(:));
mergedImage = repmat(image, 1, 1, 3);
mergedImage = reshape(mergedImage, numel(image), 3);
mergedImage = mergedImage + 0.2 * maskRGB;
colored = mask > 0;
colored = colored(:);
mergedImage(colored,:) = mergedImage(colored,:) ./ 1.2;
mergedImage = reshape(mergedImage, [size(image) 3]);