function dataPatches = averagePixelsPerPatch(data, numDivs)

patchSizeX = floor(size(data, 2) / numDivs);
patchSizeY = floor(size(data, 1) / numDivs);
totalX = patchSizeX * numDivs;
totalY = patchSizeY * numDivs;
% arrange all patches of each frame into one column
%  1st: arrange columns of patches in a 3D stack (columns stacked in 3rd
%  dimension, time is now in 4th dimension)
numFrames = size(data, 3);
dataPatches = reshape(data(1 : totalY, 1 : totalX, :), totalY, patchSizeX, ...
    numDivs, numFrames);
% clear data
%  2nd: flip 1st and 2nd dimension so that patches lign up in rows (all
%  pixels in d(:,a,b,c) belong to the same patch
dataPatches = permute(dataPatches, [2 1 3 4]);
%  3rd: arrange all pixels of each patch into one column, and all patches
%  of one frame into on plane (spanned by 1st and 2nd dimension)
dataPatches = reshape(dataPatches, patchSizeY * patchSizeX, numDivs^2, numFrames);
% average across all pixels in each frame patch
dataPatches = mean(dataPatches, 1);
% rearrange patch average values to original position within frame
dataPatches = reshape(squeeze(dataPatches), numDivs, numDivs, size(dataPatches, 3));