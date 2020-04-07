function RFs = getRFsOfBrainRegions(data, frameTimes, stimTimes, ...
    stimFrames, stimFrameTimes, numDivs)

patchSizeX = floor(size(data, 2) / numDivs);
patchSizeY = floor(size(data, 1) / numDivs);
totalX = patchSizeX * numDivs;
totalY = patchSizeY * numDivs;
% arrange all patches of each frame into one column
%  1st: arrange columns of patches in a 3D stack (columns stacked in 3rd
%  dimension, time is now in 4th dimension)
numFrames = size(data, 3);
d = reshape(data(1 : totalY, 1 : totalX, :), totalY, patchSizeX, ...
    numDivs, numFrames);
% clear data
%  2nd: flip 1st and 2nd dimension so that patches lign up in rows (all
%  pixels in d(:,a,b,c) belong to the same patch
d = permute(d, [2 1 3 4]);
%  3rd: arrange all pixels of each patch into one column, and all patches
%  of one frame into on plane (spanned by 1st and 2nd dimension)
d = reshape(d, patchSizeY * patchSizeX, numDivs^2, numFrames);
% average across all pixels in each frame patch
d = mean(d, 1);
% rearrange patch average values to original position within frame
d = reshape(squeeze(d), numDivs, numDivs, size(d, 3));

%% Find RFs for each patch in recorded frames (neuronal data)
fullRFs = cell(numDivs, numDivs);
for x = 1:numDivs
    for y = 1:numDivs
        [fullRFs(y,x), RFtimes] = whiteNoise.getReceptiveField( ...
            squeeze(d(y,x,:)), frameTimes, stimFrames, stimFrameTimes, ...
            stimTimes, {'Absolute'}, 0);
        fullRFs{y,x} = smooth3(fullRFs{y,x}, 'gaussian');
    end
end

%% Find spatial RF for each patch in recorded frames
RFs = cell(numDivs, numDivs);
maxi = 0;
for x = 1:numDivs
    for y = 1:numDivs
        % find time point with maximum values (average across all pixels)
        rf = reshape(fullRFs{y,x}, [], size(fullRFs{y,x}, 3));
        [m, ind] = max(mean(abs(rf), 1));
        if m > maxi
            maxi = m;
        end
        % average all frames within +-200 ms around max frame
        timeStep = diff(RFtimes(1:2));
        numFr = floor(0.2 / timeStep);
        RFs{y,x} = mean(fullRFs{y,x}(:,:,max(1,ind-numFr) : ...
            min(ind+numFr,size(fullRFs{y,x},3))), 3);
    end
end