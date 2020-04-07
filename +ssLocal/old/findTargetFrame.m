function [targetFrame, targetValue] = findTargetFrame(movie)
% Finds target frame within movie, which has minimum luminance and which is
% used for registration later
% targetFrame - (Y,X) array
% targetValue - mean luminance value of targetFrame
 
[h, w, nFrames] = size(movie);
meanF = smooth(mean(reshape(movie, h*w, nFrames)));
%now look in the middle third of the image frames for the minimum
fromFrame = round(nFrames*1/3);
toFrame = round(nFrames*2/3);
[targetValue, idx] = min(meanF(fromFrame:toFrame));
minFrame = fromFrame + idx;
%create a Gaussian filter for filtering registration frames
hGauss = fspecial('gaussian', [5 5], 1);
%Gaussian filter the target image
targetFrame = imfilter(movie(:,:,minFrame), hGauss, 'same', 'replicate');