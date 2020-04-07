function saveFramesAsMovie(frames, framesToAverage, frameTimes, speedUp, ...
    fileName)

yPx = size(frames,1);
xPx = size(frames, 2);
numFrames = floor(size(frames, 3) / framesToAverage);
% cut off last few frames
frames(:,:, numFrames * framesToAverage + 1 : end) = [];
frames = reshape(frames, xPx * yPx, framesToAverage, numFrames);
frames = mean(frames, 2);
frames = reshape(frames, yPx, xPx, []);

samplingRate = 1 / median(diff(frameTimes)) / framesToAverage * speedUp;
mini = min(frames(:));
maxi = prctile(frames(:), 99.999);

figure
colormap gray
mov(1:numFrames) = struct('cdata', [], 'colormap', []);
for k = 1:numFrames
    imagesc(frames(:,:,k),[mini maxi])
    axis off image
    mov(k) = getframe(gcf);
end

movie2avi(mov, [fileName '.avi'], 'fps', samplingRate, 'compression', 'none')
close gcf