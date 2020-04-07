% this script will play the .bin file and calculate the average frame

[filename, filepath] = uigetfile('*.bin', [], 'G:\Processing\M140701_MK012\2014-08-10');
[~, fn, fext] = fileparts(filename);
fullfilestem = fullfile(filepath, fn);
[sz, cl, info] = loadArrInfo(fullfilestem);
nFrames = sz(3)

clear allFrames;

% now deciding for how many frames we have enough memory
[user,sys] = memory;

fakeFrame = zeros(sz(1), sz(2), cl);
ffInfo = whos('fakeFrame');
bytesPerFrame = ffInfo.bytes;
nFrames2LoadMax = floor(sys.PhysicalMemory.Available/bytesPerFrame*0.9)

nFrames = min(nFrames, nFrames2LoadMax)

% loading all the stack (as much as possible)
allFrames = loadMovieFrames(fullfilestem, 1, nFrames);

maxIntensity = max(allFrames(:));
minIntensity = min(allFrames(:));

meanFrame = int16(mean(allFrames, 3));
saveastiff(meanFrame, [fullfilestem, '_AVG.tiff']);

% let's stretch the contrast slightly, for plotting only
% define the lowest 0.1% of the pixels to be black
low = prctile(meanFrame(:), 0.1);
% the highest 0.1% of the pixels will be white
high = prctile(meanFrame(:), 99.9);
meanFrame(meanFrame>high) = high;
meanFrame(meanFrame<low) = low;

figure('Position', [10  300  1840  760]);
subplot(1, 2, 2);
imagesc(meanFrame);
colormap gray;
axis equal tight off
title(strrep(filename, '_', '\_'));
subplot(1, 2, 1);

tic
for iFrame = 1:nFrames
    frame = double(allFrames(:,:,iFrame)-minIntensity)/double(maxIntensity-minIntensity)*64;
    if iFrame == 1
        h = image(frame);
        axis equal tight off
    else
        set(h, 'CData', frame);
    end
    if ~mod(iFrame, 10)
        title(sprintf('%d/%d', iFrame, nFrames))
    end
    drawnow;
%         pause(0.001);
end
toc

fprintf('Played %d frames at %3.0f fps\n', nFrames, nFrames/toc);


