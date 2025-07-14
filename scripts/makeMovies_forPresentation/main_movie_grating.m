folder = 'D:\';
framerate = 30;
duration = 5;
cyclesPerFrame = 5;
cyclesPerSec = 2;
numPix = 200;

numFrames = duration * framerate;
pixPerCycle = ceil(numPix / cyclesPerFrame);
pixShiftPerSec = pixPerCycle * cyclesPerSec;
pixShiftPerFrame = max(1, round(pixShiftPerSec / framerate));

x = linspace(1, 2*pi*(cyclesPerFrame+2), ...
    ceil(numPix / cyclesPerFrame * (cyclesPerFrame+2)));
y = sin(x);
A = repmat(y, numPix, 1);


v = VideoWriter(fullfile(folder, 'grating.mp4'), 'MPEG-4');
v.FrameRate = framerate;
open(v);

figure
axis image off ij
colormap gray
set(gca, 'nextplot', 'replacechildren')
for f = 1:numFrames
    shift = mod(round(pixShiftPerFrame*(f-1)), pixPerCycle);
    imagesc(A(:,(1:numPix) + shift))
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v)
close gcf