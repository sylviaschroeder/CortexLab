function neuropilRatios = extractNeuropilRatio(frames, pixelSize, ...
    neuropilInnerRing, neuropilOuterRing)

meanFrame = mean(frames, 3);
vesselMaps = {};
figure
imagesc(meanFrame)
axis image
colormap gray
hold on
while true
    [vesselROI, x, y] = roipoly;
    fill(x, y, 'w', 'FaceAlpha', 0.4)
    vesselMaps{end+1} = find(vesselROI);
    display('Press ''q'' to quit or any other key to specify more ROIs.')
    press = waitforbuttonpress;
    if press == 1 && strcmp(get(gcf, 'CurrentCharacter'), 'q')
        break
    end
end

neuropilRatios = zeros(size(frames, 3), length(vesselMaps));
surroundTraces = ssLocal.extractNeuropil(frames, vesselMaps, pixelSize, ...
    neuropilInnerRing, neuropilOuterRing);
frames = reshape(frames, [], size(frames, 3));
for roi = 1:length(vesselMaps)
    trace = mean(frames(vesselMaps{roi},:), 1)';
    neuropilRatios(:, roi) = trace ./ surroundTraces(:, roi);
end