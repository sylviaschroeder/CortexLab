function playReceptiveField(receptiveFields, RFTimes)

rf = cat(4, receptiveFields{:});
rf = mean(rf, 4);
clim = max(abs(rf(:)));
figure
for k = 1:size(rf, 3)
    imagesc(rf(:,:,k), [-clim clim])
    colorbar
    title([num2str(RFTimes(k)) ' s from frame onset'])
    pause(0.2)
end