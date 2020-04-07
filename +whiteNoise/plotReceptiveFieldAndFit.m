function plotReceptiveFieldAndFit(RFs, RFtimes, fitRFs, ...
    stimPosition, RFtypes)

clim = [-1 1] * max(abs(reshape(cat(3, RFs{:}), [], 1)));

figure('Position', [1 41 1920 1083])

for type = 1:length(RFs)
    rf = RFs{type};
    for frame = 1:length(RFtimes)
        subplot(length(RFs)+3, length(RFtimes), ...
            (type-1)*length(RFtimes)+frame)
        imagesc(stimPosition([1 2]), stimPosition([3 4]), rf(:,:,frame), clim)
        axis image
        title(sprintf('-%.2f s', RFtimes(frame)), 'FontWeight', 'normal')
        if frame == 1
            if type == 1
                pos = get(gca, 'Position');
                colorbar('WestOutside', 'Position', [0.08 pos(2) 0.01 pos(4)])
            end
            ylabel([RFtypes{type} ' RF'], 'FontWeight', 'bold')
        else
            set(gca, 'XTick', [], 'YTick', [])
        end
    end
end

if isempty(fitRFs)
    return
end

[x, y] = meshgrid(0.5:0.1:size(RFs{1},2)+0.5, 0.5:0.1:size(RFs{1},1)+0.5);
xdata = zeros(size(x,1), size(x,2), 2);
xdata(:,:,1) = x;
xdata(:,:,2) = y;

subplot(length(RFs)+3, 1, length(RFs)+(1:3))
hold on
for field = 1:length(fitRFs)
    fitRF = D2GaussFunctionRot(fitRFs(field).gaussianFit, xdata);
    line = '-';
    switch fitRFs(field).RFtype
        case 'ON'
            color = 'r';
            if fitRFs(field).RFsign == -1
                line = '--';
            end
        case 'OFF'
            color = 'b';
            if fitRFs(field).RFsign == 1
                line = '--';
            end
        case 'Absolute'
            color = 'k';
            if fitRFs(field).RFsign == -1
                line = '--';
            end
        otherwise
            color = 'm';
    end
    contour(linspace(stimPosition(1), stimPosition(2), size(x, 2)), ...
        linspace(stimPosition(3), stimPosition(4), size(x, 1)), ...
        fitRF, [1 1] / 2 * fitRFs(field).gaussianFit(1), [color line], ...
        'LineWidth', 2.0);
end
legend({fitRFs.RFtype})
axis image
set(gcf, 'renderer', 'zbuffer')
set(gca,'YDir', 'reverse');
xlabel('screen X position');
ylabel('screen Y position');
title(sprintf('Position: %.2f, %.2f', fitRFs(1).RF_pos(1), fitRFs(1).RF_pos(2)))