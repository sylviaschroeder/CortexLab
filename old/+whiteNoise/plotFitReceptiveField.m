function plotFitReceptiveField(stimGridSize, fitRFparameters, stimPosition, ...
    RFtypes, RFsigns)

[x, y] = meshgrid(0.5:0.1:stimGridSize(2)+0.5, 0.5:0.1:stimGridSize(1)+0.5);
xdata = zeros(size(x,1), size(x,2), 2);
xdata(:,:,1) = x;
xdata(:,:,2) = y;

if ~iscell(fitRFparameters)
    fitRFparameters = {fitRFparameters};
end
if ~iscell(RFtypes)
    RFtypes = {RFtypes};
end
hold on
for field = 1:length(fitRFparameters)
    fitRF = D2GaussFunctionRot(fitRFparameters{field}, xdata);
    line = '-';
    switch RFtypes{field}
        case 'ON'
            color = 'r';
            if RFsigns(field) == -1;
                line = '--';
            end
        case 'OFF'
            color = 'b';
            if RFsigns(field) == 1;
                line = '--';
            end
        otherwise
            color = 'k';
    end
    contour(linspace(stimPosition(1), stimPosition(2), size(x, 2)), ...
        linspace(stimPosition(3), stimPosition(4), size(x, 1)), ...
        fitRF, [1 1] / 2 * fitRFparameters{field}(1), [color line], ...
        'LineWidth', 2.0);
end

if nargin > 4
    if ~iscell(RFtypes)
        RFtypes = {RFtypes};
    end
    legend(RFtypes)
end
set(gcf, 'renderer', 'zbuffer')
title('2-D gaussian fit contour')
set(gca,'YDir', 'reverse');
xlabel('screen X position');
ylabel('screen Y position');