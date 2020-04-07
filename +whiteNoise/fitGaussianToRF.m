function RFinfo = fitGaussianToRF(receptiveField, RFtimes, stimPosition, ...
    makePlots, RFtype, separate)

if separate
    [spatialField, timeCourse] = rfMapMakeSeparable(receptiveField, 0);
    if max(timeCourse) < -min(timeCourse)
        spatialField = -spatialField;
        timeCourse = -timeCourse;
    end
else
    [~, coord] = max(abs(receptiveField(:)));
    [maxRow, maxCol, maxTime] = ind2sub(size(receptiveField), coord);
    spatialField = receptiveField(:,:,maxTime);
    timeCourse = squeeze(receptiveField(maxRow, maxCol, :) / ...
        receptiveField(maxRow, maxCol, maxTime));
end
RFsign = 1;
color = 'r';
if max(spatialField(:)) < -1 * min(spatialField(:))
    RFsign = -1;
    color = 'b';
end

[gfit, fitResp] = whiteNoise.fit2dGaussRF(1:size(spatialField, 2), ...
    1:size(spatialField, 1), spatialField * RFsign, 0);

RFinfo.RF_raw = receptiveField;
RFinfo.RF_separated = spatialField;
RFinfo.RF_timecourse = timeCourse;
RFinfo.gaussianFit = gfit;
RFinfo.fitRF = fitResp;
RFinfo.RFsign = RFsign;
RFinfo.stimPosition = stimPosition;
RFinfo.RFtype = RFtype;
posGrid = [gfit(2), gfit(4)]; % [x y]
posVisualDegrees = [ ...
    stimPosition(1)+(posGrid(1)-0.5)*diff(stimPosition(1:2))/size(spatialField,2), ...
    stimPosition(3)+(posGrid(2)-0.5)*diff(stimPosition(3:4))/size(spatialField,1) ];
RFinfo.RF_pos = posVisualDegrees;

if makePlots
    figure('Position', [680 230 960 870])
    subplot(2,2,1)
    scaleMax = max(abs(receptiveField(:)));
    for pixRow = 1:size(receptiveField,1)
        for pixCol = 1:size(receptiveField,2)
            
            thisTrace = squeeze(receptiveField(pixRow, pixCol, :))./scaleMax;
            plot(pixCol+(1:size(receptiveField,3))/ ...
                (size(receptiveField,3)*1.1), ...
                size(receptiveField,1)-pixRow+thisTrace, 'k');
            hold on;
            
            thisTracePos = thisTrace;
            thisTracePos(thisTrace<1/3) = NaN;
            plot(pixCol+(1:size(receptiveField,3))/ ...
                (size(receptiveField,3)*1.1), ...
                size(receptiveField,1)-pixRow+thisTracePos, 'r', ...
                'LineWidth', 2.0);
            
            thisTraceNeg = thisTrace;
            thisTraceNeg(thisTrace>-1/3) = NaN;
            plot(pixCol+(1:size(receptiveField,3))/ ...
                (size(receptiveField,3)*1.1), ...
                size(receptiveField,1)-pixRow+thisTraceNeg, 'b', ...
                'LineWidth', 2.0);
            
        end
    end
    ylim([-0.1 size(receptiveField,1)])
    title('Raw Ca-correlated average');
    set(gca, 'XTick', [], 'YTick', []);
    axis tight;
    
    subplot(2,2,2);
    xPointDist = (stimPosition(2) - stimPosition(1)) / size(spatialField, 2);
    xPoints = stimPosition(1)+xPointDist/2 : xPointDist : stimPosition(2);
    yPointDist = (stimPosition(4) - stimPosition(3)) / size(spatialField, 1);
    yPoints = stimPosition(3)+yPointDist/2 : yPointDist : stimPosition(4);
    imagesc(xPoints([1 end]), yPoints([1 end]), ...
        spatialField, [-1 1] * max(abs(spatialField(:))));
%     colormap(colormap_blueblackred);
    title('Space-time separable RF')
    xlabel('Distance from vertical meridian (in degrees)')
    ylabel('Distance from horizontal meridian (in degrees)')
    title('Spatial RF (space-time separated)')
    
    subplot(2,2,3)
    plot(RFtimes, timeCourse, 'k')
    xlabel('Time (in s)')
    ylabel('Modulation')
    title('Temporal RF (space-time separated)')    
    
    subplot(2,2,4);
    whiteNoise.plotFitReceptiveField(size(spatialField), gfit, stimPosition, color);
    
    if nargin > 4
        annotation('textbox', [0 0.9 1 0.1], 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'LineStyle', 'none', ...
        'String', [RFtype ' receptive field'])
    end
end