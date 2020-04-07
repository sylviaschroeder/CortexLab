function RFs = plotReceptiveField(receptiveField, RFtimes, stimPosition, RFtype)

RFs = cell(size(receptiveField));
for type = 1:length(receptiveField)
    % smooth receptive field
    RFs{type} = smooth3(receptiveField{type}, 'gaussian');
end
clim = [-1 1] * max(abs(reshape(cat(3, RFs{:}), [], 1)));

for type = 1:length(receptiveField)
    rf = RFs{type};
    
    screenSize = get(0, 'ScreenSize');
    cols = 8;
    rows = ceil(length(RFtimes) / cols);
    cols = ceil(length(RFtimes) / rows);
    figHeight = 0.2 * screenSize(4) * rows;
    figure('Position', [10 screenSize(4)-figHeight-100 ...
        screenSize(3)-20 figHeight])
    
    changePos = 0;
    for frame = 1:length(RFtimes)
        subplot(rows, cols, frame)
        imagesc(stimPosition([1 2]), stimPosition([3 4]), rf(:,:,frame), clim)
        axis image
        title(sprintf('%.2f s before response', RFtimes(frame)))
        if frame == 1
            pos = get(gca, 'Position');
            colorbar('WestOutside', 'Position', [0.08 pos(2) 0.01 pos(4)])
            aspect = get(gca, 'PlotBoxAspectRatio');
            height = aspect(1)/aspect(2);
            if height < 1
                height = height + (1-height)/2;
            else
                height = 1;
            end
            if pos(2) + height*pos(4) > 0.89
                changePos = 1;
                newHeight = (0.89-pos(2)) / height;
                newWidth = pos(3) * newHeight / pos(4);
            end
        else
            axis off
        end
        if changePos == 1
            pos = get(gca, 'Position');
            set(gca, 'Position', [pos([1 2]) newWidth newHeight])
        end
    end
    
    annotation('textbox', [0 0.9 1 0.1], 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'LineStyle', 'none', ...
        'String', [RFtype{type} ' receptive field'])
end