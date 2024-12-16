%% Folders
folderSource = 'Z:\ProcessedData';
folderTarget = 'C:\Users\Sylvia\University College London\Rossi, Federico - SC_Saccades_paper\Paper\SC_data';

%% Loop over datasets
db =  db_saccadeDataFromSchroederLab;
for k = 1:length(db)
    fs = fullfile(folderSource, db(k).subject, db(k).date);
    ft = fullfile(folderTarget, db(k).subject, db(k).date, '001');
    if ~isfolder(ft)
        mkdir(ft);
    end
    % copyfile(fullfile(fs, 'calcium.dff.npy'), ...
    %     fullfile(ft, '_ss_2pCalcium.dff.npy'));
    % copyfile(fullfile(fs, 'calcium.timestamps.npy'), ...
    %     fullfile(ft, '_ss_2pCalcium.timestamps.npy'));
    % copyfile(fullfile(fs, 'planes.delay.npy'), ...
    %     fullfile(ft, '_ss_2pPlanes.delay.npy'));
    % planes = readNPY(fullfile(fs, 'rois.planes.npy'));
    % writeNPY(planes + 1, fullfile(ft, '_ss_2pRois._ss_2pPlanes.npy'));
    % copyfile(fullfile(fs, 'rois.id.npy'), ...
    %     fullfile(ft, '_ss_2pRois.ids.npy'));
    % % NOTE: x and y position is given in pixels, not in microns!
    % copyfile(fullfile(fs, 'rois.xyz.npy'), ...
    %     fullfile(ft, '_ss_2pRois.xyz.npy'));
    % 
    % if isfile(fullfile(fs, 'sparse.startTime.npy'))
    %     t = readNPY(fullfile(fs, 'sparse.startTime.npy'));
    %     writeNPY((1:length(t))', ...
    %         fullfile(ft, '_ss_sparseNoise._ss_sparseNoiseID.npy'))
    %     writeNPY(t, ...
    %         fullfile(ft, '_ss_sparseNoise.times.npy'));
    %     edges = readNPY(fullfile(fs, 'sparseExp.edges.npy'));
    %     edges = edges(end-3:end)'; % due to older bug where 1st entry in edges
    %     % contained number of squares, then edges
    %     % [bottom top left right]
    %     % translate edges to [left right top bottom] where values for top and
    %     % bottom are negative if above horizon and positive if below horizon
    %     writeNPY([edges([3 4]) -edges([2 1])], ...
    %         fullfile(ft, '_ss_sparseNoiseArea.edges.npy'));
    %     frames = readNPY(fullfile(fs, 'sparse.map.npy'));
    %     frames(frames==0) = -1;
    %     frames(frames==0.5) = 0;
    %     writeNPY(frames, fullfile(ft, '_ss_sparseNoiseID.map.npy'));
    % end
    % 
    % if isfile(fullfile(fs, 'circles.startTime.npy'))
    %     t1 = readNPY(fullfile(fs, 'circles.startTime.npy'));
    %     t2 = readNPY(fullfile(fs, 'circles.endTime.npy'));
    %     writeNPY([t1 t2], fullfile(ft, '_ss_circles.intervals.npy'));
    %     copyfile(fullfile(fs, 'circles.diameter.npy'), ...
    %         fullfile(ft, '_ss_circles.diameter.npy'));
    %     copyfile(fullfile(fs, 'circles.isWhite.npy'), ...
    %         fullfile(ft, '_ss_circles.isWhite.npy'));
    %     x = readNPY(fullfile(fs, 'circles.x.npy'));
    %     y = readNPY(fullfile(fs, 'circles.y.npy'));
    %     writeNPY([x y], fullfile(ft, '_ss_circles.xyPos.npy'));
    %     copyfile(fullfile(fs, 'circlesExp.intervals.npy'), ...
    %         fullfile(ft, '_ss_recordings.circles_intervals.npy'));
    % end
    % 
    % copyfile(fullfile(fs, 'darkScreen.intervals.npy'), ...
    %     fullfile(ft, '_ss_recordings.grayScreen_intervals.npy'));

    % eye diameter and xy positions are provided in pixels (Raikhan
    % confirmed)
    t = readNPY(fullfile(fs, 'eye.timestamps.npy'));
    d = readNPY(fullfile(fs, 'eye.diameter.npy'));
    p = readNPY(fullfile(fs, 'eye.xyPos.npy'));
    % interpolate data to regular time stamps
    binSize = median(diff(t), "omitnan");
    t_new = (ceil(min(t)/binSize) : floor(max(t)/binSize))' .* binSize;
    % only use unique time stamps (needed for interpolation later)
    [~,indT] = unique(t);
    % ignore NaN values
    indT = indT(~isnan(t(indT)));
    t1 = t(indT);
    d = d(indT);
    p = p(indT,:);
    % interpolate diameters
    indVal = ~isnan(d);
    d_new = interp1(t1(indVal), d(indVal), t_new);
    indNaN = histcounts(t1(~indVal), t_new) > 0;
    d_new(indNaN) = NaN;
    % interpolate positions
    indVal = all(~isnan(p),2);
    x_new = interp1(t1(indVal), p(indVal,1), t_new);
    y_new = interp1(t1(indVal), p(indVal,2), t_new);
    p_new = [x_new, y_new];
    indNaN = histcounts(t1(~indVal), t_new) > 0;
    p_new(indNaN,:) = NaN;
    % set data to zero, if there are long gaps (>1 s) between original time 
    % stamps
    t1 = t(~isnan(t));
    gaps = find(diff(t1) > 1);
    indNaN = [];
    for g = gaps'
        j = find(t_new > t1(g) & t_new < t1(g+1));
        indNaN = [indNaN; j];
    end
    d_new(indNaN) = NaN;
    p_new(indNaN,:) = NaN;
    writeNPY(t_new, fullfile(ft, 'eye.timestamps.npy'));
    writeNPY(d_new, fullfile(ft, 'eye.diameter.npy'));
    writeNPY(p_new, fullfile(ft, 'eye.xyPos.npy'));
end