% Initialize info structure of dataset
subject = 'M141002_SS026';
expDate = '2014-10-29';
exp = 4;
info=ppbox.infoPopulate(subject, expDate, exp);
planes = 2:3;

for plane = planes
    % Load data ('meta')
    filePath = fullfile(info.folderProcessed, sprintf('%s_plane%03d_ROI', info.basename2p, plane));
    load(filePath, 'meta')
    
    % Plot cell masks
%     ssLocal.plotCellMasks(meta, 0, 1);
    
    % Get necessary data
    [stimTimes, stimSequence, stimMatrix, frameTimes, samplingRate] = ssLocal.getStimulusResponseInfo(meta);
    
    % Only consider neurons
    ind = strcmp(meta.ROI.CellClasses, 's');
    
    % Plot single trial and median responses to all stimuli
%     ssLocal.plotStimulusResponses(meta.F, samplingRate, stimMatrix, stimSequence, round(1 * samplingRate), round(4 * samplingRate));
    
    % Plot orientation tuning curves
    ssLocal.plotOrientationTuning(meta.F(:,ind), samplingRate, stimMatrix, stimSequence, 'single')
    
    % Plot orientation and direction map
    ssLocal.plotOrientationMap(meta.F(:,ind), samplingRate, stimMatrix, ...
        stimSequence, meta.ROI.CellMaps(ind), [length(meta.validY), length(meta.validX)]);

%     doPlot = 1;
%     if plane ~= planes(1)
%         doPlot = 0;
%     end
    % Plot orientation tuning in condition running vs. not running
%     ssLocal.plotConditionedOrientationTuning(meta, 'running', doPlot);
    % Plot orientation tuning in condition pupil dilated vs. constricted
%     ssLocal.plotConditionedOrientationTuning(meta, 'pupilSize', doPlot);
    % Plot orientation tuning in condition eye in centre vs, eye moved away
%     ssLocal.plotConditionedOrientationTuning(meta, 'pupilPos', doPlot);
end

% publish('main_publishResults.m','outputDir','\Results_1_twoGratings', 'maxOutputLines', 0, 'figureSnapMethod', 'getframe', 'useNewFigure', false)