function printTuningSelectivityStats(info, planes, selThreshold)

totalOri = 0;
totalDir = 0;
totalOriAndDir = 0;
ori = 0;
dir = 0;
oriAndDir = 0;
timeAfterOnset = 0;
timeAfterOffset = 0;

orientationSelectivities = [];
directionSelectivities = [];

for iPlane = planes
    filePath = fullfile(info.folderProcessed, ...
        sprintf('%s_plane%03d_ROI', info.basename2p, iPlane));
    load(filePath, 'meta')
    
    [~, stimSequence, stimMatrix, ~, samplingRate] = ...
        ssLocal.getStimulusResponseInfo(meta);
    
    traces = meta.F(:, strcmp(meta.ROI.CellClasses, 's'));
    [oriSel, dirSel] = gratings.getTuningSelectivity_raw(traces, samplingRate, ...
        stimMatrix, stimSequence, timeAfterOnset, timeAfterOffset);
    totalOri = totalOri + sum(~isnan(oriSel));
    totalDir = totalDir + sum(~isnan(dirSel));
    totalOriAndDir = totalOriAndDir + sum(~isnan(oriSel) & ~isnan(dirSel));
    ori = ori + sum(oriSel >= selThreshold);
    dir = dir + sum(dirSel >= selThreshold);
    oriAndDir = oriAndDir + sum(oriSel >= selThreshold & ...
        dirSel >= selThreshold);
    
    orientationSelectivities = [orientationSelectivities; oriSel];
    directionSelectivities = [directionSelectivities; dirSel];
end

fprintf(['Orientation selective: %d of %d\n' ...
    'Direction selective: %d of %d\nOrientation and direction selective: %d of %d\n'], ...
    ori, totalOri, dir, totalDir, oriAndDir, totalOriAndDir)

figure
hist(orientationSelectivities)
xlabel('OSI')
ylabel('# Neurons')

figure
hist(directionSelectivities)
xlabel('DSI')
ylabel('# Neurons')