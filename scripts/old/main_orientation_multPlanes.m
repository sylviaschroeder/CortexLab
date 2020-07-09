%% Load meta
% subject='M141002_SS026';
% expDate='2014-10-29';
% exp=4;
% planes=2:3;

% subject='M141003_SS027';
% expDate='2014-10-23';
% exp=2;
% planes=1:4;
% expDate='2014-10-30';
% exp=3;
% planes=2:4;

subject='M141007_SS029';
expDate='2014-11-12';
exp=4;
planes=2:3;

% subject='M150323_SS042';
% expDate = '2015-04-14';
info=ppbox.infoPopulate(subject, expDate, exp);

%% Plot tuning selectivity statistics
% selThreshold = 0.3;
% gratings.printTuningSelectivityStats(info, planes, selThreshold)

%% Plot orientation tuning of brain regions
% suffix = 'registered';
suffix = 'ROI';
% suffix = 'raw';
prefDirs = [];
prefOris = [];
for iPlane = planes
    filePath = fullfile(info.folderProcessed, ...
        sprintf('%s_plane%03d_%s', info.basename2p, iPlane, suffix));
    load(filePath, 'meta')
    traces = meta.F(:,strcmp(meta.ROI.CellClasses, 's'));
    
    [stimTimes, stimSequence, stimMatrix, frameTimes, samplingRate] = ...
        ssLocal.getStimulusResponseInfo(meta);
    
    [~, neuronDirections, neuronBoth] = ...
        gratings.getPreferredOrientations(traces, samplingRate, ...
        stimMatrix, stimSequence, 0, 0);
    prefDirs = [prefDirs; neuronDirections(neuronDirections(:,2)>0.3,1)];
    prefOris = [prefOris; neuronBoth(neuronBoth(:,2)>0.3,1)];

%     filePathMovie = fullfile(info.folderProcessed, [meta.chData(1).basename '_' suffix]);
%     [data, infoReg] = loadArr(filePathMovie);
%     dataPatches = ssLocal.averagePixelsPerPatch(data, 5);
%     gratings.plotTuningOfBrainRegions(dataPatches, samplingRate, stimMatrix, stimSequence, iPlane);
end
figure
ind = prefDirs(:,1) < 11.25;
prefDirs(ind,1) = prefDirs(ind,1)+360;
hist(prefDirs(:,1), 22.5:22.5:360)
set(gca,'XTick',45:45:360)
xlabel('Preferred direction')
ylabel('# Neurons')
figure
ind = prefOris(:,1) < 11.25;
prefOris(ind,1) = prefOris(ind,1)+180;
hist(prefOris(:,1), 22.5:22.5:180)
set(gca,'XTick',45:45:180)
xlabel('Preferred orientation')
ylabel('# Neurons')

% publish('main_orientation_multPlanes.m','outputDir','\Results_4_tuningSelectivity', 'maxOutputLines', 0, 'figureSnapMethod', 'getframe', 'useNewFigure', false)