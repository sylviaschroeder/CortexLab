%% Load database
db_sparseNoise
% db_boutons_sparseNoise

%% Folders
% loading data
folderROIData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\InfoStructs';
% saving data
folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\receptiveFields\SC neurons';
% folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\boutons\receptiveFields';

%% Parameters
lambda = 0.08; %[0.02 0.04 0.08 0.16 0.3 0.6 1.2 2.4];
RFlimits = [0.4 1.1];
crossFolds = 10;

%% Fit RFs

% TODO: remove duplicate neurons!!!

RFs = db;
gDev = gpuDevice;
reset(gDev);
for k=1:length(db)
    fprintf('Dataset %d: %s %s exp.: %d\n', k, db(k).subject, ...
        db(k).date, db(k).exp);
    folder = fullfile(folderROIData, db(k).subject, ...
        db(k).date, num2str(db(k).exp));
    file = [sprintf('%s_%d_%s',db(k).date, db(k).exp, ...
        db(k).subject) '_2P_plane%03d_ROI.mat'];
    
    for iPlane = 1:length(db(k).planes)
        % load meta
        data=load(fullfile(folder, sprintf(file,db(k).planes(iPlane))));
        meta = data.meta;
        meta.folderTL = strrep(meta.folderTL, 'ioo.ucl.ac.uk', ...
            'cortexlab.net');
        
        % get stimulus information
        if iPlane == 1
            stimTimes = ppbox.getStimTimes(meta);
            [stimFrames, stimPosition] = whiteNoise.getStimulusFrames(meta);
            stimFrameDur = mean(stimTimes.offset - stimTimes.onset) / size(stimFrames,3);
            stimFrameTimes = ((1:size(stimFrames,3))-1) .* stimFrameDur;
            allStimFrameTimes = reshape((stimTimes.onset + stimFrameTimes)', [], 1);
            RFs(k).stimPosition = stimPosition;
        end
        
        %get calcium traces
        traceTimes = ppbox.getFrameTimes(meta);
        traces = meta.F_final;
        % set data of duplicate neurons or unhealthy neurons to NaN
        traces(:, meta.ROI.isDuplicate==1 | meta.ROI.isSwitchOn == 1) = NaN;
        
        stimTimeStep = median(diff(stimFrameTimes));
        RFtimesInFrames = floor(RFlimits(1) / stimTimeStep) : ...
            ceil(RFlimits(2) / stimTimeStep);
        RFTimes = RFtimesInFrames * stimTimeStep;
        
        fprintf('  Plane %d of %d (%d cells):', iPlane, length(db(k).planes), ...
            size(calciumTraces,2))
        
%         RFs(k).times = RFTimes;
        for iCell = 1:size(calciumTraces,2)
            if all(isnan(traces(:,iCell)))
                continue
            end
            fprintf(' %d', iCell)
            [receptiveFields, ev, traceSnippets, predictions] = ...
                whiteNoise.getReceptiveField( ...
                traces(:,iCell), traceTimes, stimFrames, stimFrameTimes, ...
                stimTimes, RFtimesInFrames, lambda, crossFolds, 0);
            RFs(k).plane(iPlane).cell(iCell).receptiveField = squeeze(mean(receptiveFields,2));
            RFs(k).plane(iPlane).cell(iCell).explainedVariance = ev;
            RFs(k).plane(iPlane).cell(iCell).traces = traceSnippets;
            RFs(k).plane(iPlane).cell(iCell).predictions = predictions;
            save(fullfile(folderResults, 'receptiveFields.mat'), 'RFs', ...
                'lambda', 'stimPosition')
        end
        fprintf('\n')
    end
end











        
figFolder = fullfile(folderResults, [db(k).subject '_' db(k).date ...
    '_' num2str(db(k).exp)], num2str(db(k).planes(iPlane)));
if ~exist(figFolder, 'dir')
    mkdir(figFolder);
end

smoothedRFs = cell(size(receptiveFields));
for type = 1:length(receptiveFields)
    % smooth receptive field
    smoothedRFs{type} = smooth3(receptiveFields{type}, 'gaussian');
end

RFinds = [];
for type = 1:length(RFtypes)
    if max(abs(smoothedRFs{type}(:))) / std(smoothedRFs{type}(:)) > 5
        infRF(type) = whiteNoise.fitGaussianToRF(smoothedRFs{type}, ...
            RFtimes, stimPosition, 0, RFtypes{type}, 0);
        RFinds = [RFinds, type];
    end
end
infRF = infRF(RFinds);

whiteNoise.plotReceptiveFieldAndFit(smoothedRFs, RFtimes, ...
    infRF, stimPosition, RFtypes);

fig = gcf;
fig.PaperPositionMode = 'auto';
print(fullfile(figFolder, sprintf('ReceptiveField%03d.tiff', iCell)), ...
    '-dtiff','-r0')
close gcf

RFinfo(end+1).RFs = infRF;


whiteNoise.plotAllRFs(RFinfo, 'Absolute');
fig = gcf;
fig.PaperPositionMode = 'auto';
print(fullfile(fullfile(folderResults, [db(k).subject '_' db(k).date ...
    '_' num2str(db(k).exp)]), 'AllReceptiveFields.tiff'), '-dtiff','-r0')
close gcf