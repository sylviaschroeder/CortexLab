%% Load database
db_boutons_sparseNoise

%% Folders
% loading data
folderROIData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\InfoStructs';
% saving data
folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\boutons\receptiveFields';

%% Parameters
RFtypes = {'Absolute', 'ON', 'OFF'};

%% Fit RFs and plot data

for k=1:length(db)
    fprintf('Dataset %d: %s %s exp.: %d\n', k, db(k).subject, ...
        db(k).date, db(k).exp);
    folder = [folderROIData filesep db(k).subject filesep ...
        db(k).date filesep num2str(db(k).exp)];
    fileStart = [db(k).date '_' num2str(db(k).exp) '_' ...
        db(k).subject];
    file = [fileStart '_2P_plane%03d_ROI.mat'];
    
    RFinfo = struct([]);
    for iPlane = 1:length(db(k).planes)
        % load meta
        data=load(fullfile(folder, sprintf(file,db(k).planes(iPlane))));
        meta = data.meta;
        calciumTraces = meta.F_final;
        
        frameTimes = ppbox.getFrameTimes(meta);
        stimTimes = ppbox.getStimTimes(meta);
        [stimFrames, stimFrameTimes, stimPosition] = whiteNoise.getStimulusFrames(meta);
        
        figFolder = fullfile(folderResults, [db(k).subject '_' db(k).date ...
            '_' num2str(db(k).exp)], num2str(db(k).planes(iPlane)));
        if ~exist(figFolder, 'dir')
            mkdir(figFolder);
        end
        
        for iCell = 1:size(calciumTraces,2)
            [receptiveFields, RFtimes, trialTraces] = whiteNoise.getReceptiveField( ...
                calciumTraces(:,iCell), frameTimes, stimFrames, stimFrameTimes, ...
                stimTimes, RFtypes, 0, 'corrs');
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
        end
    end
    whiteNoise.plotAllRFs(RFinfo, 'Absolute');
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(fullfile(fullfile(folderResults, [db(k).subject '_' db(k).date ...
        '_' num2str(db(k).exp)]), 'AllReceptiveFields.tiff'), '-dtiff','-r0')
    close gcf
end