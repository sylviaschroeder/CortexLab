%% Load database
make_db_wheelTask_visualNoise

%% Folders
% loading data
folderROIData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\InfoStructs';
% saving data
folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\wheelTask\receptiveFields';

%% Parameters
RFtypes = {'Absolute', 'ON', 'OFF'};

%% Fit RFs and plot data

for k=1:length(db)
    fprintf('Dataset %d: %s %s\n', k, db(k).mouse_name, db(k).date);
    folder = fullfile(folderROIData, db(k).mouse_name, ...
        db(k).date, num2str(db(k).expts(db(k).noiseExp)));
    fileStart = [db(k).date '_' num2str(db(k).expts(db(k).noiseExp)) '_' ...
        db(k).mouse_name];
    file = [fileStart '_2P_plane%03d_ROI.mat'];
    
    RFinfo = struct([]);
    trace = [];
    for iPlane = 1:length(db(k).planesToProcess)
        % load meta
        data=load(fullfile(folder, sprintf(file,db(k).planesToProcess(iPlane))));
        meta = data.meta;
        if iPlane == 1
            trace = meta.F_final;
            frameTimes = ppbox.getFrameTimes(meta);
        else
            tr = meta.F_final;
            t = ppbox.getFrameTimes(meta);
            trace = [trace, interp1(t, tr, frameTimes, 'pchip')'];
        end
    end
    trace = (trace - mean(trace,1)) ./ std(trace, 0, 1);
    trace = mean(trace, 2);
            
    stimTimes = ppbox.getStimTimes(meta);
    [stimFrames, stimFrameTimes, stimPosition] = whiteNoise.getStimulusFrames(meta);
    
    
    figFolder = fullfile(folderResults, [db(k).mouse_name(9:end) '_' db(k).date ...
        '_' num2str(db(k).newExp)]);
    if ~exist(figFolder, 'dir')
        mkdir(figFolder);
    end
        
    [receptiveFields, RFtimes, trialTraces] = whiteNoise.getReceptiveField( ...
        trace, frameTimes, stimFrames, stimFrameTimes, ...
        stimTimes, RFtypes, 0, 'corrs');
    smoothedRFs = cell(size(receptiveFields));
    for type = 1:length(receptiveFields)
        % smooth receptive field
        smoothedRFs{type} = smooth3(receptiveFields{type}, 'gaussian');
    end
            
    RFinds = [];
    clear infRF
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
    print(fullfile(figFolder, 'pixelAverageRF.tiff'), '-dtiff','-r0')
    close gcf
    
    RFinfo = infRF;
    save(fullfile(figFolder, 'RFinfo.mat'), 'RFinfo')
end