%% Load database
db_sparseNoise

%% Folders
folderROIData = 'C:\Temp2p';

%% Parameters
% how to determine RFs
method='corrs'; % 'corrs' or 'kernels'
% to calculate F_0
window=400;
percentile=5;
% types of RFs
RFtypes = {'Absolute', 'ON', 'OFF'};

%% Loop across datasets
peaks = struct('raw', cell(1,length(db)), 'smoothed', cell(1,length(db)), ...
    'plane', cell(1,length(db)));
receptiveFields = struct('raw', cell(1,length(db)), ...
    'smoothed', cell(1,length(db)), 'plane', cell(1,length(db)));
RFinfo = struct('info', cell(1,length(db)), 'plane', cell(1,length(db)));
for k=1:length(db)
    
    %% Load data
    folder = [folderROIData filesep db(k).subject filesep ...
        db(k).date filesep num2str(db(k).exp)];
    fileStart = [db(k).date '_' num2str(db(k).exp) '_' ...
        db(k).subject];
    file = [fileStart '_2P_plane%03d_ROI.mat'];
    raw = [];
    smoothed = [];
    RFs_raw  = {};
    RFs_smoothed = {};
    planes = [];
    RFfits = struct([]);
    for iPlane=1:length(db(k).planes)
        % load meta
        data=load(fullfile(folder, sprintf(file,db(k).planes(iPlane))));
        info = orderfields(data.meta);
        neurons = find(strcmp(data.meta.ROI.CellClasses, 's'));
    
        %% Preprocess calcium traces
        % interpolate neural responses to get a single series of time points
        calciumTraces = info.F(:,neurons);
        calciumTime = ppbox.getFrameTimes(info);
        sampleInt = median(diff(calciumTime));
        % caclulate F_0
        [~,F_0] = ssLocal.removeSlowDrift(calciumTraces, calciumTime, window, percentile);
        calciumTraces = (calciumTraces - F_0) ./ F_0;
        
        %% Stimulus information
        stimTimes = ssLocal.getStimTimes(info);
        [stimFrames, stimFrameTimes, stimPosition] = whiteNoise.getStimulusFrames(info);
        
        %% Get RFs
        for iCell = 1:length(neurons)
            [recFields, RFtimes] = whiteNoise.getReceptiveField( ...
                calciumTraces(:, iCell), calciumTime, stimFrames, stimFrameTimes, ...
                stimTimes, RFtypes, 0, method);
            rfSm = cell(1, length(RFtypes));
            r = NaN(1, length(RFtypes));
            s = NaN(1, length(RFtypes));
            [r(1),indPeakRaw] = max(abs(recFields{1}(:)));
            r(1) = r(1) / std(recFields{1}(:));
            rfSm{1} = smooth3(recFields{1}, 'gaussian');
            [s(1),indPeakSmoothed] = max(abs(rfSm{1}(:)));
            s(1) = s(1) / std(rfSm{1}(:));
            for type = 2:length(recFields)
                % smooth receptive field
                rfSm{type} = smooth3(recFields{type}, 'gaussian');
                % get peaks
                r(type) = abs(recFields{type}(indPeakRaw)) / ...
                    std(recFields{type}(:));
                s(type) = abs(rfSm{type}(indPeakSmoothed)) / std(rfSm{type}(:));
            end
            raw = [raw; r];
            smoothed = [smoothed; s];
            RFs_raw = [RFs_raw; recFields];
            RFs_smoothed = [RFs_smoothed; rfSm];
            planes = [planes; db(k).planes(iPlane)];
            
            if s(1) > 5
                RFfits(end+1).RFs = whiteNoise.fitGaussianToRF(recFields{1}, ...
                    RFtimes, stimPosition, 0, RFtypes{1}, 0);
            else
                RFfits(end+1).RFs = [];
            end
        end
    end
    peaks(k).raw.Absolute = raw(:,1);
    peaks(k).raw.On = raw(:,2);
    peaks(k).raw.Off = raw(:,3);
    peaks(k).smoothed.Absolute = smoothed(:,1);
    peaks(k).smoothed.On = smoothed(:,2);
    peaks(k).smoothed.Off = smoothed(:,3);
    peaks(k).plane = planes;
    receptiveFields(k).raw.Absolute = RFs_raw(:,1);
    receptiveFields(k).raw.On = RFs_raw(:,2);
    receptiveFields(k).raw.Off = RFs_raw(:,3);
    receptiveFields(k).smoothed.Absolute = RFs_smoothed(:,1);
    receptiveFields(k).smoothed.On = RFs_smoothed(:,2);
    receptiveFields(k).smoothed.Off = RFs_smoothed(:,3);
    receptiveFields(k).plane = planes;
    RFinfo(k).info = RFfits;
    RFinfo(k).plane = planes;
end