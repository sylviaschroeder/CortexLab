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
% to estimate response amplitudes
minSR = 2; % minimum sample rate of calcium trace is minSR * sample rate of visual stimulus
responseTime = .7; % in s, estimate response amplitude within this time after frame onset
plotKernels = true;

lambda = [0.02 0.04 0.08 0.16 0.3 0.6 1.2 2.4];
RFlimits = [0.4 1.1];
crossFolds = 10;

%% Fit RFs

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
    
    respVectors = []; % [neurons x time]
    respAmplitudes = []; %[neurons x stimFrames]
    
    % for each neuron, get response amplitude for each frame
    for iPlane = 1:length(db(k).planes)
        if plotKernels
            fPlots = fullfile(folderResults, 'plots_kernels', ...
                sprintf('%s_%s', db(k).subject, db(k).date), ...
                sprintf('Plane%02d',db(k).planes(iPlane)));
            if ~isfolder(fPlots)
                mkdir(fPlots)
            end
        end
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
        sampleDur = median(diff(traceTimes));
        maxDur = stimFrameDur / minSR;
        if sampleDur > maxDur
            n = ceil(sampleDur / maxDur);
            upDur = sampleDur / n;
            upTimes = traceTimes(1):upDur:traceTimes(end);
        else
            upDur = sampleDur;
            upTimes = traceTimes;
        end
        traces = meta.F_final;
        % set data of duplicate neurons or unhealthy neurons to NaN
        traces(:, meta.ROI.isDuplicate==1 | meta.ROI.isSwitchOn == 1) = NaN;
        
        vectorLen = round(responseTime / upDur);
        vectors = NaN(vectorLen, size(traces,2));
        amplitudes = NaN(length(allStimFrameTimes), size(traces,2));
        predictions = NaN(size(traces,1), size(traces,2));
        adjR2 = NaN(1, size(traces,2));
        
        fprintf('  Plane %d of %d (%d cells)\n', iPlane, length(db(k).planes), ...
            size(traces,2))
        for iCell = 1:size(traces,2)
            if all(isnan(traces(:,iCell)))
                continue
            end
            fprintf(' %d', iCell)
            % upsample trace so that sampling rate is at least 10 times
            % that of stimulus rate
            nanOrig = isnan(traces(:,iCell));
            if sampleDur > maxDur
                upTrace = interp1(traceTimes(~nanOrig), ...
                    traces(~nanOrig,iCell), upTimes, 'spline');
                nanUp = general.findMatchingIndices(nanOrig, traceTimes, upTimes);
                upTrace(nanUp) = NaN;
            else
                upTrace = traces(:,iCell);
            end
            % get "pulse onsets", i.e. the samples corresponding to
            % stimulus frame onsets
            goodStimFrames = allStimFrameTimes < traceTimes(end) & ...
                allStimFrameTimes > traceTimes(1);
            pulses = find(hist(allStimFrameTimes(goodStimFrames), upTimes));
            % determine temporal response to all stim frames and amplitude
            % for each stim frame
            [amps, v, ~, pred, ~, r2] = models.SumOfPulses(upTrace, pulses, ...
                vectorLen, 0, 1);
            vectors(:,iCell) = v;
            amplitudes(goodStimFrames,iCell) = amps;
            predictions(:,iCell) = interp1(upTimes, pred, traceTimes, 'spline');
            adjR2(iCell) = r2;
            if plotKernels
                fig = gcf;
                fig.PaperPositionMode = 'auto';
                print(fullfile(fPlots, sprintf('Neuron%03d.jpg', iCell)), ...
                    '-djpeg','-r0')
                savefig(fig,fullfile(fPlots, sprintf('Neuron%03d.fig', iCell)), 'compact')
                close gcf
            end
        end
        fprintf('\n')
        RFs(k).plane(iPlane).tempResp = vectors;
        RFs(k).plane(iPlane).amplitudes = amplitudes;
        RFs(k).plane(iPlane).predictions = predictions;
        RFs(k).plane(iPlane).adjR2 = adjR2;
        respVectors = [respVectors, vectors];
        respAmplitudes = [respAmplitudes, amplitudes];
    end
    RFs(k).times = (0:vectorLen-1) .* upDur;
    % reshape stimulus frames to [frames x pixels] and repeat
    stimulus = repmat(reshape(stimFrames, [], size(stimFrames,3))', ...
        length(stimTimes.onset), 1);
    % exclude NaN responses
    goodNeurons = ~all(isnan(respAmplitudes),1);
    goodFrames = ~any(isnan(respAmplitudes(:,goodNeurons)), 2);
    rf = abs(stimulus(goodFrames,:))\respAmplitudes(goodFrames,goodNeurons);
    rf = permute(reshape(rf, size(stimFrames,1), size(stimFrames,2), []), [3 1 2]);
end
