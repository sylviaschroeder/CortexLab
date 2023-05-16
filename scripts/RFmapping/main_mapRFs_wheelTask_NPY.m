%% Folders
folderBase = 'C:\Users\Sylvia\OneDrive - University of Sussex\Lab\DATA\NPY\task_2p';
folderTools = 'C:\dev\toolboxes';
folderThisRepo = 'C:\dev\workspaces\CortexLab';
folderResults = 'C:\Users\Sylvia\OneDrive - University of Sussex\Lab\RESULTS\wheelTask\RFs_new';

%% Parameters
% for correcting baseline drifts of calcium traces at start of experiments
driftWin = 20; % in s, window to test whether baseline is higher than normal
driftThresh = 1.5; % in std, threshold for drift
correctWin = 150; % in s, window to fit exponential

% for receptive field estimates
% used for fitting 2 RFs (ON and OFF simultaneously), and fitting running
% kernels and RFs simultaneously
lambdasStim = logspace(-5, 0, 6);
RFlimits = [0.2 0.4];
crossFolds = 10;

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(fullfile(folderThisRepo))

%% Fit RFs and get cross-validated explained variance

subjects = dir(fullfile(folderBase, 'SS*'));
for subj = 3%1:length(subjects)
    name = subjects(subj).name;
    dates = dir(fullfile(folderBase, name, '2*'));
    for dt = 4%1:length(dates)
        date = dates(dt).name;
        folder = fullfile(folderBase, name, date);

        if ~isfile(fullfile(folder, '_ss_sparseNoise.times.npy'))
            continue
        end
        
        fResults = fullfile(folderResults, name, date);
        if ~isfolder(fResults)
            mkdir(fResults)
        end
        
        % load data
        traces = readNPY(fullfile(folder, '_ss_2pCalcium.dff.npy'));
        time = readNPY(fullfile(folder, '_ss_2pCalcium.timestamps.npy'));
        delays = readNPY(fullfile(folder, '_ss_2pPlanes.delay.npy'));
        planes = readNPY(fullfile(folder, '_ss_2pRois._ss_2pPlanes.npy'));
        stimTimes = readNPY(fullfile(folder, '_ss_sparseNoise.times.npy'));
        stimMaps = readNPY(fullfile(folder, '_ss_sparseNoiseID.map.npy'));
        stimSeq = readNPY(fullfile(folder, '_ss_sparseNoise._ss_sparseNoiseID.npy'));
        stimPos = readNPY(fullfile(folder, '_ss_sparseNoiseArea.edges.npy'));
        
        % interpolate calcium traces to align all to same time
        t_ind = time > stimTimes(1) - 10 & time < stimTimes(end) + 10;
        traces = traces(t_ind,:);
        time = time(t_ind);
        timeBin = median(diff(time));
        delta_t = median(diff(delays));
        upsample = round(timeBin / delta_t);
        timeBin = timeBin / upsample;
        t = reshape((time + (0:upsample-1) * timeBin)', [], 1);
        tr_up = NaN(length(t), size(traces,2));
        for d = 1:length(delays)
            indUnits = find(planes == d);
            for n = indUnits'
                if all(isnan(traces(:,n)))
                    continue
                end
                nanInd1 = isnan(traces(:,n));
                tr_up(:,n) = interp1(time(~nanInd1) + delays(d), ...
                    traces(~nanInd1,n), t, 'pchip');
                nanInd2 = reshape(repmat(nanInd1, 1, upsample), [], 1);
                tr_up(nanInd2,n) = NaN;
            end
        end
        
        % remove strong baseline decay at start of experiment in cells that
        % show it
        indUnits = find(nanmean(tr_up(1:round(driftWin / timeBin),:),1) > ...
            nanmean(tr_up,1) + driftThresh .* nanstd(tr_up,0,1));
        ind = round(correctWin / timeBin);
        for iUnit = 1:length(indUnits)
            y = tr_up(:, indUnits(iUnit));
            y = fillmissing(y, 'linear');
            % fit double exponential to start of trace
            f = fit((1:length(y))', y, ...
                @(a,b,c,d,e,x) a + b .* exp(-x ./ c) + d .* exp(-x ./ e), ...
                'Lower', [0 0 0 0 0], ...
                'Upper', [max(y) max(y) 500 max(y) 500], ...
                'StartPoint', [min(y) mean(y) 50 mean(y) 5]);
            % remove fit
            tr_up(:, indUnits(iUnit)) = y - f(1 : size(tr_up,1)) + f.a;
        end
        
        % more stimulus information
        stimFrames = stimMaps(stimSeq,:,:);
        stimFrameDur = median(diff(stimTimes));
        RFtimesInFrames = floor(RFlimits(1) / stimFrameDur) : ...
            ceil(RFlimits(2) / stimFrameDur);
        clear stimMaps stimSeq traces

        % map RF
        [rFields, ev, pred, t_pred] = ...
            whiteNoise.getReceptiveField( ...
            tr_up, t, stimFrames, stimTimes, ...
            RFtimesInFrames, lambdasStim, crossFolds);
        
        v = mean(ev,3); % [neuron x lamStim], average across cross-folds
        [maxEV, maxStimLam] = max(v,[],2);
        
        % CONTINUE HERE

        % test signficance of each RF
        [ev, ev_shift] = ...
            whiteNoise.receptiveFieldShiftTest( ...
            traces, time, stimFrames, stimTimes, ...
            RFtimesInFrames, runSpeed, runTime, ...
            runKernels, runWin, rFields, maxStimLam, 500);
        pvals = sum(ev_shift > ev, 2) ./ size(ev_shift,2);
        pvals(isnan(ev)) = NaN;

        % FIT ELLIPSE
        
        writeNPY(permute(rFields, [5 1 2 3 4]), fullfile(fResults, '_ss_rf.maps.npy'))
        writeNPY(maxEV, fullfile(fResults, '_ss_rf.explVars.npy'))
        writeNPY(maxEVRun, fullfile(fResults, '_ss_rf.explVarsRunning.npy'))
        writeNPY(maxEVStim, fullfile(fResults, '_ss_rf.explVarsStim.npy'))
        writeNPY(maxRunLam, fullfile(fResults, '_ss_rf.lambdasRunning.npy'))
        writeNPY(maxStimLam, fullfile(fResults, '_ss_rf.lambdasStim.npy'))
        writeNPY(rfPos, fullfile(fResults, '_ss_rfDescr.edges.npy'))
        writeNPY(RFtimesInFrames * stimFrameDur, fullfile(fResults, '_ss_rfDescr.timestamps.npy'))
        writeNPY(runKernels, fullfile(fResults, '_ss_rfRunningKernels.dff.npy'))
        writeNPY(runWin', fullfile(fResults, '_ss_rfRunningKernels.timestamps.npy'))
        
        writeNPY(pvals, fullfile(fResults, '_ss_rf.pValues.npy'))
    end
end