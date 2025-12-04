%% Important info
% Description of spike sorting outputs:
% https://docs.google.com/document/d/1OqIqqakPakHXRAwceYLwFY9gOrm8_P62XIfCTnHwstg/edit?tab=t.0#heading=h.5houj8bng5o

%% Folders
folderSourceRaw = 'Z:\RawData\';
folderSourceProcessed = 'Z:\ProcessedData';
% The Machine:
folderTarget = 'C:\Users\sylvi\OneDrive - University of Sussex\Projects\2023_OrientationColumns\DataToPublish\ephys';
% folderTarget = 'C:\Users\Sylvia\OneDrive - University of Sussex\Projects\2023_OrientationColumns\DataToPublish\ephys';
folderTools = 'C:\dev\toolboxes';
foldersRepo = 'C:\dev\workspaces\CortexLab';

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(genpath(fullfile(foldersRepo)))

%% Setup Python packages
np = py.importlib.import_module('numpy');

%% Loop over datasets
db =  db_orientationDataFromSchroederLab;
for k = 1:length(db)
    fprintf("%s %s\n", db(k).subject, db(k).date)
    ft = fullfile(folderTarget, db(k).subject, db(k).date);
    if ~isfolder(ft)
        mkdir(ft);
    end
    fs = fullfile(folderSourceProcessed, db(k).subject, db(k).date, 'exp01');

    % copy spike data
    fMeta = fullfile(folderSourceRaw, db(k).subject, db(k).date, ...
        'exp01', sprintf('%s_%s_g0', db(k).subject, db(k).date), ...
        sprintf('%s_%s_g0_imec0', db(k).subject, db(k).date));
    fCurated = fullfile(fMeta, 'ibl_sorter', 'curated');
    fAlf = fullfile(fMeta, 'ibl_sorter', 'alf');
    % From /alf folder: 
    % spike times, clusters, amplitudes, depths
    % template waveforms, waveform channels
    % channel coordinates
    % IBL documentation: https://docs.google.com/document/d/1OqIqqakPakHXRAwceYLwFY9gOrm8_P62XIfCTnHwstg/edit?tab=t.0#heading=h.5houj8bng5o
    spikeTimes = readNPY(fullfile(fAlf, 'spikes.times.npy'));
    % clusters: integers counting from 0
    spikeClusters = readNPY(fullfile(fAlf, 'spikes.clusters.npy'));
    % amplitude of each spike in V
    spikeAmplitudes = readNPY(fullfile(fAlf, 'spikes.amps.npy'));
    % Depth along probe of each spike (µm; computed from waveform center of
    % mass). 0 means probe tip, positive upwards.
    spikeDepths = readNPY(fullfile(fAlf, 'spikes.depths.npy')); 
    % Waveform from spike sorting templates (in V).
    clusterWaveforms = readNPY(fullfile(fAlf, 'clusters.waveforms.npy'));
    % Index of channels that are stored for each cluster waveform. Sorted 
    % by increasing distance from the maximum amplitude channel.
    clusterWfChannels = readNPY(fullfile(fAlf, 'clusters.waveformsChannels.npy'));
    % From /curated folder:
    % cluster_info
    % spike_clusters
    cluster_info = readtable(fullfile(fCurated, 'cluster_info.tsv'), ...
        'FileType', 'text');
    good = strcmp('good', cluster_info.group);
    clusters_good = cluster_info.cluster_id(good);
    spikeClusters_curated = readNPY(fullfile(fCurated, 'spike_clusters.npy'));
    spike_ind = ismember(spikeClusters_curated, clusters_good);
    % curate cluster waveforms and waveform channels:
    % check whether curated cluster was merged, choose waveform of original
    % cluster with highest number of spikes
    clusterWaveforms_curated = NaN([length(clusters_good), ...
        size(clusterWaveforms, [2 3])]);
    clusterWfChannels_curated = NaN(length(clusters_good), ...
        size(clusterWfChannels, 2));
    for cl = 1:length(clusters_good)
        ind = spikeClusters_curated == clusters_good(cl);
        clusters = unique(spikeClusters(ind));
        if length(clusters) > 1 % cluster was merged
            mx = 0;
            mxC = NaN;
            for j = 1:length(clusters)
                count = sum(spikeClusters(ind) == clusters(j));
                if count > mx
                    mx = count;
                    mxC = clusters(j);
                end
            end
            clusters = mxC;
        end
        clusterWaveforms_curated(cl,:,:) = clusterWaveforms(clusters+1,:,:);
        clusterWfChannels_curated(cl,:) = clusterWfChannels(clusters+1,:);
    end
    % save data
    writeNPY(spikeTimes(spike_ind), fullfile(ft, 'spikes.times.npy'))
    writeNPY(spikeClusters_curated(spike_ind), ...
        fullfile(ft, 'spikes.clusters.npy'))
    writeNPY(spikeAmplitudes(spike_ind), fullfile(ft, 'spikes.amps.npy'))
    writeNPY(spikeDepths(spike_ind), fullfile(ft, 'spikes.depths.npy'))
    writeNPY(clusters_good, fullfile(ft, 'clusters.ids.npy'))
    writeNPY(clusterWaveforms_curated, ...
        fullfile(ft, 'clusters.waveforms.npy'))
    writeNPY(clusterWfChannels_curated, ...
        fullfile(ft, 'clusters.waveformsChannels.npy'))
    copyfile(fullfile(fAlf, 'channels.localCoordinates.npy'), ...
        fullfile(ft, 'channels.localCoordinates.npy'));
    clear spikeTimes spikeClusters spikeAmplitudes spikeDepths 
    clear clusterWaveforms clusterWfChannels spikeClusters_curated
    clear clusterWaveforms_curated clusterWfChannels_curated

    % copy grating data
    if isfile(fullfile(fs, "gratings.startTime.npy"))
        descr_py = np.load(fullfile(fs, "gratingsExp.description.npy"), ...
            pyargs("allow_pickle", true));
        numExp = double(descr_py.size);
        expDescr = cell(numExp, 1);
        for ex = 1:numExp
            descr = pyrun("b = c[ind,0]", "b", c = descr_py, ind = int32(ex-1));
            expDescr{ex} = char(descr);
        end
        gratingExp = find(strcmp('Gratings', expDescr));
        % remove invalid grating experiments (e.g., where atropine was
        % used)
        if ~strcmp(db(k).gratingExperiments, "all")
            gratingExp = gratingExp(db(k).gratingExperiments);
        end

        intv = readNPY(fullfile(fs, "gratingsExp.intervals.npy"));
        if size(intv,2) == 1 % was incorrectly saved into single vector
            intv = reshape(intv, 2, [])';
        end
        intv = intv(gratingExp,:);

        t0 = readNPY(fullfile(fs, "gratings.startTime.npy"));
        t1 = readNPY(fullfile(fs, "gratings.endTime.npy"));
        dirs = double(readNPY(fullfile(fs, "gratings.direction.npy")));
        sf = readNPY(fullfile(fs, "gratings.spatialF.npy"));
        tf = readNPY(fullfile(fs, "gratings.temporalF.npy"));
        contrast = readNPY(fullfile(fs, "gratings.contrast.npy"));
        sz = [size(t0,1), size(t1,1), size(dirs,1), size(sf,1), ...
            size(tf,1), size(contrast,1)];
        if ~all(sz(1) == sz(2:end))
            fprintf("  Unequal number of rows: grating files\n")
            return
        end
        validTrials = false(size(t0));
        for j = 1:size(intv,1)
            tr = t0 >= intv(j,1) & t1 <= intv(j,2);
            nonDirFeatures = [sf(tr), tf(tr), contrast(tr)];
            uniqueStim = unique(nonDirFeatures, "rows");
            if size(uniqueStim,1) > 1
                continue
            end
            validTrials(tr) = true;
        end

        stimFeatures = [dirs sf tf contrast];
        stimFeatures = stimFeatures(validTrials,:);
        % mirror the direction angles along the vertical to make them
        % consistent with data collected at the Cortex Lab
        stimFeatures(:,1) = mod(180 - stimFeatures(:,1), 360);
        stimUnique = unique(stimFeatures, "rows");
        stimUnique = sortrows(stimUnique);
        trialStims = squeeze(all(stimFeatures == ...
            permute(stimUnique, [3 2 1]), 2));
        [trialIDs, stimIDs] = find(trialStims);
        [a, ind] = sort(trialIDs);
        stimIDs = stimIDs(ind);

        writeNPY([t0(validTrials) t1(validTrials)], ...
            fullfile(ft, "_ss_gratingsDrifting.intervals.npy"));
        writeNPY(stimIDs, ...
            fullfile(ft, "_ss_gratingsDrifting._ss_gratingsDriftingID.npy"));
        writeNPY(stimUnique(:,1), ...
            fullfile(ft, "_ss_gratingsDriftingID.directions.npy"));
        writeNPY(stimUnique(:,2), ...
            fullfile(ft, "_ss_gratingsDriftingID.spatialFrequencies.npy"));
        writeNPY(stimUnique(:,3), ...
            fullfile(ft, "_ss_gratingsDriftingID.temporalFrequencies.npy"));
        writeNPY(stimUnique(:,4), ...
            fullfile(ft, "_ss_gratingsDriftingID.contrasts.npy"));
        writeNPY(intv, fullfile(ft, "_ss_recordings.gratingsDrifting_intervals.npy"));
    end

    % copy visual noise data
    if isfile(fullfile(fs, 'sparse.startTime.npy'))
        t = readNPY(fullfile(fs, "sparse.startTime.npy"));
        map = readNPY(fullfile(fs, "sparse.map.npy"));
        edges = readNPY(fullfile(fs, "sparseExp.edges.npy"))';
        intv = readNPY(fullfile(fs, "sparseExp.intervals.npy"));
        if size(intv, 2) ~= 2
            fprintf("  Incorrect interval format: sparse noise files\n")
            return
        end
        if size(edges,2) < 4 % only number of rows and columns were saved
            fprintf("  Sparse noise: position of stimulus edges missing. -> Save NaNs.\n")
            edges = NaN(size(edges,1), 4);
        else
            if size(t,1) ~= size(map,1) || size(edges,1) ~= size(intv,1)
                fprintf("  Unequal number of rows: sparse noise files\n")
                return
            end
            if size(intv,1) > 1 && size(edges,1) == 1
                edges = reshape(edges', [], size(intv,1))';
            end
            if size(edges,2) > 6
                fprintf("  Sparse noise: stimulus edge data incorrect.\n")
                return
            end
            % original edges: [rows, columns, bottom, top, left, right]
            % or [columns, bottom, top, left, right]
            % where bottom and top are positive if above horizon;
            % transform to: [left, right, top, bottom] where bottom and
            % top are negative if above horizon
            edges = edges(:, end-3 : end);
            edges = [edges(:, [3 4]) -edges(:, [2 1])];
        end
        map(map == 0) = -1;
        map(map == 0.5) = 0;
        writeNPY(t, fullfile(ft, "_ss_sparseNoise.times.npy"));
        writeNPY((1:size(map,1))', ...
            fullfile(ft, "_ss_sparseNoise._ss_sparseNoiseID.npy"));
        writeNPY(map, fullfile(ft, "_ss_sparseNoiseID.map.npy"));
        writeNPY(edges, fullfile(ft, "_ss_sparseNoiseArea.edges.npy"));
        writeNPY(intv, fullfile(ft, "_ss_recordings.sparseNoise_intervals.npy"));
    end

    % copy circles data
    if isfile(fullfile(fs, "circles.startTime.npy"))
        t = readNPY(fullfile(fs, "circles.startTime.npy"));
        diam = readNPY(fullfile(fs, "circles.diameter.npy"));
        x = readNPY(fullfile(fs, "circles.x.npy"));
        y = readNPY(fullfile(fs, "circles.y.npy"));
        isWhite = readNPY(fullfile(fs, "circles.isWhite.npy"));
        sz = [size(t,1), size(diam,1), size(x,1), size(y,1), ...
            size(isWhite,1)];
        if ~all(sz(1) == sz(2:end))
            fprintf("  Unequal number of rows: circle files\n")
            fprintf("  %d ", sz)
            fprintf("\n  --> Drop excess entries from the end.\n")
            if sz(1) < sz(2)
                diam(sz(1)+1:end) = [];
                x(sz(1)+1:end) = [];
                y(sz(1)+1:end) = [];
                isWhite(sz(1)+1:end) = [];
            end
        else
            intv = readNPY(fullfile(fs, "circlesExp.intervals.npy"));
            if size(intv, 2) ~= 2
                fprintf("  Incorrect interval format: circles files\n")
                return
            end
        end
        writeNPY(t, fullfile(ft, "circles.times.npy"));
        writeNPY(diam, fullfile(ft, "circles.diameters.npy"));
        writeNPY(x, fullfile(ft, "circles.xPos.npy"));
        writeNPY(y, fullfile(ft, "circles.yPos.npy"));
        writeNPY(isWhite, fullfile(ft, "circles.isWhite.npy"));
        writeNPY(intv, fullfile(ft, "_ss_recordings.circles_intervals.npy"));
        clear t diam x y isWhite intv
    end
end