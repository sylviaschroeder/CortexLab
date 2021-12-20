%% Folders
folderData = 'C:\STORAGE\OneDrive - University of Sussex\Lab\DATA\DataToPublish\arousal_NYP_matlab\sc neurons ephys';
folderTools = 'C:\STORAGE\workspaces';
folderMyRepos = 'C:\dev\workspace';

%% Parameters
binSize = 0.005;
prePostStim = 0.2;

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(genpath(fullfile(folderTools, 'CircStat2012a')))
addpath(genpath(fullfile(folderTools, 'spikes')))
addpath(genpath(fullfile(folderMyRepos, 'CortexLab')))
addpath(fullfile(folderMyRepos, 'schroeder-et-al-2020'))

%% Get rasters

dirs = dir(folderData);
dirs = dirs(3:end);
subjects = {dirs([dirs.isdir]).name};
for s = 1:length(subjects)
    dirs = dir(fullfile(folderData, subjects{s}));
    dirs = dirs(3:end);
    dates = {dirs.name};
    for d = 1:length(dates)
        folder = fullfile(folderData, subjects{s}, dates{d}, '001');
        ids = readmatrix(fullfile(folder, 'clusters.uuids.csv'));
        st = readNPY(fullfile(folder, 'spike.times.npy'));
        clusters = readNPY(fullfile(folder, 'spike.clusters.npy'));
        stimTimes = readNPY(fullfile(folder, '_ss_grating.intervals.npy'));
        stimSequence = readNPY(fullfile(folder, '_ss_grating._ss_gratingID.npy'));
        directions = readNPY(fullfile(folder, '_ss_gratingID.directions.npy'));
        laserOn = readNPY(fullfile(folder, '_ss_gratingID.laserOn.npy'));
        laserOnTimes = readNPY(fullfile(folder, '_ss_gratingID.laserOnTime.npy'));
        laserOffTimes = readNPY(fullfile(folder, '_ss_gratingID.laserOffTime.npy'));
        
        lOn = mean(laserOnTimes);
        lOff = mean(laserOffTimes);
        stimDur = mean(diff(stimTimes,1,2));
        window = [-prePostStim stimDur+prePostStim];
        dirs = unique(directions(~isnan(directions)));
        cols = hsv(length(dirs));
        cols_laser = 0.8 .* cols;
        % use same ID for blank trials
        inds = find(isnan(directions) & ~laserOn);
        stimSequence(ismember(stimSequence, inds(2:end))) = inds(1);
        inds = find(isnan(directions) & laserOn);
        stimSequence(ismember(stimSequence, inds(2:end))) = inds(1);
        
        for n = 1:length(ids)
            [spAligned, trials] = ephys.alignData(st(clusters == ids(n)), ...
                stimTimes(:,1), window);
            figure
            
            % plot raster of laser-off trials
            trOff = find(ismember(stimSequence, find(~laserOn)));
            valid = ismember(trials, trOff);
            sp = spAligned(valid);
            tr = trials(valid);
                        
            subplot(2,2,3)
            ephys.plotPSTH(gca, sp, tr, window, stimSequence(trOff), ...
                [dirs; NaN], [cols; 0 0 0], {});
            
            for sp = 1:4
                subplot(2,2,sp), hold on
            end
            
            
            for stim = 1:length(dirs)
                stimID = find(directions == dirs(stim) & ~laserOn);
                stimTrials = find(stimSequence == 15);
            end
        end
    end
end