%% Folders
% UCL
% folderData = 'C:\STORAGE\OneDrive - University of Sussex\Lab\DATA\DataToPublish\arousal_NYP_matlab\sc neurons ephys';
% folderTools = 'C:\STORAGE\workspaces';
% folderMyRepos = 'C:\dev\workspace';
% Sussex
folderData = 'C:\Users\Sylvia\OneDrive - University of Sussex\Lab\DATA\DataToPublish\arousal_NYP_matlab\sc neurons ephys';
folderTools = 'C:\dev\toolboxes';
folderMyRepos = 'C:\dev\workspaces';
folderPlots = 'C:\Users\Sylvia\OneDrive - University of Sussex\Lab\RESULTS\Electrophys_SC\V1_inactivation\PSTHs';

%% Parameters
binSize = 0.005;
prePostStim = 0.3;
smoothingSigma = 0.02;

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(genpath(fullfile(folderTools, 'CircStat2012a')))
addpath(genpath(fullfile(folderTools, 'spikes')))
addpath(genpath(fullfile(folderMyRepos, 'CortexLab')))
addpath(fullfile(folderMyRepos, 'schroeder-et-al-2020'))

%% Plot rasters and PSTHs for all neurons

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
        depths = readNPY(fullfile(folder, 'spike.depths.npy'));
        sc_depth = readNPY(fullfile(folder, 'probe.scDepth.npy'));
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
        % use same ID for blank trials
        inds = find(isnan(directions) & ~laserOn);
        stimSequence(ismember(stimSequence, inds(2:end))) = inds(1);
        stimOff = setdiff(find(~laserOn), inds(2:end));
        inds = find(isnan(directions) & laserOn);
        stimSequence(ismember(stimSequence, inds(2:end))) = inds(1);
        stimOn = setdiff(find(laserOn), inds(2:end));

        fPlots = fullfile(folderPlots, subjects{s}, dates{d});
        if ~isfolder(fPlots)
            mkdir(fPlots)
        end
        
        for n = 1:length(ids)
            spDepth = mean(depths(clusters == ids(n) & ~isnan(depths)));
            if spDepth < sc_depth(1) || spDepth > sc_depth(2)
                continue
            end
            spDepth = (sc_depth(2) - spDepth) / diff(sc_depth);
            [spAligned, trials] = ephys.alignData(st(clusters == ids(n)), ...
                stimTimes(:,1), window);

            figure('Position', [680 260 1115 715])
            sgtitle(sprintf('Neuron %03d (depth: %d%% from top of SC)', ...
                ids(n), round(spDepth * 100)))
            
            % plot raster of laser-off trials
            trOff = find(ismember(stimSequence, find(~laserOn)));
            valid = ismember(trials, trOff);
            sp = spAligned(valid);
            tr = trials(valid);
            seq = NaN(size(stimSequence));
            seq(trOff) = stimSequence(trOff);
                        
            subplot(2,2,3)
            ephys.plotPSTH(gca, sp, tr, window, seq, ...
                [dirs; NaN], [cols; 0 0 0], {});
            hold on
            % plot stimulus on- and offset
            plot([0 0], ylim, 'k');
            plot([1 1] .* stimDur, ylim, 'k');
            xlim(window)
            xlabel('Time from visual stimulus (s)')
            ylabel('Direction')

            % plot traces for laser-off trials
            subplot(2,2,1)
            [traces, bins] = ephys.tracesFromPSTH(sp, tr, ...
                seq, window, binSize, smoothingSigma, false);
            plot(bins, traces(stimOff,:))
            colororder(gca, [cols; 0 0 0])
            xlim(window)
            ylabel('Firing rate (sp/s)')
            title('Laser Off')

            % plot raster of laser-on trials
            trOn = find(ismember(stimSequence, find(laserOn)));
            valid = ismember(trials, trOn);
            sp = spAligned(valid);
            tr = trials(valid);
            seq = NaN(size(stimSequence));
            seq(trOn) = stimSequence(trOn);
                        
            subplot(2,2,4)
            ephys.plotPSTH(gca, sp, tr, window, seq, ...
                [dirs; NaN], [cols; 0 0 0], {});
            hold on
            % plot stimulus on- and offset
            plot([0 0], ylim, 'k');
            plot([1 1] .* stimDur, ylim, 'k');
            % plot laser on- and offset
            plot([1 1] .* lOn, ylim, 'c');
            plot([1 1] .* (stimDur + lOff), ylim, 'c');
            xlim(window)
            xlabel('Time from visual stimulus (s)')
            ylabel('Direction')

            % plot traces for laser-on trials
            subplot(2,2,2)
            [traces, bins] = ephys.tracesFromPSTH(sp, tr, ...
                seq, window, binSize, smoothingSigma, false);
            plot(bins, traces(stimOn,:))
            colororder(gca, [cols; 0 0 0])
            xlim(window)
            ylabel('Firing rate (sp/s)')
            title('Laser On')

            % set y-axes to same scale for trace plots
            subplot(2,2,1)
            y1 = ylim;
            subplot(2,2,2)
            y2 = ylim;
            y = max(y1(2), y2(2));
            for sp = 1:2
                subplot(2,2,sp)
                ylim([0 y])
            end

            saveas(gcf, fullfile(fPlots, sprintf('neuron%03d.png', ids(n))))
            close gcf
        end
    end
end

%% Plot PSTHs for example neurons and stimuli
subject = 'SS061';
date = '2016-05-11';
units =   [80 91 190 194 379];
stimuli = [ 9 10   9   5  10];
cols = [0 1 1; 0 0 0; 0.75 0.25 0.25; 1 0 0];
% cols = [0 0.5 0.5; 0 0 0];

folder = fullfile(folderData, subject, date, '001');
st = readNPY(fullfile(folder, 'spike.times.npy'));
clusters = readNPY(fullfile(folder, 'spike.clusters.npy'));
depths = readNPY(fullfile(folder, 'spike.depths.npy'));
sc_depth = readNPY(fullfile(folder, 'probe.scDepth.npy'));
stimTimes = readNPY(fullfile(folder, '_ss_grating.intervals.npy'));
stimSequence = readNPY(fullfile(folder, '_ss_grating._ss_gratingID.npy'));
directions = readNPY(fullfile(folder, '_ss_gratingID.directions.npy'));

stimDur = mean(diff(stimTimes,1,2));
window = [-prePostStim stimDur+prePostStim];

fPlots = 'C:\Users\Sylvia\OneDrive - University of Sussex\Funding\2022-01_BBSRC Standard\Figures\Plots\V1 inactivation';
if ~isfolder(fPlots)
    mkdir(fPlots)
end

for n = 1:length(units)
    spDepth = mean(depths(clusters == units(n) & ~isnan(depths)));
    spDepth = (sc_depth(2) - spDepth) / diff(sc_depth);
    [spAligned, trials] = ephys.alignData(st(clusters == units(n)), ...
        stimTimes(:,1), window);

    stim = [13 26 stimuli(n) stimuli(n)+13]; % blank /w laser, blank /wo laser, 
    % grating /w laser, grating /wo laser, 
%     stim = [stimuli(n) stimuli(n)+13]; % grating /w laser, grating /wo laser

    % get traces for specific stimuli
    trStim = find(ismember(stimSequence, stim));
    valid = ismember(trials, trStim);
    sp = spAligned(valid);
    tr = trials(valid);
    seq = NaN(size(stimSequence));
    seq(trStim) = stimSequence(trStim);
    for s = 1:length(stim)
        seq(seq == stim(s)) = s;
    end
    [traces, bins] = ephys.tracesFromPSTH(sp, tr, ...
        seq, window, binSize, smoothingSigma, false);
%     % subtract prestimulus activity
%     ind = bins < 0;
%     preStim = mean(traces(:,ind), 2);
%     traces = traces - preStim;

    figure
    sgtitle(sprintf('Neuron %03d (depth: %d%% from top of SC)', ...
        units(n), round(spDepth * 100)))

    subplot(2,1,1)
    % mark time of stim presentation (simultaneous with laser)
    maxi = max(ceil(traces(:)/10) * 10);
    fill([0 2 2 0], [0 0 maxi maxi], 'k', 'EdgeColor', 'none', ...
        'FaceColor', 'k', 'FaceAlpha', 0.2)
    hold on

    % plot traces
    h = zeros(1, size(cols,1));
    for s = 1:length(h)
        h(s) = plot(bins, traces(s,:), 'Color', cols(s,:), 'LineWidth', 1);
    end
    xlim(window)
    ylim([0 maxi])
    set(gca, 'box', 'off')
    ylabel('Firing rate (sp/s)')
    legend(h, {'no stim + laser', 'no stim', 'stim + laser', 'stim'})
%     legend(h, {'stim + laser', 'stim only'})

    subplot(2,1,2)
    % mark time of stim presentation (simultaneous with laser)
    maxi = length(trStim)+2;
    fill([0 2 2 0], [0 0 maxi maxi], 'k', 'EdgeColor', 'none', ...
        'FaceColor', 'k', 'FaceAlpha', 0.2)
    hold on

    % plot raster
    ephys.plotPSTH(gca, sp, tr, window, seq, ...
        stim, cols, {});
    set(gca, 'YDir', 'normal', 'box', 'off')
    xlim(window)
    ylim([0 maxi-1])
    xlabel('Time from stimulus (s)')
    ylabel('Direction')

    savefig(gcf, fullfile(fPlots, sprintf('%s_%s_neuron%03d.fig', subject, date, units(n))), 'compact')
    saveas(gcf, fullfile(fPlots, sprintf('%s_%s_neuron%03d.png', subject, date, units(n))))
    close gcf
end