% TODO: check firing rate -> SS093 24.05.2018 K1 unit 134!

%% Parameters
limits = [-0.2 0.6];
moveTimeThresh = 0.2;
minTrials = 3;
binSize = 0.001;
sigma = 0.01;

%% Folders
folderData = '\\zubjects.cortexlab.net\Subjects';
folderPlots = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\task_ephys\PSTHs';

folderScript = 'C:\dev\workspace\CortexLab';
folderTools = 'C:\STORAGE\workspaces';

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')));
addpath(genpath(fullfile(folderTools, 'spikes')));
addpath(genpath(fullfile(folderTools, 'wheelAnalysis')));
addpath(genpath(fullfile(folderScript)));

%% Define datasets
subject = 'SS093';
date = '2018-05-24';

folderAlf = fullfile(folderData, subject, date, 'alf');
colors = lines(8);

%% Collect relevant event times etc for data set
% stim times, movement times, go cue, ...
taskData = io.getTaskInfo(folderAlf);
tlData = io.getTimelineData(folderAlf);

taskTime = [min(taskData.stimOnIntv(:)) - 20, taskData.stimOnIntv(end,2) + 60];
contrasts = [taskData.contrastLeft taskData.contrastRight];
uniCon = {unique(contrasts(:,1)), unique(contrasts(:,2))};
cIDs = zeros(size(contrasts));
[~,~,cIDs(:,1)] = histcounts(contrasts(:,1), [uniCon{1}', uniCon{1}(end)+1]);
[~,~,cIDs(:,2)] = histcounts(contrasts(:,2), [uniCon{2}', uniCon{2}(end)+1]);
groupsLeft = (cIDs(:,1)-1).*length(uniCon{2}) + cIDs(:,2);
namesLeft = num2str([reshape(repmat(uniCon{1},1,length(uniCon{2}))',[],1), ...
    repmat(uniCon{2},length(uniCon{1}),1)], '%.2f  %.2f\n');
groupsRight = (cIDs(:,2)-1).*length(uniCon{1}) + cIDs(:,1);
namesRight = num2str([repmat(uniCon{1},length(uniCon{2}),1), ...
    reshape(repmat(uniCon{2},1,length(uniCon{1}))',[],1)], '%.2f  %.2f\n');
colWeights = linspace(0, 0.7, length(uniCon{2}));
stimColors = permute(lines(length(uniCon{1})), [2 3 1]); 
stimColors = stimColors.*(1-colWeights) + zeros(3,1,length(uniCon{1})).*colWeights;
stimColors = reshape(stimColors, 3, [])';
stimColors(1,:) = 0;

%% Data for each probe
probeInfo = readtable(fullfile(folderAlf, '_ss_probes.location.csv'));
probeHems = probeInfo.hemisphere;
probeLoc = probeInfo.location;
probeNames = probeInfo.probe;
for probe = 1:size(probeInfo,1)
    folderProbe = fullfile(folderAlf, probeNames{probe});
    spikeData = io.getSpikeData(folderProbe);
    units = spikeData.IDs(spikeData.groups == 2); % only look at SUA
    if strcmp('right', probeHems(probe))
        stimGroups = groupsLeft;
        stimNames = namesLeft;
        numContra = length(uniCon{1});
        numIpsi = length(uniCon{2});
    else
        stimGroups = groupsRight;
        stimNames = namesRight;
        numContra = length(uniCon{2});
        numIpsi = length(uniCon{1});
    end
    lgdInds = [1 numIpsi.*(1:numContra)];

    folderProbePlots = fullfile(folderPlots, ...
        sprintf('%s_%s_%s', subject, date, probeNames{probe}));
    if ~isfolder(folderProbePlots)
        mkdir(folderProbePlots)
    end
    %% Plot PSTH and firing rate for each neuron
    % For each neuron, rasters for stim onset (separated by stim ID), go cue,
    % movement onset (separated by before/after go cue, left/right), feedback
    % (separated by pos/neg), stim offset (separated by stim ID)
    for n = 1:length(units)
        st = spikeData.times(spikeData.clusters==units(n));
        % align spikes and events to stimulus onset
        [st_al, st_trial] = ephys.alignData(st, ...
            taskData.stimOnIntv(:,1), limits);
        if isempty(st_al)
            continue
        end
        [moveOn_al, moveOn_trial] = ephys.alignData(tlData.wheelMoves(:,1), ...
            taskData.stimOnIntv(:,1), limits);
        [moveOff_al, moveOff_trial] = ephys.alignData(tlData.wheelMoves(:,2), ...
            taskData.stimOnIntv(:,1), limits);
        [beep_al, beep_trial] = ephys.alignData(taskData.goCueTimes, ...
            taskData.stimOnIntv(:,1), limits);
        % determine firing rate averaged across groups of trials, neglect
        % trials with movements
        % [put this into function?!]
        badTrials = any(tlData.wheelMoves(:,1)' < taskData.stimOnIntv(:,1) & ...
            tlData.wheelMoves(:,2)' > taskData.stimOnIntv(:,1), 2); % stim onset during movement
        earlyMoveTrials = tlData.wheelMoves(:,1)' - taskData.stimOnIntv(:,1);
        earlyMoveTrials(earlyMoveTrials < 0) = Inf;
        badTrials = badTrials | min(earlyMoveTrials, [], 2) < moveTimeThresh;
        goodTrials = find(~badTrials);
        valid = ismember(st_trial, goodTrials);
        stGr = stimGroups;
        stGr(badTrials) = NaN;
        [traces, bins] = ephys.tracesFromPSTH(st_al(valid), st_trial(valid), ...
            stGr, limits, binSize, sigma);
        if size(traces,1) < max(stimGroups)
            traces = padarray(traces, max(stimGroups) - size(traces,1), ...
                NaN, 'post');
        end
        
        figure('Position', [100 40 900 950])
        subplot(3,1,1:2)
        hold on
        plot([0 0], [0, max(st_trial)+1], 'k')
        plot([1 1].*moveTimeThresh, [0, max(st_trial)+1], 'k--')
        lgd = ephys.plotPSTH(gca, st_al, st_trial, limits, ...
            stimGroups, stimNames, ...
            {moveOn_al, moveOn_trial, '>', colors(1,:), 'move start'; ...
            moveOff_al, moveOff_trial, '<', colors(2,:), 'move end'; ...
            beep_al, beep_trial, 'o', colors(3,:), 'go cue'});
        lgd.Location = 'NorthEastOutside';
        xticklabels('')
        ylabel('Contrasts')
        title(sprintf('Unit %d, depth: %d um, %s %s', units(n), ...
            round(mean(spikeData.depths(spikeData.clusters==units(n)))), ...
            probeHems{probe}, probeLoc{probe}))
        
        yDist = max(traces(:)) * 0.9;
        subplot(3,1,3)
        hold on
        pos = get(gca, 'Position');
        pos(4) = 0.39 - pos(2);
        set(gca, 'Position', pos)
        h = zeros(size(traces,1),1);
        for c = 1:numContra
            for j = numIpsi:-1:1
                gr = (c-1)*numContra + j;
                h(gr) = plot(bins, traces(gr,:) - (c-1)*yDist, ...
                    'Color', stimColors(gr,:), 'LineWidth', 2);
            end
        end
        yt = yticks;
        yt = [(-numContra+1:0).*yDist max([yt(end) round(yDist/2)])];
        yticks(yt)
        yticklabels([zeros(1,numContra) yt(end)])
        yl = [yt(1) yt(end)];
        ylim(yl)
        plot([0 0], yl, 'k')
        plot([1 1].*moveTimeThresh, yl, 'k--')
        legend(h, stimNames, 'Location', 'NorthEastOutside')
        xlabel('Time from stimulus onset')
        ylabel('Firing rate (sp/s)')
        
        saveas(gcf, fullfile(folderProbePlots, ...
            sprintf('unit%04d_stimOnset.jpg', units(n))))
        close(gcf)
    end
end