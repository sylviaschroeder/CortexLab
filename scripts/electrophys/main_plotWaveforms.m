%% Folders
folderTools = 'C:\dev\toolboxes';
folderThisRepo = 'C:\dev\workspaces\CortexLab';
folderData = 'C:\Users\Sylvia\OneDrive - University of Sussex\Projects\Ephys_task_left and right SC\Data';

%% Examples
examples = {'SS088', '2018-02-01', 'K1', 275;...
%     'SS093', '2018-05-25', 'K1', 289; ...
%     'SS093', '2018-05-25', 'K2', 193; ...
    'SS093', '2018-05-24', 'K1', 392};

%% Add paths
addpath(genpath(fullfile(folderTools, 'spikes')))
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(fullfile(folderThisRepo))

%% Waveforms
for d = 1:size(examples,1)
    folder = fullfile(folderData, examples{d,1}, examples{d,2}, examples{d,3});
    ids = readNPY(fullfile(folder, 'clusters.ids.npy'));
    clu = readNPY(fullfile(folder, 'spikes.clusters.npy'));
    amps = readNPY(fullfile(folder, 'spikes.amps.npy'));
    depths = readNPY(fullfile(folder, 'spikes.depths.npy'));
    surface = readNPY(fullfile(folder, '_ss_colliculusTop.depth.npy'));
    templates = readNPY(fullfile(folder, 'spikes.templates.npy'));
    wfs = readNPY(fullfile(folder, 'templates.waveforms.npy')); % [templates x time x chans]
    wfChan = readNPY(fullfile(folder, 'templates.waveformsChannels.npy')); % [templates x chans]
    coord = readNPY(fullfile(folder, 'channels.localCoordinates.npy'));

    allTempls = unique(templates);
    for ex = 1:length(examples{d,4})
        unit = examples{d,4}(ex);
        tem = unique(templates(clu == unit));
        amp = mean(amps(clu == unit));
        depth = mean(depths(clu == unit));
        if length(tem) > 1
            % find template occurring most often
            num = NaN(1, length(tem));
            for t = 1:length(tem)
                num(t) = sum(templates == tem(t));
            end
            [~,t] = max(num);
            tem = tem(t);
        end

        % find channel distance
        ind = find(allTempls == tem);
        ycoords = unique(coord(wfChan(ind,:),2));
        mini = min(ycoords);
        maxi = max(ycoords);
        chanDist = median(diff(ycoords));

        figure('Position', [1120 50 746 946]);
        hold on
        arrayfun(@(x)plot(coord(wfChan(ind,x),1)+0.3*(1:size(wfs,2))', ...
            coord(wfChan(ind,x),2)+10.*smooth(squeeze(wfs(ind,:,x))), 'k'), 1:size(wfs,3));
        ylim([mini - 0.5*chanDist, maxi + 0.5 * chanDist])
        ylabel('Depth (um)')
        set(gca, 'XTick', [])
        title(sprintf('%s, %s, probe %s, neuron %d (ampl.: %.1f uV, depth: %d um)', ...
            examples{d,1}, examples{d,2}, examples{d,3}, unit, amp, round(depth)))
    end
end