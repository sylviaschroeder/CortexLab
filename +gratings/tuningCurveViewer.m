function tuningCurveViewer(spikeTimes, spikeClus, cluIDs, cluYpos, ...
    stimMatrix, time, directions, blanks, laserOn, running)

stimOffsets = [0 0];

binSize = diff(time([1 2]));
offsetSamples = round(stimOffsets ./ binSize);

dirs = unique(directions(:,1));
laserConds = [0 1];
runConds = [0 1];

stimIDPerCond = NaN(length(dirs)+1,2);

for l = 1:2
    for s = 1:length(dirs)
        stimID = directions(directions(:,1)==dirs(s),2);
        stimID = stimID(laserOn(stimID) == laserConds(l));
        stimIDPerCond(s,l) = stimID;
    end
    stimID = blanks(laserOn(blanks) == laserConds(l));
    stimIDPerCond(end,l) = stimID;
end

timePlus = [time(1)-binSize time time(end)+binSize];
[yPosSorted,ind] = sort(cluYpos,'descend');
neurons = cluIDs(ind);
for iCell = 1:length(neurons)
    st = spikeTimes(spikeClus==neurons(iCell));
    if isempty(st)
        continue
    end
    spikeCount = hist(st, timePlus);
    spikeCount([1 end]) = [];
    respPerTrial = ssLocal.getTracesPerStimulus(spikeCount(:), stimMatrix, offsetSamples);
    response = squeeze(mean(respPerTrial,4)) ./ binSize; % [stimulus x repetition]
    if nansum(response(:)) == 0
        continue
    end
    means = nan(length(dirs)+1,2,2);
    sems = nan(length(dirs)+1,2,2);
    for s = 1:size(stimIDPerCond,1)
        for l = 1:2
            for r = 1:2
                stim = stimIDPerCond(s,l);
                reps = running(stim,:)==runConds(r);
                means(s,l,r) = mean(response(stim,reps));
                sems(s,l,r) = std(response(stim,reps)) / sqrt(sum(reps));
            end
        end
    end
    maxi = max(means(:)+sems(:));
    mini = min(means(:)-sems(:));
    range = maxi-mini;
    maxi = maxi+.1*range;
    mini = mini-.1*range;
    
    figure('Position', [374 662 1335 420])
    for l=1:2
        h = [0 0];
        subplot(1,2,l)
        hold on
        h(1) = errorbar(dirs, means(1:end-1,l,1), sems(1:end-1,l,1), 'ko-');
        errorbar(400, means(end,l,1), sems(end,l,1), 'ko')
        h(2) = errorbar(dirs, means(1:end-1,l,2), sems(1:end-1,l,2), 'ro-');
        errorbar(400, means(end,l,2), sems(end,l,2), 'ro')
        legend(h,'not running','running')
        set(gca,'XTick',[dirs;400],'XTickLabel',[cellstr(num2str(dirs)); 'blank'])
        xlim([-10 410])
        ylim([mini maxi])
        xlabel('Direction')
        ylabel('Mean firing rate')
    end
    subplot(1,2,1)
    title('No laser (V1 normal)')
    subplot(1,2,2)
    title('With laser (V1 inactivated)')
    annotation('textbox', [.02 .85 .09 .08], 'String', ['Cell ' ...
        num2str(neurons(iCell))], 'FontSize', 15, 'FontWeight', 'bold', ...
        'LineStyle', 'none')
    depth = cluYpos(cluIDs==neurons(iCell));
    annotation('textbox', [.02 .75 .09 .08], 'String', ['Depth: ' ...
        num2str(depth)], 'FontSize', 10, 'LineStyle', 'none')
    
    pause
%     close gcf
end