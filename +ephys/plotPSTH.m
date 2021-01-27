function lgd = plotPSTH(gca, spikes_al, spikes_trials, limits, groups, groupNames, ...
    markers)

axes(gca)
hold on

groups_uni = unique(groups);

k = 0;
groupBorders = zeros(length(groups_uni),1);
h_markers = NaN(size(markers,1),1);
% plot each group of trials
for gr = 1:length(groups_uni)
    grTrials = find(groups == groups_uni(gr));
    groupBorders(gr) = k + length(grTrials);
    if mod(gr,2) == 0
        fill([limits flip(limits)], k + [0 0 length(grTrials) length(grTrials)], ...
            'k', 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.1)
    end
    % plot each trial/line
    for j = 1:length(grTrials)
        k = k + 1;
        % plot markers
        for m = 1:size(markers,1)
            ind = markers{m,2} == grTrials(j);
            mrk = markers{m,1}(ind);
            if isempty(mrk)
                continue
            end
            h = plot(mrk, ones(length(mrk),1).*k, markers{m,3}, ...
                'MarkerFaceColor', markers{m,4}, 'MarkerEdgeColor', 'none');
            if isnan(h_markers(m))
                h_markers(m) = h;
            end
        end
        % plot spikes
        ind = spikes_trials == grTrials(j);
        st = spikes_al(ind);
        plot(st, ones(length(st),1).*k, 'k.', 'MarkerSize', 4)
    end
end
axis([limits 0 groupBorders(end)+1])
set(gca, 'YDir', 'reverse')
yticks(groupBorders - diff([0;groupBorders])./2)
yticklabels(groupNames)
lgd = legend(h_markers, markers(:,end));
hold off