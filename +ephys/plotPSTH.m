function lgd = plotPSTH(gca, spikes_al, spikes_trials, limits, groups, ...
    groupNames, groupColors, markers)

% Plot raster of grouped events (spikes) and add other (rare) events (e.g.
% behavioural events like wheel movements)
%
% lgd               Legend object (for marker names)
%
% gca               Axes where plot will appear
% spikes_al         [st x 1]; times of all spikes aligned to their relative
%                   events
% spikes_trials     [st x 1]; ID of event (trial) that each spike was
%                   aligned to
% limits            [pre post]; time spans before and after each aligned
%                   event that will appear in plot
% groups            [trial x 1]; group ID of each alignment event (trial)
% groupNames        [group x 1]; names of groups
% groupColors       [group x 3]; colours for dots of each group
% markers           {ev x 5}; other events whose times will be plotted into
%                   the raster (movement onsets); each row contains: (1)
%                   times of ev, (2) ID of ev (trial), (3) marker symbol
%                   (e.g. '+'), (4) [R G B] colour of marker, (5) name of
%                   ev (e.g. 'move start')

axes(gca)
hold on

groups_uni = unique(groups);

if isempty(groupColors)
    groupColors = zeros(length(groups_uni), 3);
end

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
        plot(st, ones(length(st),1).*k, '.', 'MarkerSize', 4, ...
            'Color', groupColors(gr,:))
    end
end
axis([limits 0 groupBorders(end)+1])
set(gca, 'YDir', 'reverse')
yticks(groupBorders - diff([0;groupBorders])./2)
yticklabels(groupNames)
if ~isempty(markers)
    lgd = legend(h_markers, markers(:,end));
end
hold off