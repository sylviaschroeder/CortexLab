function figHandles = plotStimRespVsNonVis(params, nonVisual)

% params        output of nonVis.modelResponsePerTrial
%               .alphaEachTrial     [trial x stimulus]
%               .alphas             [stimulus x 1]
%               .beta               double
%               .gamma              double
% nonVisual     [trial x stimulus]

numStim = size(params.alphaEachTrial,2);

colors = jet(numStim);
mini = min(nonVisual(:));
maxi = max(nonVisual(:));
range = maxi-mini;
mini = mini-0.05*range;
maxi = maxi+0.05*range;

% full model
figHandles = [0 0];
figHandles(1) = figure;
hold on
for iStim = 1:numStim
    [nvSorted, ind] = sort(nonVisual(:,iStim));
    plot(nvSorted, params.alphaEachTrial(ind,iStim), 'o-', 'Color', colors(iStim,:))
    if isfield(params, 'alphas')
        plot([mini maxi], params.alphas(iStim) + ...
            (params.alphas(iStim)*params.beta + params.gamma) .* [mini maxi], ...
            'Color', colors(iStim,:), 'LineWidth', 2)
    end
end
% plot([mini maxi], [1 1]*abs(params.gamma)/abs(params.beta), 'k:', 'LineWidth', 2)
xlim([mini maxi])
xlabel('Nonvisual signal (per trial)')
ylabel('Stimulus response (per trial)')

colormap(colors)
c = colorbar;
c.Label.String = 'Stimulus';
c.Ticks = (.5:numStim) / numStim;
c.TickLabels = 1:numStim;
c.Direction = 'reverse';

% each stimulus is fit with separate line
if isfield(params, 'alphasUncond')
    figHandles(2) = figure;
    hold on
    for iStim = 1:numStim
        [nvSorted, ind] = sort(nonVisual(:,iStim));
        plot(nvSorted, params.alphaEachTrial(ind,iStim), 'o-', 'Color', colors(iStim,:))
        plot([mini maxi], params.alphasUncond(iStim) + ...
            params.gammasUncond(iStim) .* [mini maxi], ...
            'Color', colors(iStim,:), 'LineWidth', 2)
    end
    xlim([mini maxi])
    xlabel('Nonvisual signal (per trial)')
    ylabel('Stimulus response (per trial)')
    
    colormap(colors)
    c = colorbar;
    c.Label.String = 'Stimulus';
    c.Ticks = (.5:numStim) / numStim;
    c.TickLabels = 1:numStim;
    c.Direction = 'reverse';
end