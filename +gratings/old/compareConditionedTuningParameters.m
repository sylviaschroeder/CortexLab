function compareConditionedTuningParameters(fitParameters, baselines, labels, ...
    adjustedRsquared)
% fitParameters     [neuron x condition x parameter]
% labels            {1 x 2}; strings with condition names
% adjustedRsquared  [neuron x condition]; contains adjusted R_squared of
%                   fit for each neuron and each condition (plus all data 
%                   not divided into conditions)

Rthreshold = 0.3;
parameters = {'Pref. direction', 'Amplitude at preferred direction', ...
    'Amplitude at null direction', 'Offset', 'Tuning width'};

for par = 2:4
    figure('Position', [860 680 1040 420])
    
    subplot(1,2,1)
    [f, gof] = fit(fitParameters(:,1,par), fitParameters(:,2,par), 'a*x');
    plot(fitParameters(:,1,par), fitParameters(:,2,par), 'k.')
    hold on
    mini = min(reshape(fitParameters(:,[1 2],par),[],1));
    maxi = max(reshape(fitParameters(:,[1 2],par),[],1));
    plot([mini maxi], f([mini maxi]), 'r')
    axis([mini maxi mini maxi])
    a = coeffvalues(f);
    confidenceInt = confint(f);
    text(mini+0.05*(maxi-mini), mini+0.9*(maxi-mini), ...
        sprintf('y = a * x\na = %.2f (%.2f %.2f)\nadj R^2 = %.3f', ...
        a, confidenceInt(1), confidenceInt(2), gof.adjrsquare))
    title(sprintf('%s (all data used)', parameters{par}))
    xlabel(labels{1})
    ylabel(labels{2})
    
    subplot(1,2,2)
    ind = all(adjustedRsquared(:,1:2) >= Rthreshold, 2);
    [f, gof] = fit(fitParameters(ind,1,par), fitParameters(ind,2,par), 'a*x');
    plot(fitParameters(ind,1,par), fitParameters(ind,2,par), 'k.')
    hold on
    mini = min(reshape(fitParameters(:,[1 2],par),[],1));
    maxi = max(reshape(fitParameters(:,[1 2],par),[],1));
    plot([mini maxi], f([mini maxi]), 'r')
    axis([mini maxi mini maxi])
    a = coeffvalues(f);
    confidenceInt = confint(f);
    text(mini+0.05*(maxi-mini), mini+0.9*(maxi-mini), ...
        sprintf('y = a * x\na = %.2f (%.2f %.2f)\nadj R^2 = %.3f', ...
        a, confidenceInt(1), confidenceInt(2), gof.adjrsquare))
    title(sprintf('%s (only good fits used)', parameters{par}))
    xlabel(labels{1})
    ylabel(labels{2})
end

% compare various combinations of parameters
newParameters = cat(3, sum(fitParameters(:,:,[2 4]), 3), ... % Rp + R0
    min(fitParameters(:,:,2) ./ fitParameters(:,:,3), ...
    ones(size(fitParameters,1), size(fitParameters,2)) * 10), ... % Rp / Rn
    (fitParameters(:,:,2)-fitParameters(:,:,3)) ./ ...
    sum(fitParameters(:,:,[2 3]),3), ...                     % (Rp-Rn) / (Rp+Rn)
    fitParameters(:,:,2) ./ bsxfun(@plus, fitParameters(:,:,4), baselines'));           % Rp / R0
parNames = {'Max. responses (Rp+R0)', ...
    'Direction selectivity (Rp/Rn)', ...
    'Direction selectivity (Rp-Rn)/(Rp+Rn)', ...
    'Tuning strength (Rp/(R0+baseline))'};
for par = 1:size(newParameters,3)
    figure('Position', [860 680 1040 420])
    
    subplot(1,2,1)
    [f, gof] = fit(newParameters(:,1,par), newParameters(:,2,par), 'a*x');
    plot(newParameters(:,1,par), newParameters(:,2,par), 'k.')
    hold on
    mini = min(reshape(newParameters(:,[1 2],par),[],1));
    maxi = max(reshape(newParameters(:,[1 2],par),[],1));
    plot([mini maxi], f([mini maxi]), 'r')
    axis([mini maxi mini maxi])
    a = coeffvalues(f);
    confidenceInt = confint(f);
    text(mini+0.05*(maxi-mini), mini+0.9*(maxi-mini), ...
        sprintf('y = a * x\na = %.2f (%.2f %.2f)\nadj R^2 = %.3f', ...
        a, confidenceInt(1), confidenceInt(2), gof.adjrsquare))
    title(sprintf('%s (all data used)', parNames{par}))
    xlabel(labels{1})
    ylabel(labels{2})
    
    subplot(1,2,2)
    ind = all(adjustedRsquared(:,1:2) >= Rthreshold, 2);
    [f, gof] = fit(newParameters(ind,1,par), newParameters(ind,2,par), 'a*x');
    plot(newParameters(ind,1,par), newParameters(ind,2,par), 'k.')
    hold on
    mini = min(reshape(newParameters(:,[1 2],par),[],1));
    maxi = max(reshape(newParameters(:,[1 2],par),[],1));
    plot([mini maxi], f([mini maxi]), 'r')
    axis([mini maxi mini maxi])
    a = coeffvalues(f);
    confidenceInt = confint(f);
    text(mini+0.05*(maxi-mini), mini+0.9*(maxi-mini), ...
        sprintf('y = a * x\na = %.2f (%.2f %.2f)\nadj R^2 = %.3f', ...
        a, confidenceInt(1), confidenceInt(2), gof.adjrsquare))
    title(sprintf('%s (only good fits used)', parNames{par}))
    xlabel(labels{1})
    ylabel(labels{2})
end