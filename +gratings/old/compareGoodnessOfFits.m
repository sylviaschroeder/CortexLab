function compareGoodnessOfFits(adjustedRsquared, labels)

% adjustedRsquared  [neuron x fits]
% labels            {1 x fits}

figure('Position', [10 680 1900 420])
plots = sum(1:size(adjustedRsquared,2)-1);
k = 1;
maxi = max(adjustedRsquared(:));
mid = ceil(plots/2);
for type1 = 1:size(adjustedRsquared,2)-1
    for type2 = type1+1:size(adjustedRsquared,2)
        subplot(1,plots,k)
        plot(adjustedRsquared(:,type1), adjustedRsquared(:,type2), 'k.')
        hold on
        plot([0 maxi],[0 maxi], 'k:')
        axis([0 maxi 0 maxi])
        xlabel(labels{type1})
        ylabel(labels{type2})
        
        if k == mid
            title('Adjusted R^2 values of tuning curve fits')
        end
        
        k = k + 1;
    end
end