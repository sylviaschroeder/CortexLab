dirs = 0:30:330;
xlimits = [-10 370];
stimEnd = round(results.plane(1).stimDuration / ...
    median(diff(results.plane(1).kernelTime)));
folder = 'C:\RESULTS\boutons\nonVisualEffects\tuningCurves\M160923_SS069_2016-10-13';
for ip = 1:length(results.plane)
    for ic = 1:length(results.plane(ip).kernelFit)
        if isempty(results.plane(ip).kernelFit(ic).kernel)
            continue
        end
        m = mean(results.plane(ip).kernelFit(ic).alphaEachTrial,1);
        s = std(results.plane(ip).kernelFit(ic).alphaEachTrial,0,1) ./ ...
            sqrt(size(results.plane(ip).kernelFit(ic).alphaEachTrial,1));
        
        kernelSign = sign(sum(results.plane(ip).kernelFit(ic).kernel(1:stimEnd)));
        m = m * kernelSign;
        s = s * kernelSign;
        
        figure
        hold on
%         fill(xlimits([1 2 2 1]), [[1 1].*(m(end)+s(end)), ...
%             [1 1].*(m(end)-s(end))], 'k', 'EdgeColor', 'none', ...
%             'FaceColor', 'k', 'FaceAlpha', .3)
%         plot(xlimits, [1 1].*m(end), 'k:', 'LineWidth', 2)
        errorbar([dirs 360], m([1:end-1 1]), s([1:end-1 1]), 'k', 'LineWidth', 2)
        plot(xlimits, [0 0], 'k:')
        xlim(xlimits)
        xlabel('Direction')
        ylabel('\DeltaF/F')
        title(sprintf('Plane %d, bouton %d', results.planes(ip), ic))
        
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(fullfile(folder, sprintf('plane%02d_ROI%03d.jpg', ...
            results.planes(ip), ic)), '-djpeg','-r0')
        close gcf
        
%         figure
%         plot(results.plane(ip).kernelFit(ic).kernel)
    end
end