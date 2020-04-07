db_ephys_opticTract

samplingRate = 30000;

%% Parameters
binSizeFlicker = 0.001;

% colormaps
red = [1 0 .5];
% green = [0 1 .5];
blue = [0 .5 1];
% orange = [1 .5 0];
black = [1 1 1].*0.5;
grad = linspace(0,1,100)';
reds = red.*flip(grad) + [1 1 1].*grad;
% greens = green.*flip(grad) + [1 1 1].*grad;
blacks = black.*flip(grad) + [1 1 1].*grad;
% cm_ON = [greens; flip(reds(1:end-1,:),1)];
cm_ON = [blacks; flip(reds(1:end-1,:),1)];
blues = blue.*flip(grad) + [1 1 1].*grad;
% oranges = orange.*flip(grad) + [1 1 1].*grad;
% cm_OFF = [oranges; flip(blues(1:end-1,:),1)];
cm_OFF = [blacks; flip(blues(1:end-1,:),1)];
colormaps = {cm_ON, cm_OFF};

%% Folders
% file location depending on PCs
folderPC = 'C:\STORAGE\OneDrive - University College London';
folderResults = fullfile(folderPC, 'Lab\RESULTS\receptiveFields\axons');
plotFolder = fullfile(folderPC, 'Lab\RESULTS\OpticTract');

%% Plot RFs
titles = {'ON field','OFF field'};
data = load(fullfile(folderResults, 'receptiveFields.mat'));
RFs = data.RFs;
for k = 1:length(RFs)
    pos = RFs(k).stimPosition;
    good = find(RFs(k).OTgood{1});
    for iCell = 1:length(RFs(k).explainedVariances)
        rf = RFs(k).receptiveFields(:,:,:,:,iCell);
        rf(:,:,:,2) = -rf(:,:,:,2);
        squW = diff(pos(1:2)) / size(rf,2);
        squH = diff(pos(3:4)) / size(rf,1);
        [mx,mxTime] = max(max(reshape(permute(abs(rf),[1 2 4 3]),[],size(rf,3)),[],1));
        % fit Gaussian
        fitPars = whiteNoise.fit2dGaussRF(mean(rf(:,:,mxTime,:),4), false);
        [x0, y0] = meshgrid(linspace(0.5, size(rf,2)+0.5, 100), ...
            linspace(0.5, size(rf,1)+0.5, 100));
        fitRF = whiteNoise.D2GaussFunctionRot(fitPars, cat(3, x0, y0));
        [x1,y1] = meshgrid(linspace(pos(1), pos(2), 100), ...
            linspace(pos(3), pos(4), 100));
        figure('Position',[1930 690 765 420])
        for f = 1:2
            for t = 1:length(RFs(k).RFTimes)
                subplot(2,length(RFs(k).RFTimes),(f-1)*2+t)
                imagesc([pos(1)+squW/2 pos(2)-squW/2], [pos(3)+squH/2 pos(4)-squH/2], ...
                    rf(:,:,end+1-t,f),[-mx mx])
                hold on
                contour(x1, y1, fitRF, [1 1].*fitPars(1)/2, 'k', 'LineWidth', 1)
                axis image
                set(gca, 'box', 'off', 'XTick', pos(1:2), 'YTick', [pos(3) 0 pos(4)], ...
                    'YTickLabel', [-pos(3) 0 -pos(4)])
                colormap(gca, colormaps{f})
                if t == 1
                    ylabel(titles{f})
                else
                    p = get(gca, 'Position');
                    colorbar('Position',[p(1)+p(3)+.01,p(2),0.03,p(4)])
                end
                if f == 1
                    title(sprintf('%.2f s', RFs(k).RFTimes(end+1-t)))
                end
            end
        end
        annotation('textbox', [0 .95 1 .03], 'String', sprintf('%s, %s, unit %d', ...
            RFs(k).subject, RFs(k).date, RFs(k).OTunits{1}(good(iCell))), ...
            'FontSize', 10, 'FontWeight', 'bold', 'LineStyle', 'none', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
    end
end

%% Plot correlograms with running
data = load(fullfile(plotFolder, 'runningCorrelation_dark', 'correlations_runningFiltered.mat'));
corrs = data.corrs;
maxLag = 30;

for k = 1:length(corrs)
    lags = corrs(k).lags;
    for n = 1:length(corrs(k).units)
        if corrs(k).units(n).runningTime < 0.05
            continue
        end
        null = prctile(corrs(k).units(n).nullCrossCorr, [2.5 97.5], 1);
        real = corrs(k).units(n).crosscorr;
        mini = floor(min([null(:);real(:)]) * 10) / 10;
        maxi = ceil(max([null(:);real(:)]) * 10) / 10;
        figure
        hold on
        fill([lags flip(lags)], [null(1,:) flip(null(2,:))], 'k', ...
            'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.2)
        plot(lags, real, 'k', 'LineWidth', 1)
        plot(lags([1 end]), [0 0], 'k:')
        plot([0 0], [mini maxi], 'k:')
        xlim([-maxLag maxLag])
        ylim([mini maxi])
        set(gca, 'YTick', mini:0.1:maxi, 'XTick', [-maxLag 0 maxLag])
        xlabel('Time lag (s)')
        ylabel('Cross-correlation')
        title(sprintf('%s, %s, unit %d', ...
            corrs(k).subject, corrs(k).date, corrs(k).units(n).ID))
    end
end

%% Scatter: corr. with running versus p-value
data = load(fullfile(plotFolder, 'runningCorrelation_dark', 'correlations_runningFiltered.mat'));
corrs = data.corrs;

pVals = [];
rhos =[];
for k = 1:length(corrs)
    [~,zeroInd] = min(abs(corrs(k).lags));
    for n = 1:length(corrs(k).units)
        if corrs(k).units(n).runningTime < 0.05
            continue
        end
        rho = corrs(k).units(n).crosscorr(zeroInd);
        rhos(end+1,1) = rho;
        null = corrs(k).units(n).nullCrossCorr(:,zeroInd);
        p = sum(null > rho) / length(null);
        if p > 0.5
            p = 1-p;
        end
        pVals(end+1,1) = 2 * p;
    end
end
pVals(pVals==0) = 1e-4;

figure
scatter(rhos, pVals, 'k', 'filled')
hold on
plot([-.4 0.4],[1 1].*0.05, 'k')
set(gca, 'YScale', 'log', 'YDir', 'reverse', 'YMinorTick', false, 'XTick', [-.4 0 .4])
xlabel('Correlation with running')
ylabel('p-value')
title(['n = ' num2str(length(rhos))])