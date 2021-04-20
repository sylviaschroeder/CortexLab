%% Folders
folderBase = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\DataToPublish\arousal_NYP_matlab';
folderTools = 'C:\STORAGE\workspaces';
folderThisRepo = 'C:\dev\workspace\schroeder-et-al-2020';
folderPlots = 'C:\STORAGE\OneDrive - University College London\Lab\Results\boutons\nonVisualEffects\plots_stimResponses';

%% Parameters
buffer = 2; % sec (before and after stim)

%% Examples
examples = {'SS076', '2017-10-04'};
exBoutons = [24 58 175];

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(fullfile(folderThisRepo))

%% Plots
for ex = 1:size(examples,1)
    folder = fullfile(folderPlots, sprintf('%s_%s', examples{1}, examples{2}));
    if ~isfolder(folder)
        mkdir(folder)
    end
    planes = readNPY(fullfile(folderBase, 'boutons', examples{ex,1}, ...
        examples{ex,2}, '001\_ss_2pRois._ss_2pPlanes.npy'));
    ids = readNPY(fullfile(folderBase, 'boutons', examples{ex,1}, ...
        examples{ex,2}, '001\_ss_2pRois.ids.npy'));
    planeDelays = readNPY(fullfile(folderBase, 'boutons', examples{ex,1}, ...
        examples{ex,2}, '001\_ss_2pPlanes.delay.npy'));
    traces = readNPY(fullfile(folderBase, 'boutons', examples{ex,1}, ...
        examples{ex,2}, '001\_ss_2pCalcium.dff.npy'));
    time = readNPY(fullfile(folderBase, 'boutons', examples{ex,1}, ...
        examples{ex,2}, '001\_ss_2pCalcium.timestamps.npy'));
    
    stimIntervals = readNPY(fullfile(folderBase, 'boutons', examples{ex,1}, ...
        examples{ex,2}, '001\_ss_grating.intervals.npy'));
    stimSequence = readNPY(fullfile(folderBase, 'boutons', examples{ex,1}, ...
        examples{ex,2}, '001\_ss_grating._ss_gratingID.npy'));
    directions = readNPY(fullfile(folderBase, 'boutons', examples{ex,1}, ...
        examples{ex,2}, '001\_ss_gratingID.directions.npy'));
    largePupil = readNPY(fullfile(folderBase, 'boutons', examples{ex,1}, ...
        examples{ex,2}, '001\_ss_gratingTrials.largePupil.npy'));
    
    for iBout = 1:length(ids)
        t = time + planeDelays(planes(iBout));
        stimMatrix = exp.buildStimMatrix(stimSequence, stimIntervals, t);
        directions(isnan(directions)) = [];
        timeBin = median(diff(time));
        repetitions = sum(stimSequence == 1);
        stimDurInFrames = round(sum(stimMatrix(1,:)) / repetitions);
        stimDur = stimDurInFrames * timeBin;
        offset = ceil(buffer / timeBin);
        resp = squeeze(exp.getTracesPerStimulus(traces(:,iBout), stimMatrix, ...
            [1 1] .* offset)); % [stimulus x trial x time]
        
        ind = ~largePupil;
        ind = repmat(ind, 1, 1, size(resp,3));
        temp = resp(1:end-1,:,:);
        temp(~ind) = NaN;
        meanSmall = squeeze(nanmean(temp, 2));
        semSmall = squeeze(nanstd(temp, 0, 2) ./ sqrt(sum(~ind,2)));
        temp = resp(1:end-1,:,:);
        temp(ind) = NaN;
        meanLarge = squeeze(nanmean(temp, 2));
        semLarge = squeeze(nanstd(temp, 0, 2) ./ sqrt(sum(ind,2)));
        
        mini = min([meanSmall(:) - semSmall(:); meanLarge(:) - semLarge(:)]);
        maxi = max([meanSmall(:) + semSmall(:); meanLarge(:) + semLarge(:)]);
        rng = maxi - mini;
        mini = mini - 0.05*rng;
        maxi = maxi + 0.05*rng;
        xDist = .5;
        traceDur = stimDur + 2*buffer;
        respTime = (-offset:stimDurInFrames+offset-1) .* timeBin;
        
        if iBout == 1
            maxResps = NaN(length(ids), length(respTime), 2);
            baseInds = respTime < 0;
            respInds = respTime >= 0 & respTime <= stimDur;
        end
        amps = nanmean(meanSmall(:,respInds),2) - ...
            nanmean(meanSmall(:,baseInds),2);
        if any(amps > 0)
            [~,maxStim] = max(amps);
            maxResps(iBout,:,1) = meanSmall(maxStim,:);
            maxResps(iBout,:,2) = meanLarge(maxStim,:);
        end
        
        figure('Position',[2 574 1915 420])
        hold on
        h = [0 0];
        x0 = 0;
        for st = 1:size(meanSmall,1)
            plot([0 0] + x0, [mini maxi], 'k:')
            plot([stimDur stimDur] + x0, [mini maxi], 'k:')
            plot(respTime([1 end]) + x0, [0 0], 'k:')
            fill([respTime flip(respTime)] + x0, ...
                [meanSmall(st,:)-semSmall(st,:), flip(meanSmall(st,:)+semSmall(st,:))], ...
                'k', 'FaceColor', 'k', ...
                'FaceAlpha', 0.1, 'EdgeColor', 'none')
            h1 = plot(respTime + x0, meanSmall(st,:), 'Color', 'k', 'LineWidth', 2);
            fill([respTime flip(respTime)] + x0, ...
                [meanLarge(st,:)-semLarge(st,:), flip(meanLarge(st,:)+semLarge(st,:))], ...
                'k', 'FaceColor', 'r', ...
                'FaceAlpha', 0.1, 'EdgeColor', 'none')
            h2 = plot(respTime + x0, meanLarge(st,:), ...
                'Color', 'r', 'LineWidth', 2);
            x0 = x0 + traceDur + xDist;
            if st==1
                h(1) = h2;
                h(2) = h1;
            end
        end
        axis tight
        set(gca, 'XTick', [0 stimDur])
        legend(h, {'small pupil','large pupil'})
        xlabel('Stimuli')
        ylabel('\DeltaF/F')
        title(sprintf('Bouton %d',ids(iBout)))
        
        if ismember(ids(iBout), exBoutons)
            print(fullfile(folder, sprintf('Plane%02d_ID%03d.eps', ...
                planes(iBout), ids(iBout))), '-depsc', '-tiff', '-painters')
        end
%         saveas(gcf, fullfile(folder, sprintf('Plane%02d_ID%03d.jpg', ...
%             planes(iBout), ids(iBout))))
        close(gcf)
    end
end

maxResps(all(all(isnan(maxResps),2),3),:,:) = [];
maxResps = maxResps - mean(mean(maxResps(:,baseInds,:),2),3);
maxResps = maxResps ./ max(maxResps(:,respInds,1),[],2);

m = squeeze(mean(maxResps,1));
s = squeeze(std(maxResps,0,1) ./ sqrt(size(maxResps,1)));
mini = min(m(:) - s(:));
maxi = max(m(:) + s(:));
rng = maxi - mini;
mini = mini - 0.05*rng;
maxi = maxi + 0.05*rng;

h = [0 0];
figure
hold on
plot([0 0], [mini maxi], 'k:')
plot([stimDur stimDur], [mini maxi], 'k:')
plot(respTime([1 end]), [0 0], 'k:')
fill([respTime flip(respTime)], [m(:,1) - s(:,1); flip(m(:,1) + s(:,1))]', ...
    'k', 'FaceColor', 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
h(1) = plot(respTime, m(:,1), 'Color', 'k', 'LineWidth', 2);
fill([respTime flip(respTime)], [m(:,2) - s(:,2); flip(m(:,2) + s(:,2))]', ...
    'k', 'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
h(2) = plot(respTime, m(:,2), 'Color', 'r', 'LineWidth', 2);
axis tight
legend(h, {'small pupil','large pupil'})
xlabel('Time from stim onset (in s)')
ylabel('\DeltaF/F (norm.)')
title('Population average (no suppressed boutons)')

print(fullfile(folder, 'Population.eps'), '-depsc', '-tiff', '-painters')
close(gcf)