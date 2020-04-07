%% Folders
folderData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\electrophys\';
folderEye = '\\zserver.cortexlab.net\Data\EyeCamera';
folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\Electrophys_SC\nonVisualEffects\';
folderPlots = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\Electrophys_SC\rasterplots';

%% Datasets
db_ephys_driftingGratings

iSets = [1 2 3 6 7 8 9];
% iSet = 1; % SS061
% iSet = 9; % [2 3 6 7 8 9]

%% Parameters
% rasterWindow = [-2 2];
rasterWindows = [[-2 2]; repmat([-1 1],6,1)];
% stimOffset = [0 0];
stimOffsets = [[0 0]; repmat([.1 -.1],6,1)];

% load tuning structure
data = load(fullfile(folderResults, 'pupil', 'tuning_prefDirSigmaDIFixed_isTuned.mat'));
tuning = data.tuning;


for j = 1:length(iSets)
    iSet = iSets(j);
    %% Get data
    % load data (time, spikeCounts, cellIDs, depths, stimMatrix,
    % directions, blanks, stimTimes, laserOn, nSpikes, sampleRate, spikeWidths,
    % timeToEphys, waveforms)
    load(fullfile(folderData, db(iSet).subject, db(iSet).date, ...
        sprintf('%02d_data.mat', db(iSet).exp)));
    notBlanks = ~ismember(1:size(stimMatrix,1), blanks)';
    sampleRate = 1 / median(diff(time));
    
    % load pupil data
    data = load(fullfile(folderEye, db(iSet).subject, db(iSet).date, ...
        sprintf('%02d_eye_processed.mat', db(iSet).exp)), 'results');
    pupil = data.results;
    nonVisData = nonVis.getPupilDiam(pupil);
    data = load(fullfile(folderEye, db(iSet).subject, db(iSet).date, ...
        sprintf('%02d_eyeTime.mat', db(iSet).exp)));
    nonVisTime = data.eyeTime;
    
    % get responses aligned to stimulus onset
    winSamples = round(rasterWindows(j,:) .* sampleRate);
    stimDur = mean((stimTimes.offset+stimOffsets(j,2)) - (stimTimes.onset+stimOffsets(j,1)));
    rasterTime = (winSamples(1) : round(stimDur*sampleRate)+winSamples(2)) ./ sampleRate;
    responses = ssLocal.getTracesPerStimulus(spikeCounts, stimMatrix, ...
        [-winSamples(1), winSamples(2)]); % [neurons x stimuli x repetitions x time]
    stimOffsetSamples = round(stimOffsets(j,:) .* sampleRate);
    responses(:,:,:,[1:stimOffsetSamples(1), end+stimOffsetSamples(2)+1:end]) = [];
    responses = responses(:, ~laserOn, :, :);

    % get median pupil size for each trial
    nanInds = isnan(nonVisData);
    nonVisual = interp1(nonVisTime(~nanInds), nonVisData(~nanInds), ...
        time, 'pchip');
    nanTimes = nonVisTime(isnan(nonVisData));
    nanTimes = hist(nanTimes, time)>0;
    nonVisual(nanTimes) = NaN;
    nonVisual = squeeze(ssLocal.getTracesPerStimulus( ...
        nonVisual, stimMatrix, [1 0])); % [stimuli x repetitions x time]
    % nonVisual(:,:,[1:stimOffsetSamples(1), end+stimOffsetSamples(2)+1:end]) = [];
    ind = sum(isnan(nonVisual),3) > size(nonVisual,3)/2; % trials where
    % nonvisual signal is NaN more than half of the times
    nonVisual = nanmean(nonVisual, 3); % [stimuli x repetitions]
    nonVisual(ind) = NaN;
    nonVisual = nonVisual(~laserOn,:);

    % get size of modulation due to pupil for each neuron
    tunedNeurons = find(tuning(iSet).isTuned);
    suppressed = tuning(iSet).isSuppressed(tunedNeurons);
    responses = responses(tunedNeurons,:,:,:);
    modulation = NaN(length(tunedNeurons),1);
    prefStims = NaN(length(tunedNeurons),1);
    %stim directions
    stimInds = find(~laserOn & notBlanks);
    dirs = [directions(ismember(directions(:,2),stimInds), 1); 360];
    for k = 1:length(tunedNeurons)
        small = sum(tuning(iSet).cond(1).cell(tunedNeurons(k)).parameters([2 4]));
        large = sum(tuning(iSet).cond(2).cell(tunedNeurons(k)).parameters([2 4]));
        m = min([small, large]);
        if m < 0
            small = small - m;
            large = large - m;
        end
        modulation(k) = (large - small) / (abs(large) + abs(small));
        [~,s] = min(abs(tuning(iSet).cond(1).cell(tunedNeurons(k)).parameters(1)-dirs));
        prefStims(k) = s;
    end
    prefStims(prefStims == length(dirs)) = 1;
    % enh = find(suppressed==-1);
    % [mEnh, sortedEnhanced] = sort(modulation(enh), 'descend');
    % prefEnh = prefStims(enh(sortedEnhanced));
    % sup = find(suppressed==1);
    % [mSup, sortedSuppressed] = sort(modulation(sup), 'descend');
    % prefSup = prefStims(sup(sortedSuppressed));
    [m, sorted] = sort(modulation, 'descend');

    % make plots
    pos = {[1 41 1920 1083], [1921 41 1920 1083]};
    labels = {'small', 'large'};
    
    folder = fullfile(folderPlots, [db(iSet).subject '_' db(iSet).date '_' num2str(db(iSet).exp)]);
    if ~isdir(folder)
        mkdir(folder)
    end
    
    % meanResps = mean(mean(responses(:,:,:,-winSamples(1)+1:end-winSamples(2)),4),3);
    % Single neuron rasterplots
    for n = 1:length(sorted)
        if isnan(m(n))
            continue
        end
        r = permute(responses(sorted(n), prefStims(sorted(n)), :, :), [3 4 1 2]); % [repetitions x time]
        p = nonVisual(prefStims(sorted(n)), :); % [1 x repetitions]
        [pupil, order] = sort(p, 'ascend');
        r = r(order,:);
        r(isnan(pupil),:) = [];
        order(isnan(pupil)) = [];
        pupil(isnan(pupil)) = [];
        r = abs(double(r>0)-1);
        r(isnan(r)) = 1;
        
        figure('Position', [2030 900 1695 205])
        subplot(1, 12, 1:11)
        imagesc(rasterTime([1 end]), [1 size(r,1)], r)
        colormap gray
        hold on
        plot([0 0], [.5 size(r,1)+.5], 'r')
        plot([1 1].*stimDur, [.5 size(r,1)+.5], 'r')
        set(gca, 'box', 'off', 'YTick', [])
        xlabel('Time from stimulus onset (in s)')
        title(sprintf('Neuron %d, modulation: %.2f, suppressed: %d', ...
            tunedNeurons(sorted(n)), m(n), suppressed(sorted(n))))
        ax1 = gca;
        subplot(1, 12 , 12)
        imagesc(pupil(:))
        set(gca, 'XTick', [], 'YTick', [1 length(pupil)], 'YTickLabel', ...
            round(pupil([1 end]).*10)./10)
        title('Pupil diameter')
        pos = get(gca, 'Position');
        set(gca, 'Position', [pos(1) 0.22 .02 .65])
        set(ax1, 'Position', [.1 .22 .7 .65])
        
        f = gcf;
        f.PaperPositionMode = 'auto';
        print(fullfile(folder, sprintf('%02d_cell%03d.tiff', n, tunedNeurons(sorted(n)))), '-dtiff', '-r0')
        close gcf
    end
end

% Population rasterplots
% for st =1:size(responses,2)
%     j = find(stimInds(st)==directions(:,2));
%     if isempty(j)
%         continue
%     end
%     smallPupil = ~isnan(tuning(iSet).cond(1).cell(1).responses(j,:));
%     conditions = {find(smallPupil), find(~smallPupil)};
%     for c = 1:2
%         figure('Position', pos{c})
%         colormap gray
%         for rep = 1:length(conditions{c})
%             subplot(2, length(conditions{c}), rep)
%             r = double(permute(responses(sorted(supprSorted==-1), st, ...
%                 conditions{c}(rep), :), [1 4 2 3]) > 0);
%             r = abs(r - 1);
%             imagesc(rasterTime([1 end]), [1 size(r,1)], r)
%             hold on
%             plot([0 0], [.5 size(r,1)+.5],'r')
%             plot([stimDur stimDur], [.5 size(r,1)+.5],'r')
%             
%             subplot(2, length(conditions{c}), length(conditions{c})+rep)
%             r = double(permute(responses(sorted(supprSorted==1), st, ...
%                 conditions{c}(rep), :), [1 4 2 3]) > 0);
%             r = abs(r - 1);
%             imagesc(rasterTime([1 end]), [1 size(r,1)], r)
%             hold on
%             plot([0 0], [.5 size(r,1)+.5],'r')
%             plot([stimDur stimDur], [.5 size(r,1)+.5],'r')
%         end
%         annotation('textbox', [0 0.96 1 0.03], 'String', ...
%         sprintf('Stimulus %d, %s pupil', st, labels{c}), ...
%         'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
%         'HorizontalAlignment', 'center')
%     end
% end