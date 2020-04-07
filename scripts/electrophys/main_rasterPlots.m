%% Folders
folderData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\electrophys\';
folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\Electrophys_SC\nonVisualEffects\';

%% Datasets
db_ephys_driftingGratings
% iSet = 1; % SS061
iSet = 9; % [2 3 6 7 8 9]

%% Parameters
% rasterWindow = [-2 2];
rasterWindow = [-1 1];
stimOffset = [.1 -.1];

%% Get data
% load data (time, spikeCounts, cellIDs, depths, stimMatrix,
% directions, blanks, stimTimes, laserOn, nSpikes, sampleRate, spikeWidths,
% timeToEphys, waveforms)
load(fullfile(folderData, db(iSet).subject, db(iSet).date, ...
    sprintf('%02d_data.mat', db(iSet).exp)));
notBlanks = ~ismember(1:size(stimMatrix,1), blanks)';
sampleRate = 1 / median(diff(time));

% load tuning structure
data = load(fullfile(folderResults, 'pupil', 'tuning_prefDirSigmaDIFixed_isTuned.mat'));
tuning = data.tuning;

winSamples = round(rasterWindow .* sampleRate);
stimDur = mean((stimTimes.offset+stimOffset(2)) - (stimTimes.onset+stimOffset(1)));
rasterTime = (winSamples(1) : round(stimDur*sampleRate)+winSamples(2)) ./ sampleRate;
responses = ssLocal.getTracesPerStimulus(spikeCounts, stimMatrix, ...
    [-winSamples(1), winSamples(2)]); % [neurons x stimuli x repetitions x time]
stimOffsetSamples = round(stimOffset .* sampleRate);
responses(:,:,:,[1:stimOffsetSamples(1), end+stimOffsetSamples(2)+1:end]) = [];
responses = responses(:, ~laserOn, :, :);

tunedNeurons = find(tuning(iSet).isTuned);
suppressed = tuning(iSet).isSuppressed(tunedNeurons);
responses = responses(tunedNeurons,:,:,:);
modulation = NaN(length(tunedNeurons),1);
prefStims = NaN(length(tunedNeurons),1);
%stim directions
stimInds = find(~laserOn & notBlanks);
dirs = directions(ismember(directions(:,2),stimInds), 1);
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
% enh = find(suppressed==-1);
% [mEnh, sortedEnhanced] = sort(modulation(enh), 'descend');
% prefEnh = prefStims(enh(sortedEnhanced));
% sup = find(suppressed==1);
% [mSup, sortedSuppressed] = sort(modulation(sup), 'descend');
% prefSup = prefStims(sup(sortedSuppressed));
[m, sorted] = sort(modulation, 'descend');


pos = {[1 41 1920 1083], [1921 41 1920 1083]};
labels = {'small', 'large'};

% meanResps = mean(mean(responses(:,:,:,-winSamples(1)+1:end-winSamples(2)),4),3);
% Single neuron rasterplots
for n = 1:length(sorted)
    % j is wrong!!! need to look at tuning.cond.cell.directions, not at
    % directions!
    j = find(stimInds(prefStims(sorted(n)))==directions(:,2));
    smallPupil = ~isnan(tuning(iSet).cond(1).cell(1).responses(j,:));
    conditions = {find(smallPupil), find(~smallPupil)};
    r = NaN(length(smallPupil)+2, size(responses,4));
    i = 0;
    for c = 1:2
        r(i+(1:length(conditions{c})),:) = ...
            responses(sorted(n), prefStims(sorted(n)), conditions{c}, :);
        i = i + length(conditions{c}) + 2;
    end
    r = abs(double(r>0)-1);
    r(isnan(r)) = 1;
    figure('Position', [2030 900 1695 205])
    imagesc(rasterTime([1 end]), [1 size(r,1)], r)
    colormap gray
    hold on
    plot([0 0], [.5 size(r,1)+.5], 'r')
    plot([1 1].*stimDur, [.5 size(r,1)+.5], 'r')
    plot(rasterTime([1 end]), [1 1].*(length(conditions{1})+1.5), 'k:')
    set(gca, 'box', 'off', 'YTick', [])
    xlabel('Time from stimulus onset (in s)')
    title(sprintf('Neuron %d, modulation: %.2f, suppressed: %d', ...
        tunedNeurons(sorted(n)), m(n), suppressed(sorted(n))))
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