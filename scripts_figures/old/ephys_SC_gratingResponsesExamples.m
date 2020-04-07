folderData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\electrophys\';
folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\Electrophys_SC\nonVisualEffects\';

set = 1; % 'SS061', '20160511', exp 1
neuron = 28;

beforeStim = 0.5;
afterStim = 0.5;
winStd = 10;

%% Plot PSTHs of example neuron:
% (1) large pupil with V1, (2) large pupil without V1, (3) small pupil with
% V1, (4) small pupil without V1

% load tuning data
data = load(fullfile(folderResults, 'pupil', 'tuning_prefDirSigmaDIFixed_lineFit.mat'));
tuning = data.tuning;

% load spike data (time, spikeCounts, cellIDs, dephts, stimMatrix,
% directions, blanks, stimTimes, laserOn, laserTimesRelative,
% nSpikes, sampleRate, spikeWidths, stimSequence, timelineToEphys,
% waveforms)
spData = load(fullfile(folderData, tuning(set).subject, tuning(set).date, ...
    sprintf('%02d_data.mat', tuning(set).exp)));

% get spike counts of neuron
sr0 = 1/diff(spData.time([1 2]));
responses = squeeze(ssLocal.getTracesPerStimulus(spData.spikeCounts(:,neuron), ...
    spData.stimMatrix, round([beforeStim afterStim].*sr0))); % [stim x rep x time]
t0 = ((1:size(responses,3))-1-round(beforeStim*sr0)) ./ sr0;

% get preferred stimulus
[~,prefStim] = min(min(abs(tuning(set).cond(1).cell(neuron).parameters(1) - ...
    spData.directions), abs(tuning(set).cond(1).cell(neuron).parameters(1)-360 ...
    - spData.directions)));
prefStims = find(tuning(set).cond(1).cell(neuron).directions == spData.directions(prefStim));

% get responses to pref stim
resp = zeros(2,size(responses,2),size(responses,3));
resp(1,:,:) = responses(spData.directions==spData.directions(prefStim) & ...
    spData.laserOn==1, :, :);
resp(2,:,:) = responses(spData.directions==spData.directions(prefStim) & ...
    spData.laserOn==0, :, :);
resp = reshape(resp,[],size(resp,3));

% get conditions (1: small pupil, no laser, 2: large pupil, no laser, 3:
% small pupil, laser, 4: large pupil, laser)
conditions = NaN(2, size(tuning(set).cond(1).cell(neuron).responses,2));
for c = 1:4
    conditions(~isnan(tuning(set).cond(c).cell(neuron).responses(prefStims,:))) = c;
end
conditions = reshape(conditions,[],1);

figure
hold on
cols = 'krcm';
win = normpdf(-4*winStd:4*winStd,0,winStd);
h = zeros(1,4);
for c = 1:4
    r = mean(resp(conditions==c,:),1);
    r = conv(r, win, 'same') * sr0;
    h(c) = plot(t0, r, cols(c));
end
yLimits = get(gca, 'YLim');
plot([0 0], yLimits, 'k:')
plot([1 1] .* mean(spData.stimTimes.offset-spData.stimTimes.onset), ...
    yLimits, 'k:')
legend(h, {tuning(set).cond.name})
title(sprintf('%s %s, exp %d, neuron %d, resp to pref stim',tuning(set).subject, ...
    tuning(set).date, tuning(set).exp, neuron))
xlabel('Time from stim onset (s)')
ylabel('Firing rate (sp/s)')