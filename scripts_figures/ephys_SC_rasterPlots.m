%% Define data
subject = 'SS061';
date = '2016-05-11';
exp = 2;
dataset = 1;

units =  [2, 4, 7, 34, 35, 51, 95];
exStim = [0, 0, 0,  3,  0, -4,  0];

%% Folders
folderData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\electrophys\';
folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\Electrophys_SC\nonVisualEffects\';

%% Parameters
% rasterWindow = [-2 2];
rasterWindow = [-1 1];
laserConds = [false, true];
smoothStd = 0.1;

%% Get data
% load data (time, spikeCounts, cellIDs, depths, stimMatrix,
% directions, blanks, stimTimes, laserOn, nSpikes, sampleRate, spikeWidths,
% timeToEphys, waveforms)
load(fullfile(folderData, subject, date, sprintf('%02d_data.mat', exp)));
sampleRate = 1 / median(diff(time));

% load tuning structure
data = load(fullfile(folderResults, 'pupil', 'tuning_prefDirSigmaDIFixed_isTuned.mat'));
tuning = data.tuning;

winSamples = round(rasterWindow .* sampleRate);
stimDur = mean(stimTimes.offset - stimTimes.onset);
responses = ssLocal.getTracesPerStimulus(spikeCounts, stimMatrix, ...
    [-winSamples(1), winSamples(2)]); % [neurons x stimuli x repetitions x time]
rasterTime = ((0:size(responses,4)-1) + winSamples(1)) ./ sampleRate;

win = normpdf(-4*smoothStd*sampleRate : 4*smoothStd*sampleRate, ...
    0, smoothStd*sampleRate);
w = round(0.5*sampleRate);

%% Plot PSTHs
cols = 'krcm';
for k = 1:length(units)
    rasters = cell(1, 4);
    resp = squeeze(responses(units(k),:,:,:));
    [~,pref] = min(abs(tuning(dataset).cond(1).cell(units(k)).parameters(1) - ...
        [directions(1:length(directions)/2-1);360]));
    if pref == length(directions)/2 % if closest to 360, set pref stim to 0
        pref = 1;
    end
    pref = pref + exStim(k);
    pref = [pref, pref + length(directions)/2];
    for l = 1:2
        stim = pref(tuning(dataset).laserOn(pref)==laserConds(l));
        for b = 1:2
            rasters{b+(l-1)*2} = permute(resp(stim, ...
                tuning(dataset).nonVisual(stim,:)==b,:), [2 3 1]);
        end
    end
    
    figure
    ax = [0 0];
    subplot(3,1,1:2)
    hold on
    y = 0;
    for c = 1:4
        for rep = 1:size(rasters{c},1)
            t = rasterTime(rasters{c}(rep,:)>0);
            plot([repmat(t,2,1);NaN(1,length(t))], ...
                repmat([y+1 y NaN]',1,length(t)), cols(c))
            y = y+1;
        end
        y = y+1;
    end
    plot([0 0], [0 y], 'k')
    plot([stimDur stimDur], [0 y], 'k')
    xlim(rasterTime([1 end]))
    ylim([0 y])
    set(gca, 'YDir', 'reverse', 'YTick', [])
    title(sprintf('%s, %s, exp %d, unit: %d', subject, date, exp, units(k)))
    ax(1) = gca;
    
    subplot(3,1,3)
    hold on
    maxi = 0;
    for c = 1:4
        r = mean(rasters{c},1) .* sampleRate;
        r = conv([ones(1,w).*mean(r(1:w)), r, ones(1,w).*mean(r(end-w))], ...
            win, 'same');
        r([1:w,end-w+1:end]) = [];
        plot(rasterTime, r, cols(c))
        maxi = max([maxi,r]);
    end
    plot([0 0],[0 maxi],'k')
    plot([stimDur stimDur],[0 maxi], 'k')
    ylim([0 maxi])
    ax(2) = gca;
    
    linkaxes(ax, 'x')
    xlim(rasterTime([1 end]))
    xlabel('Time from stimulus onset (s)')
    ylabel('Firing rate (spikes/s')
    leg = legend('--','+-','-+','++','Location','NorthWest');
    leg.Title.String = 'Laser/Pupil';
end