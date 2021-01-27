%% Parameters
samplingRate = 2500; % in ....imec.lf.meta file
numChans = 385;
window = [-.05 .2];

%% Folders
folderData = '\\zubjects.cortexlab.net\Subjects';
% folderPlots = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\task_ephys\PSTHs';

folderScript = 'C:\dev\workspace\CortexLab';
folderTools = 'C:\STORAGE\workspaces';

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')));
addpath(genpath(fullfile(folderTools, 'spikes')));
addpath(genpath(fullfile(folderScript)));

%% Define datasets
subject = 'SS093';
date = '2018-05-24';
probe = 'K1';

%% Load stimulus data
folderAlf = fullfile(folderData, subject, date, 'alf');
noise = io.getVisNoiseInfo(folderAlf);
noiseStim = (noise.frames(:,:,noise.stimOrder) + 1) ./ 2;
noiseDeriv = diff(cat(3, zeros(size(noiseStim,1), size(noiseStim,2)), ...
    noiseStim), 1, 3); % 1 if switching to white from black, -1 if switching to black from white, 0, if value does not change
noiseDeriv = reshape(noiseDeriv, [], size(noiseStim,3));

%% Load LFP data
folderAlign = fullfile(folderData, subject, date, 'alignments');
folderLFP = fullfile(folderData, subject, date, sprintf('ephys_%s', probe));
fileLFP = dir(fullfile(folderLFP, '*.lf.bin'));

lfp = memmapfile(fullfile(folderLFP, fileLFP.name), ...
    'Format',  {'int16', [numChans fileLFP.bytes/(numChans*2)], 'x'});
t_lfp = (1:fileLFP.bytes/(numChans*2)) ./ samplingRate;

% correct time to align with master if necessary
fileAlign = dir(fullfile(folderAlign, sprintf('correct_ephys_%s_*.npy', probe)));
if ~isempyt(fileAlign)
    probeToMaster = readNPY(fullfile(folderAlign, fileAlign.name));
    t_lfp = applyCorrection(t_lfp, probeToMaster);
end

start = find(t_lfp > noise.times(1)-5, 1);
stop = find(t_lfp > noise.times(end)+1, 1);
lfpMedian = median(lfp.Data.x(:,start:stop), 2);

% find reference/noise channels
lfpStd = std(double(lfp.Data.x(:,start+(0:100000))), 0, 2);
[stdSorted,order] = sort(lfpStd);
[~,jumpInStd] = max(diff(stdSorted(1:100)));
noiseChans = sort(order(1:jumpInStd));
goodChans = setdiff(1:numChans, noiseChans);

%% Get stimulus triggered LFP
winSamples = window(1)*samplingRate : window(2)*samplingRate;
lfpEvoked = zeros(numChans, size(noiseDeriv,1), length(winSamples));

for pix = 1:size(noiseDeriv,1)
    stimEv = noiseDeriv(pix,:) ~= 0;
    t = noise.times(stimEv);
    t_ind = zeros(size(t));
    for j = 1:length(t)
        a = find(t_lfp > t(j), 1);
        if abs(t_lfp(a-1) - t(j)) < abs(t_lfp(a) - t(j))
            a = a-1;
        end
        t_ind(j) = a;
    end
    t = t_ind + winSamples;
    ev = lfp.Data.x(:, t) - lfpMedian;
    ev = reshape(ev, numChans, size(t,1), size(t,2));
    lfpEvoked(:,pix,:) = mean(double(ev) .* noiseDeriv(pix,stimEv), 2);
end

[~, maxChan] = max(max(max(abs(lfpEvoked), [], 3), [], 2), [], 1);
[~, maxT] = max(max(abs(lfpEvoked(maxChan,:,:)), [], 2), [], 3);
[~, maxPix] = max(abs(lfpEvoked(maxChan,:,maxT)), [], 2);
maxRF = reshape(squeeze(lfpEvoked(maxChan,:,maxT)), size(noiseStim,1), size(noiseStim,2));
maxTimeCourse = smooth(squeeze(lfpEvoked(maxChan, maxPix, :)), 20);
[~, maxT] = max(abs(maxTimeCourse));
ampProfile = smooth(lfpEvoked(goodChans,maxPix,maxT),5);

figure
plot(goodChans, ampProfile(goodChans), 'k')

maxTimeCourses = squeeze(lfpEvoked(:, maxPix, :));