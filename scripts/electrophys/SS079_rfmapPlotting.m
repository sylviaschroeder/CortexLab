%% Define dataset
subject = 'SS080';
date = '2017-09-06';
exp = 2;
TLexp = 1;

%% Define folders
protocolFolder = '\\ZSERVER.cortexlab.net\Data\trodes';
hardwareFolder = '\\ZSERVER.cortexlab.net\Data\expInfo';
subjectsFolder = '\\ZSERVER.cortexlab.net\Data\Subjects';

%% load the dataset
sp = loadAllKsDir(subject, date);

%% plot drift map
incl = find(sp.spikeAmps > 40);
incl = incl(1:10:end);
figure
plotDriftmap(sp.st(incl), sp.spikeAmps(incl), sp.spikeDepths(incl));
axis tight
set(gca, 'box', 'off')

%% load stimulus parameters
data = load(fullfile(protocolFolder, subject, date, num2str(exp), ...
    'Protocol.mat'));
pars = data.Protocol;
stimFile = str2func(strtok(pars.xfile, '.'));
% load myScreenInfo
load(fullfile(hardwareFolder, subject, date, num2str(exp), ...
    sprintf('%s_%d_%s_hardwareInfo.mat', date, exp, subject)));
myScreenInfo.windowPtr = NaN;
 
% call x-file to create stimuli
SS = stimFile(myScreenInfo, pars.pars);
stimFrames = cat(3, SS.ImageTextures{:});
 
framesPerImage = pars.pars(6,1);
frameTimes = (0 : size(stimFrames, 3)-1) * framesPerImage / myScreenInfo.FrameRate;
 
rows = size(stimFrames,1);
cols = size(stimFrames,2);

ind = find(stimFrames == 1);
t = ceil(ind / (rows * cols));
ind = mod(ind, (rows * cols));
ind(ind == 0) = rows * cols;
x_wh = ceil(ind / rows);
y_wh = mod(ind, rows);
y_wh(y_wh == 0) = rows;
time_wh = frameTimes(t)';

ind = find(stimFrames == -1);
t = ceil(ind / (rows * cols));
ind = mod(ind, (rows * cols));
ind(ind == 0) = rows * cols;
x_bl = ceil(ind / rows);
y_bl = mod(ind, rows);
y_bl(y_bl == 0) = rows;
time_bl = frameTimes(t)';
 
alignDir = fullfile(subjectsFolder, subject, date, 'alignments');
bTLtoMaster = readNPY(fullfile(alignDir, ...
    sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, sp.name)));
stimOnTL = readNPY(fullfile(alignDir, ...
    sprintf('mpep_%d_onsets_in_timeline_%d.npy', exp, TLexp)));
stimOffTL = readNPY(fullfile(alignDir, ...
    sprintf('mpep_%d_offsets_in_timeline_%d.npy', exp, TLexp)));
 
stimOn = applyCorrection(stimOnTL, bTLtoMaster);
stimOff = applyCorrection(stimOffTL, bTLtoMaster);
    
% stimPosition = [pars.pars(2), pars.pars(3), pars.pars(4), pars.pars(5)] / 10;

stimPositions_wh = [repmat(y_wh, 3, 1) repmat(x_wh, 3, 1)];
stimTimes_wh = reshape(bsxfun(@plus, stimOn', time_wh), [], 1);
stimPositions_bl = [repmat(y_bl, 3, 1) repmat(x_bl, 3, 1)];
stimTimes_bl = reshape(bsxfun(@plus, stimOn', time_bl), [], 1);

%% pick some spikes

% choose only those templates in a part of the probe of interest
inRange = find(sp.templateYpos>1650 & sp.templateYpos<1800)-1; 

% sort them by amplitude, largest first
[~,ii] = sort(sp.tempAmps(inRange+1), 'descend');
inRange = inRange(ii);

% pick one with largest amplitude

for n = 1:length(inRange)
    st = sp.st(sp.clu==inRange(n));
    
    %% compute the RF
    params.makePlots = true;
    params.useSVD = true;
    % params.useSVD = false;
    params.countWindow = [-0.05 0.3];
    params.binSize = 0.01;
    
    rfMap_wh = sparseNoiseRF(st, stimTimes_wh, stimPositions_wh, params);
    set(gcf, 'Position', [1923 2 958 1114])
    annotation('textbox', [0 .95 1 .04], 'String', 'White RF', 'FontSize', 15, ...
        'FontWeight', 'bold', 'LineStyle', 'none', 'HorizontalAlignment', 'center')
    rfMap_bl = sparseNoiseRF(st, stimTimes_bl, stimPositions_bl, params);
    set(gcf, 'Position', [2882 2 958 1114])
    annotation('textbox', [0 .95 1 .04], 'String', 'Black RF', 'FontSize', 15, ...
        'FontWeight', 'bold', 'LineStyle', 'none', 'HorizontalAlignment', 'center')
    
    %% plot the waveform of this template, if desired
    
    templateToPlot = inRange(n)+1;
    
    figure; hold on;
    thisTemp = squeeze(sp.tempsUnW(templateToPlot,:,:));
    [~,peakChan] = max(max(abs(thisTemp),[],1),[],2);
    plotChans = peakChan-15:peakChan+15;
    arrayfun(@(x)plot(sp.xcoords(x)+0.3*(1:size(thisTemp,1))', sp.ycoords(x)+5*thisTemp(:,x), 'k'), plotChans)
    set(gcf, 'Position', [1536 42 381 1074])
    
    pause
    
    close all
end