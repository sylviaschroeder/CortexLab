
% script to load data for 


%% file locations

mouseName = 'SS074';
thisDate = '2017-01-05';
ephysTag = 'SC';
tlExpNum = 2;
orienExpNum = 3;
masterTimebase = 'SC';

ksDir = fullfile('\\basket.cortexlab.net\data\nick\', mouseName, thisDate, ['ephys_' ephysTag]);

% should be \\zserver\Data\Subjects\etc
expRoot = fileparts(dat.expPath(mouseName, thisDate, 1, 'main', 'master'));
rawDir = fullfile(expRoot, ['ephys_' ephysTag]);

% both of these should be \\zserver\Data\expInfo
tlFile = dat.expFilePath(mouseName, thisDate, tlExpNum, 'timeline', 'master');
protocolFile = dat.expFilePath(mouseName, thisDate, orienExpNum, 'parameters', 'master');

alignDir = fullfile(expRoot, 'alignments');

%% load sync information  


if ~strcmp(ephysTag, masterTimebase)
    bEphysToMaster = readNPY(fullfile(alignDir, ...
        sprintf('correct_ephys_%s_to_ephys_%s.npy', ephysTag, masterTimebase)));
else % this one is master, so use a dummy conversion
    bEphysToMaster = [1; 0];
end

bTLtoMaster = readNPY(fullfile(alignDir, ...
    sprintf('correct_timeline_%d_to_ephys_%s.npy', tlExpNum, masterTimebase)));

stimOnTL = readNPY(fullfile(alignDir, ...
    sprintf('mpep_%d_onsets_in_timeline_%d.npy', orienExpNum, tlExpNum)));
    
stimOffTL = readNPY(fullfile(alignDir, ...
    sprintf('mpep_%d_offsets_in_timeline_%d.npy', orienExpNum, tlExpNum)));

stimOn = applyCorrection(stimOnTL, bTLtoMaster);
stimOff = applyCorrection(stimOffTL, bTLtoMaster);


%% load protocol information

load(protocolFile);
Protocol = parameters.Protocol;

stimIDs = zeros(1, numel(Protocol.seqnums));
for q = 1:size(Protocol.seqnums,1)
    stimIDs(Protocol.seqnums(q,:)) = q;
end

ori = Protocol.pars(strcmp(Protocol.parnames, 'ori'),stimIDs);
laserAmp = Protocol.pars(strcmp(Protocol.parnames, 'amp1'),stimIDs);
contrast = Protocol.pars(strcmp(Protocol.parnames, 'cr'),stimIDs);

%% load spikes

s = loadKSdir(ksDir);

s.st = applyCorrection(s.st, bEphysToMaster);

%% compute basic things from spikes
[spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] = ...
    templatePositionsAmplitudes(s.temps, s.winv, s.ycoords, s.spikeTemplates, s.tempScalingAmps);

spikeDurs = templateDuration(s.spikeTemplates+1);

%% make some psth's
inclSpikes = spikeAmps>50; 
cluByDepth = ceil(spikeDepths/80)*80;

inclTrials = laserAmp==0;

%for convenience in this function, zero contrast trials will be labeled as
%having orientation of -30 degrees
oriC = ori; oriC(contrast==0) = -30; 

trGroups = oriC(inclTrials);
evts = stimOn(inclTrials);

psthViewer(s.st(inclSpikes), cluByDepth(inclSpikes), evts, [-0.25 1.25], trGroups);

%% psth across depth

eventTimes = stimOn; eventName = 'stimOnset';
inclSpikes = spikeAmps>50;

win = [-0.25 1.25];
timeBinSize = 0.005;

bslWin = [-0.25 0];
depthBinSize = 40;

f = figure; 

[timeBins, depthBins, allP] = psthByDepth(s.st(inclSpikes), ...
    spikeDepths(inclSpikes), depthBinSize, timeBinSize, eventTimes, win, bslWin);
plotPSTHbyDepth(timeBins, depthBins, allP, eventName, 'norm', [])
