%% Folders
folderTools = 'C:\Users\Flora\Github';
%folderTools = 'C:\STORAGE\workspaces';
folderData = 'Z:';
%folderData = '\\zubjects.cortexlab.net\Subjects';
folderScript = 'CortexLab\scripts\passiveAnalysis';

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')));
addpath(genpath(fullfile(folderTools, 'spikes')));
addpath(genpath(fullfile(folderTools, 'kilotrodeRig')));
addpath(genpath(fullfile(folderTools, 'alyx-matlab')));
addpath(genpath(fullfile(folderTools, folderScript)));

%% Define dataset and prepare folders
mouseName='SS087'; 
thisDate='2017-12-12';

timelineExpNum=2;
taskblock=3;
sparseNoiseblock=4;
passiveblock=5; 

root = fullfile(folderData, mouseName, thisDate);

% make alf folder within root
alfDir = fullfile(root, 'alf');
if ~exist(alfDir, 'dir')
   mkdir(alfDir)
end

% make align folder
alignDir = fullfile(root, 'alignments');

%% Basic info
[tags, hasEphys] = getEphysTags(mouseName,thisDate);
bTLtoMaster = readNPY(fullfile(alignDir, ...
    sprintf('correct_timeline_%d_to_ephys_%s.npy', timelineExpNum, tags{1})));

%% write alf for all probes using kilotrode rig code 
% SS: sp is never used! but will be usefull if we want to save spikes
% differently (all probes in the same folder)
% sp = loadAllKsDir(mouseName, thisDate); % this code theoretically does the alignment as well

writeEphysALF(mouseName, thisDate);

%% convert sparse noise
iN = 'sparseNoise'; 
[block, stimArrayTimes] = alignsource(root, mouseName, thisDate, ...
    sparseNoiseblock, timelineExpNum);
% NOTE: computeSparseNoiseSignals has the stimulus position hard coded. It
% should however be read out from the experiment definition, because this
% may change! I wouldn't trust the exact position now. I will change this
% in the future.
[stimTimeInds, stimPositions, stimArray] = ...
        computeSparseNoiseSignals(block);
stimTimes = cellfun(@(x)stimArrayTimes(x), stimTimeInds, 'uni', false);

alf.writeEventseries(alfDir, sprintf('_ibl_%s',iN), stimTimes{1}, [], []);
% SS: why only use stimPositions{1}? stimPositions are given in [y,x], but
% IBL standard is [x,y]!
writeNPY(stimPositions{1}, fullfile(alfDir, sprintf('_ibl_%s.xy.npy',iN)));
clear block
clear stimArrayTimes
clear iN

%% convert passive
iN = '_ss_passiveTrials';
[block,~] = alignsource(root,mouseName,thisDate,passiveblock,timelineExpNum);
cond = [block.trial.condition];
% timings
blockFlipsTimes = block.stimWindowUpdateTimes; 
pdFlipTimes = readNPY(fullfile(root,sprintf('alignments\\block_%d_sw_in_timeline_%d.npy',passiveblock,timelineExpNum)));
co = readNPY(fullfile(root,sprintf('alignments\\correct_block_%d_to_timeline_%d.npy',passiveblock,timelineExpNum)));
co = co';
s.coeff = co';
s.blockToPdTimeFrame = @(t)t*co(1) + co(2);
s.pdToBlockTImeFrame = @(t)(t - co(2))/co(1);
%offset to place every block flip after corresponding photodiode flip
lag = -max(blockFlipsTimes*co(1) - pdFlipTimes);
toPDTimeFrameLag = @(t)t*co(1) + lag;

%% 
if isfield(block.trial, 'stimulusCueStartedTime')
    s.stimOnTimes = follows(...
        toPDTimeFrameLag([block.trial.stimulusCueStartedTime]), pdFlipTimes);
    s.stimOnTimes = applyCorrection(s.stimOnTimes, bTLtoMaster);
end

if isfield(block.trial, 'stimulusCueEndedTime')
    s.stimOffTimes = follows(...
        toPDTimeFrameLag([block.trial.stimulusCueEndedTime]), pdFlipTimes);
    s.stimOffTimes = applyCorrection(s.stimOffTimes, bTLtoMaster);
end

%% fit rest of the times in the block
blockfit = robustfit([block.trial.stimulusCueStartedTime], s.stimOnTimes);
fitblocktimes = @(t)t*blockfit(2) + blockfit(1);

s.trialStartedTime=fitblocktimes([block.trial.trialStartedTime]);
s.onsetToneTime=reshape(fitblocktimes([block.trial.onsetToneSoundPlayedTime]),2,[]); % this has two arrays _ I think it measures on and offsets!
s.feedbackOnTime=fitblocktimes([block.trial.feedbackStartedTime]); 
s.feedbackOffTime=fitblocktimes([block.trial.feedbackEndedTime]);


%% select relevant trials 
%ix_onset_tone = find([cond.interactiveOnsetToneRelAmp]>0.0001);
%ix_positive_feedback = find([block.trial.feedbackType]==1); % this is basically the same sound as the valve 
%ix_negative_feedback = find([cond.negFeedbackSoundAmp]>0); % trials with negative feedback
%ix_stimulus = find([cond.interactiveOnsetToneRelAmp]<0.0001 & [block.trial.feedbackType]==-1 & [cond.negFeedbackSoundAmp]==0);

%% WRITE NPY
% ITI start
alf.writeEventseries(alfDir,sprintf('%s',iN), s.trialStartedTime, [], []);
%visual stimulus info 
contr=[cond.visCueContrast];
stimDist=[cond.distBetweenTargets];
writeNPY(contr(1,:),fullfile(alfDir, sprintf('%s.contrastLeft.npy',iN)));
writeNPY(contr(2,:),fullfile(alfDir, sprintf('%s.contrastRight.npy',iN)));
writeNPY(stimDist,fullfile(alfDir,sprintf('%s.stimuliDistance.npy',iN)));

%%
alf.writeInterval(alfDir,sprintf('%s.stimOn',iN),s.stimOnTimes,s.stimOffTimes, [],[]);
%%
%audio information from onset tone
onset_tone=[cond.interactiveOnsetToneRelAmp];
writeNPY(onset_tone,fullfile(alfDir,sprintf('%s.goCue_amp.npy',iN)));
alf.writeInterval(alfDir,sprintf('%s.goCue',iN), s.onsetToneTime(1,:), s.onsetToneTime(2,:), [], []);

% feedback type (1 when it is valve click -1 is white noise tone)
writeNPY([block.trial.feedbackType], fullfile(alfDir, sprintf('%s.feedbackType.npy',iN)));
alf.writeInterval(alfDir,sprintf('%s.feedbackOn',iN), s.feedbackOnTime,s.feedbackOffTime, [], []);
% white noise amplitude
writeNPY([cond.negFeedbackSoundAmp], fullfile(alfDir, sprintf('%s.feedback_amp.npy',iN)));

%% %%%%%%%% write choice world %%%%%%%%
iN='_ss_ChoiceWorldTrials';
%iNIBL='_ibl_trials';
[block,~]=alignsource(root,mouseName,thisDate,taskblock,timelineExpNum);
cond=[block.trial.condition];
% timings
blockFlipsTimes = block.stimWindowUpdateTimes; 
pdFlipTimes=readNPY(fullfile(root,sprintf('alignments\\block_%d_sw_in_timeline_%d.npy',taskblock,timelineExpNum)));
co=readNPY(fullfile(root,sprintf('alignments\\correct_block_%d_to_timeline_%d.npy',taskblock,timelineExpNum)));
co=co';
s.coeff = co';
s.blockToPdTimeFrame = @(t)t*co(1) + co(2);
s.pdToBlockTImeFrame = @(t)(t - co(2))/co(1);
%offset to place every block flip after corresponding photodiode flip
lag = -max(blockFlipsTimes*co(1) - pdFlipTimes);
toPDTimeFrameLag = @(t)t*co(1) + lag;

%% 
if isfield(block.trial, 'stimulusCueStartedTime')
    s.stimOnTimes = follows(...
        toPDTimeFrameLag([block.trial.stimulusCueStartedTime]), pdFlipTimes);
end

if isfield(block.trial, 'stimulusCueEndedTime') % the last trial is not always completed, hence off times can be shorter array than on times
    s.stimOffTimes = follows(...
        toPDTimeFrameLag([block.trial.stimulusCueEndedTime]), pdFlipTimes);
end

%% fit rest of the times in the block

blockfit = robustfit([block.trial.stimulusCueStartedTime],s.stimOnTimes);
fitblocktimes = @(t)t*blockfit(2) + blockfit(1);
s.trialStartedTime=fitblocktimes([block.trial.trialStartedTime]);
s.onsetToneTime=reshape(fitblocktimes([block.trial.onsetToneSoundPlayedTime]),2,[]); % this has two arrays _ I think it measures on and offsets!
s.feedbackOnTime=fitblocktimes([block.trial.feedbackStartedTime]); 
s.feedbackOffTime=fitblocktimes([block.trial.feedbackEndedTime]);
s.interactiveStartedTime=fitblocktimes([block.trial.interactiveStartedTime]);
s.interactiveEndedTime=fitblocktimes([block.trial.interactiveStartedTime]);


%visual stimulus info 
contr=[cond.visCueContrast];
writeNPY(contr(1,:),fullfile(alfDir, sprintf('%s.contrastLeft.npy',iN)));
writeNPY(contr(2,:),fullfile(alfDir, sprintf('%s.contrastRight.npy',iN)));
alf.writeInterval(alfDir,sprintf('%s.stimOn',iN),s.stimOnTimes,s.stimOffTimes, [],[]);
alf.writeInterval(alfDir,sprintf('%s.goCue',iN), s.onsetToneTime(1,:), s.onsetToneTime(2,:), [], []);
writeNPY([block.trial.feedbackType], fullfile(alfDir, sprintf('%s.feedbackType.npy',iN)));

alf.writeInterval(alfDir,sprintf('%s.feedbackOn',iN), s.feedbackOnTime,s.feedbackOffTime, [], []);


alf.writeEventseries(alfDir,sprintf('%s.interactiveStart',iN), s.interactiveStartedTime, [], []);

alf.writeEventseries(alfDir,sprintf('%s',iN), s.trialStartedTime, [], []); % that is supposed to be the trial start times 


%% %%%% FUNCTIONS %%%%%%%%
function [block,stimArrayTimes]=alignsource(root,mouseName,thisDate,rfExpNum,tlExpNum)
    alignDir = fullfile(root, 'alignments');
    load(fullfile(root,sprintf('%d\\%s_%d_%s_Block.mat',rfExpNum,thisDate,rfExpNum,mouseName)));
    % read the file that contants linear shift parameters for time all alignments 
    tlToMasterFile = dir(fullfile(alignDir, ...
        sprintf('correct_timeline_%d_to_ephys_*.npy', tlExpNum)));
    bTLtoMaster = readNPY(fullfile(alignDir,tlToMasterFile.name));
    % align timeline for that experiment
    stimArrayTimes = readNPY(fullfile(alignDir, ...
        sprintf('block_%d_sw_in_timeline_%d.npy', rfExpNum, tlExpNum)));

    stimArrayTimes = applyCorrection(stimArrayTimes, bTLtoMaster);
end 

function t = follows(a, b)
    n = numel(a);
    t = zeros(size(a));
    ti = t;
    for ii = 1:n
        ti(ii) = find(b > a(ii), 1);
        t(ii) = b(ti(ii));
    end

    d = t - a;
    range = max(d) - min(d);
    assert((range/mean(d)) < 4.2, 'delta range is much larger than the mean');
end
