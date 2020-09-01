clc; clear all; 
addpath(genpath('C:\Users\Flora\Github\npy-matlab'));
addpath(genpath('C:\Users\Flora\Github\spikes'));
addpath(genpath('C:\Users\Flora\Github\kilotrodeRig'));
addpath(genpath('C:\Users\Flora\Documents\Python Scripts\CompetitiveSC_scripts\matlab_readers')); 
% load sync information

mouseName='SS087'; 
thisDate='2017-12-12';

timelineExpNum=2;taskblock=3;sparseNoiseblock=4;passiveblock=5; 

%
root=sprintf('Z:\\%s\\%s',mouseName,thisDate);
% make alf folder within root
alfDir=fullfile(root,'alf');
if ~exist(alfDir, 'dir')
   mkdir(alfDir)
end
%% write alf for all probes using kilotrode rig code 
[tags, hasEphys] = getEphysTags(mouseName,thisDate);
sp = loadAllKsDir(mouseName, thisDate); % this code theoretically does the alignment as well
%%
writeEphysALF(mouseName, thisDate);

%% write numpy of good clusters after phy as well as mua
% no need to do this as clusters.groups.npy already contains this info
% fid=tdfread(fullfile(root,sprintf('ephys_%s\\sorting\\cluster_group.tsv',tags{1})));
% 
% goods=find(fid.group=='g');
% goodclus=fid.cluster_id(goods);
% 
% muas=find(fid.group=='m');
% muaclus=fid.cluster_id(muas);

%% convert sparse noise
iN='sparseNoise'; 
[block,stimArrayTimes]=alignsource(root,mouseName,thisDate,sparseNoiseblock,timelineExpNum);
[stimTimeInds, stimPositions, stimArray] = ...
        computeSparseNoiseSignals(block);
stimTimes = cellfun(@(x)stimArrayTimes(x), stimTimeInds, 'uni', false); 

alf.writeEventseries(alfDir,sprintf('%s',iN), stimTimes{1}, [], []);
writeNPY(stimPositions{1}, fullfile(alfDir, sprintf('%s.positions.npy',iN)));
clear block
clear stimArrayTimes
clear iN


%% convert passive

iN='passive';
[block,~]=alignsource(root,mouseName,thisDate,passiveblock,timelineExpNum);
cond=[block.trial.condition];
% timings
blockFlipsTimes = block.stimWindowUpdateTimes; 
pdFlipTimes=readNPY(fullfile(root,sprintf('alignments\\block_%d_sw_in_timeline_%d.npy',passiveblock,timelineExpNum)));
co=readNPY(fullfile(root,sprintf('alignments\\correct_block_%d_to_timeline_%d.npy',passiveblock,timelineExpNum)));
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

if isfield(block.trial, 'stimulusCueEndedTime')
    s.stimOffTimes = follows(...
        toPDTimeFrameLag([block.trial.stimulusCueEndedTime]), pdFlipTimes);
end

%% fit rest of the times in the block
blockfit=robustfit([block.trial.stimulusCueStartedTime],s.stimOnTimes);
fitblocktimes= @(t)t*blockfit(2) + blockfit(1);
s.trialStartedTime=fitblocktimes([block.trial.trialStartedTime]);
s.onsetToneTime=fitblocktimes([block.trial.onsetToneSoundPlayedTime]); 
s.feedbackOnTime=fitblocktimes([block.trial.feedbackStartedTime]); 
s.feedbackOffTime=fitblocktimes([block.trial.feedbackEndedTime]);

%% select relevant trials 
ix_onset_tone=find([cond.interactiveOnsetToneRelAmp]>0.0001);
ix_positive_feedback=find([block.trial.feedbackType]==1); % this is basically the same sound as the valve 
ix_negative_feedback=find([cond.negFeedbackSoundAmp]>0); % trials with negative feedback
ix_stimulus=find([cond.interactiveOnsetToneRelAmp]<0.0001 & [block.trial.feedbackType]==-1 & [cond.negFeedbackSoundAmp]==0);
%% WRITE NPY
% ITI start
alf.writeEventseries(alfDir,sprintf('%s.visualITI_Start',iN), s.trialStartedTime(ix_stimulus), [], []);
alf.writeEventseries(alfDir,sprintf('%s.onsetToneITI_Start',iN), s.trialStartedTime(ix_onset_tone), [], []);
alf.writeEventseries(alfDir,sprintf('%s.valveITI_Start',iN), s.trialStartedTime(ix_positive_feedback), [], []);
alf.writeEventseries(alfDir,sprintf('%s.whiteNoiseITI_Start',iN), s.trialStartedTime(ix_negative_feedback), [], []);

%visual stimulus info 
contr=[cond.visCueContrast];
stimDist=[cond.distBetweenTargets];
writeNPY(contr(1,ix_stimulus),fullfile(alfDir, sprintf('%s.LeftContrast.npy',iN)));
writeNPY(contr(2,ix_stimulus),fullfile(alfDir, sprintf('%s.RightContrast.npy',iN)));
writeNPY(stimDist(ix_stimulus),fullfile(alfDir,sprintf('%s.stimuliDistance.npy',iN)));
alf.writeEventseries(alfDir,sprintf('%s.stimON',iN), s.stimOnTimes(ix_stimulus) , [], []);
alf.writeEventseries(alfDir,sprintf('%s.stimOFF',iN), s.stimOffTimes(ix_stimulus), [], []);

%audio information from onset tone
onset_tone=[cond.interactiveOnsetToneRelAmp];
writeNPY(onset_tone(ix_onset_tone),fullfile(alfDir,sprintf('%s.onsetTone_amp.npy',iN)));
alf.writeEventseries(alfDir,sprintf('%s.onsetTone',iN), s.onsetToneTime(ix_onset_tone), [], []);

% feedback type (1 when it is valve click -1 is white noise tone)
writeNPY([block.trial.feedbackType], fullfile(alfDir, sprintf('%s.feedbackType.npy',iN)));

% valve click times 
alf.writeEventseries(alfDir,sprintf('%s.valve_clickON',iN), s.feedbackOnTime(ix_positive_feedback), [], []);
alf.writeEventseries(alfDir,sprintf('%s.valve_clickOFF',iN), s.feedbackOffTime(ix_positive_feedback), [], []);

% negative feedback sounds
alf.writeEventseries(alfDir,sprintf('%s.whiteNoiseON',iN), s.feedbackOnTime(ix_negative_feedback), [], []);
alf.writeEventseries(alfDir,sprintf('%s.whiteNoiseOFF',iN), s.feedbackOffTime(ix_negative_feedback), [], []);

%alf.writeEventseries(alfDir,sprintf('%s.feedbackON',iN), s.feedbackOnTime, [], []);
%alf.writeEventseries(alfDir,sprintf('%s.feedbackOFF',iN), s.feedbackOffTime, [], []);
%% %%%%%%%% write choice world %%%%%%%%
iN='ChoiceWorld';
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
blockfit=robustfit([block.trial.stimulusCueStartedTime],s.stimOnTimes);
fitblocktimes= @(t)t*blockfit(2) + blockfit(1);
s.trialStartedTime=fitblocktimes([block.trial.trialStartedTime]);
s.onsetToneTime=reshape(fitblocktimes([block.trial.onsetToneSoundPlayedTime]),2,[]); % this has two arrays _ I think it measures on and offsets!
s.feedbackOnTime=fitblocktimes([block.trial.feedbackStartedTime]); 
s.feedbackOffTime=fitblocktimes([block.trial.feedbackEndedTime]);
s.interactiveStartedTime=fitblocktimes([block.trial.interactiveStartedTime]);
s.interactiveEndedTime=fitblocktimes([block.trial.interactiveStartedTime]);


%visual stimulus info 
contr=[cond.visCueContrast];
writeNPY(contr(1,:),fullfile(alfDir, sprintf('%s.LeftContrast.npy',iN)));
writeNPY(contr(2,:),fullfile(alfDir, sprintf('%s.RightContrast.npy',iN)));
alf.writeEventseries(alfDir,sprintf('%s.stimON',iN), s.stimOnTimes, [], []);
alf.writeEventseries(alfDir,sprintf('%s.stimOFF',iN), s.stimOffTimes, [], []);

alf.writeEventseries(alfDir,sprintf('%s.onsetToneON',iN), s.onsetToneTime(1,:), [], []);
alf.writeEventseries(alfDir,sprintf('%s.onsetToneOFF',iN), s.onsetToneTime(2,:), [], []);

writeNPY([block.trial.feedbackType], fullfile(alfDir, sprintf('%s.feedbackType.npy',iN)));

alf.writeEventseries(alfDir,sprintf('%s.feedbackON',iN), s.feedbackOnTime, [], []);
alf.writeEventseries(alfDir,sprintf('%s.feedbackOFF',iN), s.feedbackOffTime, [], []);

alf.writeEventseries(alfDir,sprintf('%s.interactiveStart',iN), s.interactiveStartedTime, [], []);
alf.writeEventseries(alfDir,sprintf('%s.interactiveEnd',iN), s.interactiveEndedTime, [], []);


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

