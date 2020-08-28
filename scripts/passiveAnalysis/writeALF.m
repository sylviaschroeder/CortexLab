clc; clear all; 
addpath(genpath('C:\Users\Flora\Github\npy-matlab'));
addpath(genpath('C:\Users\Flora\Github\spikes'));
addpath(genpath('C:\Users\Flora\Github\kilotrodeRig'));
addpath(genpath('C:\Users\Flora\Documents\Python Scripts\CompetitiveSC_scripts\matlab_readers')); 
% load sync information

mouseName='SS087'; 
thisDate='2017-12-12';

sparseNoiseblock=4;
passiveblock=5; 
timelineExpNum=2;
%%
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
[block,stimArrayTimes]=alignsource(root,mouseName,thisDate,passiveblock,timelineExpNum);
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


%% select relevant trials 
ix_onset_tone=find([cond.interactiveOnsetToneRelAmp]>0.0001);
ix_positive_feedback=find([block.trial.feedbackType]==1); % this is basically the same sound as the valve 
ix_negative_feedback=find([cond.negFeedbackSoundAmp]>0); % trials with negative feedback
ix_stimulus=find([cond.interactiveOnsetToneRelAmp]<0.0001 & [block.trial.feedbackType]==-1 & [cond.negFeedbackSoundAmp]==0);
%%
if isfield(block.trial, 'stimulusCueStartedTime')
    s.stimOnTimes = follows(...
        toPDTimeFrameLag([block.trial.stimulusCueStartedTime]), pdFlipTimes);
end

if isfield(block.trial, 'stimulusCueEndedTime')
    s.stimOffTimes = follows(...
        toPDTimeFrameLag([block.trial.stimulusCueEndedTime]), pdFlipTimes);
end

if isfield(block.trial, 'onsetToneSoundPlayedTime')
    s.onsetToneTime = follows(...
        toPDTimeFrameLag([block.trial.onsetToneSoundPlayedTime]), pdFlipTimes);
end

if isfield(block.trial, 'feedbackStartedTime')
    s.feedbackOnTime = follows(...
        toPDTimeFrameLag([block.trial.feedbackStartedTime]), pdFlipTimes);
end

if isfield(block.trial, 'feedbackEndedTime')
    s.feedbackOffTime = follows(...
        toPDTimeFrameLag([block.trial.feedbackEndedTime]), pdFlipTimes);
end

%%
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