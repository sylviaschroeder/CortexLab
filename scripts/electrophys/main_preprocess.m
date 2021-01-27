%% Define dataset
% db_ephys_opticTract
k=1;
db(1).subject = 'SS093';
db(1).date = '2018-05-28';

subject = db(k).subject;
date = db(k).date;

%% Define folders
folderTools = 'C:\STORAGE\workspaces';
folderScript = 'C:\dev\workspace\CortexLab';
subjectsFolder = '\\zubjects.cortexlab.net\Subjects';
% subjectsFolder = '\\basket.cortexlab.net\data\nick';
% subjectsFolder = 'J:\Ephys';
alignDir = fullfile(subjectsFolder, subject, date, 'alignments');
% driftPlotFolder = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\OpticTract\driftPlots';
% driftPlotFolder = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\Electrophys_SC\driftPlots';
% spikeDepthFolder = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\Electrophys_SC\spikeDepths';
% waveformFolder = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\Electrophys_SC\spikeWaveForms';

if ~isfolder(alignDir)
    mkdir(alignDir)
end

%% Add paths
addpath('C:\STORAGE\workspaces\Rigbox')
addpath('C:\STORAGE\workspaces\Rigbox\cortexlab')
addpath(genpath('C:\STORAGE\workspaces\Rigbox\cb-tools'))
addpath(genpath(fullfile(folderTools, 'spikes')));
addpath(genpath(fullfile(folderTools, 'npy-matlab')));
addpath(genpath(fullfile(folderTools, 'kilotrodeRig')));
addpath(genpath(fullfile(folderScript)));

%% Get basic info

% if ~exist(alignDir, 'dir')
%     mkdir(alignDir);
% end
[tags, hasEphys] = getEphysTags(subject, date);
% determine what exp nums exist
[expNums, blocks, hasBlock, pars, isMpep, tl, hasTimeline] = ...
    dat.whichExpNums(subject, date);
TLexp = expNums(hasTimeline);
TLexp = TLexp(end);
useFlipper = true;
% useFlipper = false;

%% align times
% 
% % for any ephys, load the sync data
% if hasEphys
%     for t = 1:length(tags)
%         if isempty(tags{t})
%             [~, pdFlips, allET] = loadSync(subject, date);
%         else
%             [~, pdFlips, allET] = loadSync(subject, date, tags{t});
%         end
%         if useFlipper
%             ephysFlips{t} = allET{5}{1};
%         else
%             ephysFlips{t} = pdFlips;
%         end
%     end
% end
% 
% % synchronize multiple ephys to each other
% if hasEphys
%     if length(tags)>1
%         for t2 = 2:length(tags)
%             fprintf(1, 'correct ephys %s to %s\n', tags{t2}, tags{1});
%             [~, b] = makeCorrection(ephysFlips{1}, ephysFlips{t2}, false);
%             writeNPY(b, fullfile(alignDir, sprintf('correct_ephys_%s_to_ephys_%s.npy', tags{t2}, tags{1})));
%         end
%     end
% end
% 
% % detect sync events from timelines
% tlFlips = {};
% for e = 1:length(expNums)
%     if hasTimeline(e)
%         Timeline = tl{e};
%         tt = Timeline.rawDAQTimestamps;
%         if useFlipper
%             evTrace = Timeline.rawDAQData(:, strcmp({Timeline.hw.inputs.name}, 'flipper'));
%             evT = schmittTimes(tt, evTrace, [3 4]); % all flips, both up and down
%         else
%             evTrace = Timeline.rawDAQData(:, strcmp({Timeline.hw.inputs.name}, 'photoDiode'));
%             evT = schmittTimes(tt, evTrace, [3 4]); % all flips, both up and down
%             evT = evT([true; diff(evT)>0.2]);
%         end
%         tlFlips{e} = evT;
%     end
% end
% 
% % match up ephys and timeline events:
% % algorithm here is to go through each timeline available, figure out
% % whether the events in timeline align with any of those in the ephys. If
% % so, we have a conversion of events in that timeline into ephys
% %
% % Only align to the first ephys recording, since the other ones are aligned
% % to that 
% if hasEphys
%     ef = ephysFlips{1};
%     if useFlipper && ef(1)<0.001
%         % this happens when the flipper was in the high state to begin with
%         % - a turning on event is registered at the first sample. But here
%         % we need to drop it. 
%         ef = ef(2:end);
%     end
%     for e = 1:length(expNums)
%         if hasTimeline(e)
%             fprintf('trying to correct timeline %d to ephys\n', expNums(e));
%             %Timeline = tl{e};
%             tlT = tlFlips{e};
%             
%             success=false;
%             if length(tlT)==length(ef)
%                 % easy case: the two are exactly coextensive
%                 [~,b] = makeCorrection(ef, tlT, false);
%                 success = true;
%             elseif length(tlT)<length(ef) && ~isempty(tlT)
%                 [~,b,success] = findCorrection(ef, tlT, false);
%             elseif length(tlT)>length(ef) && ~isempty(tlT)
%                 [~,a,success] = findCorrection(tlT, ef, false);
%                 b = [1/a(1); -a(2)/a(1)];
%             end
%             if success
% %                 writeNPY(b, fullfile(alignDir, ...
% %                     sprintf('correct_timeline_%d_to_ephys_%s.npy', ...
% %                     e, tags{1})));
%                 writeNPY(b, fullfile(alignDir, ...
%                     sprintf('correct_timeline_%d_to_ephys_%s.npy', ...
%                     expNums(e), tags{1})));
%                 fprintf('success\n');
%             else
%                 fprintf('could not correct timeline to ephys\n');
%             end
%         end
%     end
% end
% 
% % match up blocks and mpeps to timeline in order:
% % want to connect each block or mpep with part of a timeline. So go through
% % each of these in order, looking through the timelines in sequence (of
% % what hasn't already been matched) looking for a match. 
% lastTimes = zeros(1,length(expNums));
% for e = 1:length(expNums)
%     if hasBlock(e)
%         for eTL = 1:length(expNums)
%             if hasTimeline(eTL)
%                 fprintf('trying to correct block %d to timeline %d\n', expNums(e), expNums(eTL));
%                 if useFlipper
%                     % didn't get photodiode flips above, so get them now
%                     Timeline = tl{eTL};
%                     tt = Timeline.rawDAQTimestamps;
%                     evTrace = Timeline.rawDAQData(:, strcmp({Timeline.hw.inputs.name}, 'photoDiode'));
%                     pdT = schmittTimes(tt, evTrace, [3 4]);
%                 else
%                     pdT = tlFlips{eTL};
%                 end
%                 block = blocks{e};
%                 sw = block.stimWindowUpdateTimes; 
% %                 sw = sw(2:end); % sometimes need this? Why? how did sw
%                 % get an extra event at the beginning? 
%                 
%                 success = false;
%                 if length(sw)<=length(pdT) && length(sw)>1
%                     [~,b,success,actualTimes] = findCorrection(pdT, sw, false);
%                 end
%                 if success
% %                     writeNPY(b, fullfile(alignDir, ...
% %                         sprintf('correct_block_%d_to_timeline_%d.npy', ...
% %                         e, eTL)));
% %                     writeNPY(actualTimes, fullfile(alignDir, ...
% %                         sprintf('block_%d_sw_in_timeline_%d.npy', ...
% %                         e, eTL)));
%                     writeNPY(b, fullfile(alignDir, ...
%                         sprintf('correct_block_%d_to_timeline_%d.npy', ...
%                         expNums(e), expNums(eTL))));
%                     writeNPY(actualTimes, fullfile(alignDir, ...
%                         sprintf('block_%d_sw_in_timeline_%d.npy', ...
%                         expNums(e), expNums(eTL))));
%                     fprintf('  success\n');
%                     lastTimes(eTL) = actualTimes(end);
%                 else
%                     fprintf('  could not correct block %d to timeline %d\n', expNums(e), expNums(eTL));
%                 end
%             end
%         end
%     elseif isMpep(e)
%         for eTL = 1:length(expNums)
%             if hasTimeline(eTL)
%                 fprintf('trying to correct mpep %d to timeline %d\n', expNums(e), expNums(eTL));
%                 p = pars{e}.Protocol;
%                 nStims = numel(p.seqnums);
%                 % An mpep stimulus has constant flips, with a gap in between
%                 % stimuli. We'd like to know how many stimuli were shown, how long
%                 % each lasted, and how much gap there was in between. But we can
%                 % only really get the number of stimuli. 
% %                 minBetween = 0.2; % sec, assume at least this much time in between stimuli
% %                 maxBetweenFlips = 2/60; 
% %                 if any(strcmp(p.parnames, 'dur'))
% %                     estimatedDur = min(p.pars(strcmp(p.parnames, 'dur'),:))/10; % sec
% %                     minDur = estimatedDur*0.75;
% %                 else
% %                     estimatedDur = [];
% %                     minDur = 0.2;
% %                 end
% %                 nStims = numel(p.seqnums);
% %                 pdT = tlFlips{eTL};
% %                 pdT = pdT(pdT>lastTimes(eTL)); 
% %                 
% %                 dpdt = diff([0; pdT]);
% %                 
% %                 possibleStarts = find(dpdt>minBetween);
% %                 possibleEnds = [possibleStarts(2:end)-1; length(pdT)]; 
% %                 durs = pdT(possibleEnds)-pdT(possibleStarts);
% %                 
% %                 % gaps will be the gap *after* a stimulus
% %                 gaps = pdT(possibleStarts(2:end))-pdT(possibleEnds(1:end-1));
%                 
%                 % dang, need a better system for mpep. The problem in
%                 % Radnitz 2017-01-13 is that the photodiode was picking up
%                 % some stuff around the sync square, which led to it being
%                 % too bright in the down state to flip down for just one
%                 % stimulus. Unfortunately that screws the whole thing. 
%                 Timeline = tl{eTL}; Fs = Timeline.hw.daqSampleRate;
%                 tt = Timeline.rawDAQTimestamps;
%                 tpd = Timeline.rawDAQData(:,strcmp({Timeline.hw.inputs.name}, 'photoDiode'));
%                 
%                 % before 13.02.2019 ---------------------------------------
% %                 sig = conv(diff([0; tpd]).^2, ones(1,16/1000*Fs), 'same');
%                 % after ---------------------------------------------------
%                 sig = medfilt1(abs(diff([0; tpd])),11);
%                 % ---------------------------------------------------------
%                 figure; 
%                 plot(tt(tt>lastTimes(eTL)), sig(tt>lastTimes(eTL)));
%                 title(sprintf('expect %d stims', nStims));
%                 
%                 mpepStart = input('start of mpep? ');
%                 mpepEnd = input('end of mpep? ');
%                 thresh = input('lower/upper thresh? ');
%                 
%                 [flipTimes, flipsUp, flipsDown] = schmittTimes(tt, sig, thresh);
%                 
%                 flipsUp = flipsUp(flipsUp>mpepStart & flipsUp<mpepEnd);
%                 flipsDown = flipsDown(flipsDown>mpepStart & flipsDown<mpepEnd);
%                 
%                 % ONLY FOR FLICKERING STIMULI
%                 hold on
%                 plot(flipsUp,ones(size(flipsUp)).*.1,'>g')
%                 plot(flipsDown,ones(size(flipsDown)).*.1,'<r')
%                 % --------------------------------------------------------
%                 
%                 if ~strcmp(p.xfile, ...
%                         'stimFlickeringCheckerboard_blackWhite.x')
%                     skippedFrames = (flipsUp(2:end)-flipsDown(1:end-1))<0.05; % assume stimuli are longer than 100ms
%                     flipsUp = flipsUp(~[false; skippedFrames]);
%                     flipsDown = flipsDown(~[skippedFrames; false]);
%                 end
%                 
%                 if nStims==length(flipsUp) || strcmp(p.xfile, ...
%                         'stimFlickeringCheckerboard_blackWhite.x')
%                     fprintf(1, 'success\n');
%                     stimOnsets = flipsUp;
%                     stimOffsets = flipsDown;
%                     
%                     writeNPY(stimOnsets, fullfile(alignDir, ...
%                         sprintf('mpep_%d_onsets_in_timeline_%d.npy', ...
%                         expNums(e), expNums(eTL))));
%                     writeNPY(stimOffsets, fullfile(alignDir, ...
%                         sprintf('mpep_%d_offsets_in_timeline_%d.npy', ...
%                         expNums(e), expNums(eTL))));
%                     
%                     lastTimes(eTL) = stimOffsets(end);
%                 else
%                     fprintf(1, 'unsuccessful - found %d, expected %d\n', length(flipsUp), nStims);
%                 end
%             end
%         end
%     end
% end

%% load the dataset

sp = loadAllKsDir(subject, date);

%% plot drift map, save spike depth and waveforms
samplingRate = sp(1).sample_rate;
bTLtoMaster = readNPY(fullfile(alignDir, ...
    sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, tags{1})));
bEphysToMaster = readNPY(fullfile(alignDir, ...
    sprintf('correct_ephys_%s_to_ephys_%s.npy', sp(2).name, ...
    sp(1).name)));
for probe = 2:length(sp)
    ksDir = getKSdir(subject, date, sp(probe).name);
    [~,~,sd] = ksDriftmap(ksDir);
    sp(probe).spikeDepths = sd;  % more accurate version of depths based on PCs.
    if all(sp(probe).cgs == 3)
        incl = find(sp(probe).spikeAmps > 40);
        str = 'not sorted';
        strSave = 'unsorted';
    else
        incl = find(ismember(sp(probe).clu,sp(probe).cids));
%         incl = find(ismember(sp(probe).clu,sp(probe).cids(sp(probe).cgs==2)));
        str = 'only SUA and MUA';
        strSave = 'sorted';
%         % save spike depths
%         save(fullfile(spikeDepthFolder, sprintf('%s_%s_%s.mat', ...
%             subject, date, sp(probe).name)), 'sd');
        
%         if probe == 2
%             params.dataDir = fullfile(subjectsFolder, subject, ...
%                 date, ['ephys_' sp(probe).name]);
%             params.fileName = sprintf('%s_%s_%s_g0_t0.imec.ap_CAR.bin', ...
%                 subject, date, sp(probe).name);
%             params.dataType = 'int16';
%             params.nCh = 385;
%             params.wfWin = round([-0.001 .002] .* samplingRate);
%             params.nWf = 2000;
% %             ind = ismember(sp(probe).clu,sp(probe).cids(sp(probe).cgs==2)) & ...
% %                 sp(probe).st < max(sp(probe).st)-1; % only include SUA
%             ind = ismember(sp(probe).clu,sp(probe).cids(sp(probe).cgs==2)); % only include SUA
%             % for spikeTimes: first convert spike times (of probe 2) to time of
%             % probe 1, then multiply by sampling rate
%             spikeTimes = round((sp(probe).st(ind) - bEphysToMaster(2)) ./ ...
%                 bEphysToMaster(1) .* samplingRate);
%             clusters = sp(probe).clu(ind);
%             clUnique = unique(clusters);
%             k = 0;
%             batchSize = 40;
%             while k < length(clUnique)
%                 ind = ismember(clusters, ...
%                     clUnique(k+1 : min(k+batchSize,length(clUnique))));
%                 params.spikeTimes = spikeTimes(ind);
%                 params.spikeClusters = clusters(ind);
%             
%                 w = getWaveForms(params);
%                 w = rmfield(w, 'waveForms');
%                 if k == 0
%                     wf = w;
%                     fields = fieldnames(wf);
%                 else
%                     for f = 1:length(fields)
%                         wf.(fields{f}) = [wf.(fields{f}); w.(fields{f})];
%                     end
%                 end
%                 k = k + batchSize;
%             end
%             save(fullfile(waveformFolder, sprintf('%s_%s_%s.mat', ...
%                 subject, date, sp(probe).name)), 'wf');
%             clear wf
%         end
    end
    skip = round(length(incl)/1000000);
    incl = incl(1:skip:end);
    figure('Position', [1 33 1175 1083])
    plotDriftmap(sp(probe).st(incl), sp(probe).spikeAmps(incl), sd(incl));
    hold on
    yLims = [min(sd) max(sd)];
    yLims(2) = yLims(2) + 0.02*range(yLims);
    for exp = TLexp+1:length(isMpep)
        if isMpep(exp)
            t1 = readNPY(fullfile(alignDir, sprintf('mpep_%d_onsets_in_timeLine_%d.npy', ...
                exp, TLexp)));
            t1 = applyCorrection(t1, bTLtoMaster);
            tend = readNPY(fullfile(alignDir, sprintf('mpep_%d_offsets_in_timeLine_%d.npy', ...
                exp, TLexp)));
            tend = applyCorrection(tend, bTLtoMaster);
            plot([1 1].*t1(1), yLims, 'Color', [0 0.5 0])
            plot([1 1].*tend(end), yLims, 'Color', [0.5 0 0])
            text(t1(1), double(yLims(end)-0.01*range(yLims)), [' ' ...
                pars{exp}.Protocol.xfile], 'Color', [0 0.5 0], 'Interpreter', 'none')
        elseif hasBlock(exp)
            file = fullfile(alignDir, sprintf('block_%d_sw_in_timeline_%d.npy', exp, TLexp));
            if ~isfile(file)
                continue
            end
            sw = readNPY(file);
            limits = applyCorrection(sw([1 end]), bTLtoMaster);
            plot([1 1].*limits(1), yLims, 'Color', [0 0.5 0])
            plot([1 1].*limits(2), yLims, 'Color', [0.5 0 0])
            if isfield(pars{exp}, 'defFunction')
                [~,prot] = fileparts(pars{exp}.defFunction);
            else
                prot = pars{exp}.type;
            end
            text(limits(1), double(yLims(end)-0.01*range(yLims)), [' ' ...
                prot], 'Color', [0 0.5 0], 'Interpreter', 'none')
        end
    end
    axis tight
    set(gca, 'box', 'off')
    title(sprintf('%s %s probe %s (%s)', subject, date, sp(probe).name, str))
    tmp = strfind(ksDir, filesep);
    folder = ksDir(1:tmp(end)-1);
    savefig(gcf, fullfile(folder, ...
        sprintf('%s_%s_%s_%s.fig', subject, date, sp(probe).name, strSave)), ...
        'compact')
    saveas(gcf, fullfile(folder, ...
        sprintf('%s_%s_%s_%s.png', subject, date, sp(probe).name, strSave)))
%     savefig(gcf, fullfile(driftPlotFolder, ...
%         sprintf('%s_%s_%s_%s.fig', subject, date, sp(probe).name, strSave)), ...
%         'compact')
end

% probe = 2;
% incl = find(sp(probe).spikeAmps > 40);
% skip = round(length(incl)/1000000);
% incl = incl(1:skip:end);
% figure('Position', [1 33 1175 1083])
% plotDriftmap(sp(probe).st(incl), sp(probe).spikeAmps(incl), sp(probe).spikeDepths(incl));
% axis tight
% set(gca, 'box', 'off')
% title(sprintf('%s %s probe %s', subject, date, sp(2).name))
% savefig(gcf, fullfile(driftPlotFolder, ...
%     sprintf('%s_%s_%s.fig', subject, date, sp(2).name)), 'compact')

%% align eye video
% alignVideo(subject, date, TLexp, 'eye')
% % resulting times are stored in Zserver\Data\EyeCamera; times need to be
% % aligned to ephys times: eyeTimes.*timelineToEphys(1)+timelineToEphys(2)

%% Plot photodiode in ephys time to determine times of darkness
% bTLtoMaster = readNPY(fullfile(alignDir, ...
%     sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, tags{1})));
% photodiode = tl{TLexp}.rawDAQData(:,strcmp({tl{TLexp}.hw.inputs.name}, 'photoDiode'));
% tlTime = tl{TLexp}.rawDAQTimestamps;
% t = applyCorrection(tlTime, bTLtoMaster);
% figure('Position', [2 678 1917 420])
% plot(t, photodiode, 'k')
% xlabel('Time (aligned to ephys)')
% ylabel('photodiode')
% xlim([0 t(end)])
