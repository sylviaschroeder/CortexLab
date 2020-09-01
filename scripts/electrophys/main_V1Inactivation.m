%% Dataset
db_ephys_V1

%% Define folders
protocolFolder = '\\ZUBJECTS.cortexlab.net\Subjects';
subjectsFolder = 'J:\Ephys';
folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\Electrophys_SC\V1_inactivation\';
plotFolder = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\Electrophys_SC\V1_inactivation';

%% Parameters
binSizeGrating = 0.005;
runFitTwice = 1;
doSave = 1;

%% Constants
tags.stimDur = {'dur'};
tags.stimOnset = {'ton'};
tags.stimOffset = {'toff'};
tags.laserOn = {'tstart1'};
tags.laserOff = {'tend1'};
tags.ori = {'ori1','ori'};
tags.contrast = {'c1', 'cg'};
tags.amplitude = {'amp1'};

%% Get data for each neuron (no fitting!): tuning curve with vs. without laser
binSizeGrating = 0.01;
prePostStim = 0.2;

responses = cell(1, length(db)); % each entry: [neuron x stim x rep x time]
stimuli = cell(1, length(db)); % each entry: 1st row: direction, 2nd row: 
% blank, 3rd row: laser on
cellIDs = cell(1, length(db));
depths = cell(1, length(db));
cellClass = cell(1, length(db)); % 1: MUA, 2: SUA
binsGrating = cell(1, length(db));
stimBins = cell(1, length(db));

for iExp = 1:length(db)
    fprintf('\nProcessing: %s %s\n', db(iExp).subject, db(iExp).date);
    alignDir = fullfile(subjectsFolder, db(iExp).subject, db(iExp).date, 'alignments');
    
    % load spike data
    probe = db(iExp).probe;
    tag = '';
    if ~isempty(db(iExp).probeNames)
        tag = ['_' db(iExp).probeNames{probe}];
    end
    sp = loadAllKsDir(db(iExp).subject, db(iExp).date);
    [expNums, ~, ~, ~, ~, ~, hasTimeline] = ...
        dat.whichExpNums(db(iExp).subject, db(iExp).date);
    TLexp = expNums(hasTimeline);
    TLexp = TLexp(end);
    bTLtoMaster = readNPY(fullfile(alignDir, ...
        sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, db(iExp).probeNames{1})));
    rawF = fullfile(subjectsFolder, db(iExp).subject, ...
        db(iExp).date, ['ephys' tag], 'sorting');
    if isfolder(rawF)
        [~, ~, spikeDepths] = ksDriftmap(rawF);
    end
    
    % load stimulus info
    data = load(fullfile(protocolFolder, db(iExp).subject, db(iExp).date, ...
        num2str(db(iExp).expOri), sprintf('%s_%d_%s_parameters.mat', ...
        db(iExp).date, db(iExp).expOri, db(iExp).subject)));
    parsGratings = data.parameters.Protocol;
    stimOnTL = readNPY(fullfile(alignDir, ...
        sprintf('mpep_%d_onsets_in_timeline_%d.npy', db(iExp).expOri, TLexp)));
    stimOffTL = readNPY(fullfile(alignDir, ...
        sprintf('mpep_%d_offsets_in_timeline_%d.npy', db(iExp).expOri, TLexp)));
    stimOn = applyCorrection(stimOnTL, bTLtoMaster);
    stimOff = applyCorrection(stimOffTL, bTLtoMaster);
    stimDur = mean(stimOff - stimOn);
    window = [-prePostStim stimDur+prePostStim];
    binBorders = window(1) : binSizeGrating : window(2);
    binsGrating{iExp} = binBorders(1:end-1) + binSizeGrating/2;
    stimBins{iExp} = binsGrating{iExp}>0 & binsGrating{iExp}<stimDur;
    stimIDs = repmat((1:parsGratings.npfilestimuli)',parsGratings.nrepeats,1);
    [~,order] = sort(parsGratings.seqnums(:));
    stimSeq = stimIDs(order);
    blank = parsGratings.pars(strcmp(parsGratings.parnames,'c1'),:) == 0;
    directions = parsGratings.pars(strcmp(parsGratings.parnames,'ori1'),:);
    directions(blank) = NaN;
    laserOn = parsGratings.pars(strcmp(parsGratings.parnames,'amp1'),:) > 0;
    stimuli{iExp} = [directions; blank; laserOn];
    repeats = setdiff(1:parsGratings.nrepeats, db(iExp).excludeReps);
    
    % get response for each cell
    units = sp(probe).cids;
    included = true(1, length(units));
    resp = NaN(length(units), length(directions), length(repeats), ...
        length(binsGrating{iExp}));
    depths{iExp} = NaN(length(units),1);
    for iCell = 1:length(units)
        spInd = sp(probe).clu == units(iCell);
        depths{iExp}(iCell) = nanmean(spikeDepths(spInd));
        if depths{iExp}(iCell)<db(iExp).V1depth(1) || depths{iExp}(iCell)>db(iExp).V1depth(2)
            included(iCell) = false;
            continue
        end
        for stim = 1:length(directions)
            stimInd = find(stimSeq == stim);
            stimInd = stimInd(repeats);
            binnedArray = timestampsToBinned(sp(probe).st(spInd), ...
                stimOn(stimInd), binSizeGrating, window);
            resp(iCell,stim,:,:) = binnedArray ./ binSizeGrating;
        end
    end
    responses{iExp} = resp(included,:,:,:);
    cellIDs{iExp} = units(included);
    depths{iExp} = depths{iExp}(included);
    cellClass{iExp} = sp(probe).cgs(included);
end
save(fullfile(folderResults, 'gratingResponses.mat'), 'responses', ...
    'stimuli', 'cellIDs', 'depths', 'cellClass', 'binsGrating', 'stimBins')

%% Plot inactivation index
layerBorders = [.12, .33, .47, .73];
inactivationInds = cell(1,length(responses));
depthsNorm = cell(1,length(responses));
maxResponses = cell(1,length(responses));
for iExp = 1:length(responses)
    % find pref. stim. and resp. to pref. stim. for each neuron (only when
    % laser is off)
    offStims = find(stimuli{iExp}(3,:) == 0 & stimuli{iExp}(2,:) == 0);
    binOnOff = [find(stimBins{iExp},1,'first'), find(stimBins{iExp},1,'last')];
    stimOff = binsGrating{iExp}(binOnOff(2));
    binOnOff(2) = find(binsGrating{iExp} > stimOff-0.1, 1);
    meanResps = mean(mean(responses{iExp}(:,:,:,binOnOff(1):binOnOff(2)), 4), 3); % [neuron x stim]
    [maxResps, prefStims] = max(meanResps(:,offStims),[],2);
    prefStims = offStims(prefStims)';
    onStims = NaN(1,max(offStims));
    stims = stimuli{iExp};
    stims(isnan(stims)) = 0; % set orientation of blanks to 0
    for st = 1:length(offStims)
        onStims(offStims(st)) = find(all(stims == [stims(1:2,offStims(st));1],1),1);
    end
    ind = sub2ind(size(meanResps), 1:size(meanResps,1), onStims(prefStims));
    maxResps_laser = meanResps(ind)';
    inactivationInds{iExp} = maxResps_laser ./ maxResps;
    maxResponses{iExp} = maxResps;
    
    % normalise depth
    depthsNorm{iExp} = (db(iExp).V1depth(2)-depths{iExp}) ./ diff(db(iExp).V1depth);
end

inact = cat(1,inactivationInds{:});
inact(inact>2) = 2;

mr = cat(1, maxResponses{:});
d = cat(1, depthsNorm{:});

ind = mr > 2;
figure
plot(inact(ind), d(ind), 'o')
hold on
plot([0 2],repmat(layerBorders,2,1),'k')
set(gca,'XTick',[0 1 2])
set(gca, 'YDir', 'reverse')
xlabel('Resp. @ preferred stimulus: laser on / laser off')
ylabel('Depth (normalised)')
title(sprintf('Inactivation in V1 (n = %d)', sum(ind)))

% Plot PSTH for each neuron
for iExp = 1:length(responses)
    for iCell = 1:size(responses{iExp},1)
        if maxResponses{iExp}(iCell) < 2
            continue
        end
        psth = squeeze(mean(responses{iExp}(iCell,:,:,:),3));
        figure
        imagesc(binsGrating{iExp}([1 end]), [1 size(stimuli{1},2)], psth)
        xlabel('Time from stimulus onset (s)')
        ylabel('Stimulus')
        title(sprintf('%s %s cell %d: %.2f laser resp, %.2f depth', ...
            db(iExp).subject, db(iExp).date, iCell, ...
            inactivationInds{iExp}(iCell), depthsNorm{iExp}(iCell)))
    end
end

% 2nd version of PSTHs (for paper)
k = 2;
cells = [2 34 63 28 55 18 75 73];
cm = colormap('gray');
cm = flip(cm);
for j = 1:length(cells)
    ind = cells(j);
    psth = squeeze(mean(responses{k}(ind,:,:,:),3));
    maxi = max(psth(:));
    figure
    subplot(2,1,1)
    imagesc(binsGrating{iExp}([1 end]), [1 size(stimuli{k},2)/2], ...
        psth(1:size(psth,1)/2,:),[0 maxi])
    ylim([0.5 12.5])
    title(sprintf('%s %s cell %d: %.2f laser resp, %.2f depth', ...
        db(k).subject, db(k).date, ind, inactivationInds{k}(ind), ...
        depthsNorm{k}(ind)))
    ylabel('laser on')
    colorbar
    subplot(2,1,2)
    imagesc(binsGrating{iExp}([1 end]), [1 size(stimuli{k},2)/2], ...
        psth(size(psth,1)/2+1:end,:),[0 maxi])
    ylim([0.5 12.5])
    xlabel('Time from stimulus onset')
    ylabel('laser off')
    colormap(cm)
    colorbar
end

%% Plot tuning curves (NOT fitted) for each neuron
conds = [false, true];
colors = 'kc';
labels = {'no laser','with laser'};
xOffsets = [-2 2];
classes = {'MUA','SUA'};
for iExp = 1:length(db)
    folder = fullfile(plotFolder, sprintf('%s_%s', db(iExp).subject, db(iExp).date));
    if ~isfolder(folder)
        mkdir(folder)
    end
    for iCell = 1:size(responses{iExp},1)
        resp = nanmean(squeeze(responses{iExp}(iCell,:,:,stimBins{iExp})),3); % [stim x rep]
        m = nanmean(resp,2);
        sd = nanstd(resp,0,2)./sqrt(size(resp,2));
        figure
        hold on
        h = [0 0];
        for c = 1:2
            gratings = find(stimuli{iExp}(3,:)==conds(c) & stimuli{iExp}(2,:)==0);
            blank = stimuli{iExp}(3,:)==conds(c) & stimuli{iExp}(2,:)==1;
            dirs = stimuli{iExp}(1,gratings);
            dirs = [dirs, dirs(1)+360];
            h(c) = errorbar(dirs+xOffsets(c), m(gratings([1:end 1])), ...
                sd(gratings([1:end 1])), [colors(c) 'o-'], 'LineWidth', 2);
            fill(dirs([1 end end 1]), [[1 1].*(m(blank)+sd(blank)), ...
                [1 1].*(m(blank)-sd(blank))], 'k', 'EdgeColor', 'none', ...
                'FaceColor', colors(c), 'FaceAlpha', 0.3)
            plot(dirs([1 end]), [1 1].*m(blank), [colors(c) '--'], 'LineWidth', 2)
        end
        set(gca,'box','off','XTick',0:90:360)
        xlim([-5 365])
        xlabel('Direction')
        ylabel('Firing rate (sp/s)')
        legend(h, labels)
        title(sprintf('Neuron %d (%s), depth: %.0f um', cellIDs{iExp}(iCell), ...
            classes{cellClass{iExp}(iCell)}, depths{iExp}(iCell)))
        
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(fullfile(folder, sprintf('%04d_neuron%04d.jpg', ...
            round(depths{iExp}(iCell)), cellIDs{iExp}(iCell))), '-djpeg','-r0')
        close(fig)
    end
end

%% Compute tuning curves and linear fits for running and not running
% parameters: (1) pref. dir., (2) ampl. at pref. dir, (3) direction index,
%             (4) offset of tuning curve, (5) tuning width
parameterSets = {[1 3 5], 'constant'};
fileNames = {'tuning_prefDirSigmaDIFixed.mat', 'tuning_constantFit.mat'};
classes = {'MUA','SUA'};

if ~exist(folderResults, 'dir')
    mkdir(folderResults);
end
x = 0:359;
for iPars = 1:length(parameterSets)
    fprintf('Parameter set %d\n', iPars)
    fixedPars = parameterSets{iPars};
    tuning = struct([]);
    
    for iExp = 1:length(db)
        fprintf('  Dataset %d: %s %s exp.: %d\n', iExp, db(iExp).subject, ...
            db(iExp).date, db(iExp).expOri);
        tuning(iExp).subject = db(iExp).subject;
        tuning(iExp).date = db(iExp).date;
        tuning(iExp).exp = db(iExp).expOri;
        alignDir = fullfile(subjectsFolder, db(iExp).subject, db(iExp).date, 'alignments');
        
        % load spike data
        probe = db(iExp).probe;
        sp = loadAllKsDir(db(iExp).subject, db(iExp).date);
        [expNums, ~, ~, ~, ~, ~, hasTimeline] = ...
            dat.whichExpNums(db(iExp).subject, db(iExp).date);
        TLexp = expNums(hasTimeline);
        TLexp = TLexp(end);
        bTLtoMaster = readNPY(fullfile(alignDir, ...
            sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, db(iExp).probeName)));
        rawF = fullfile(subjectsFolder, db(iExp).subject, ...
            db(iExp).date, 'ephys', 'sorting');
        if isfolder(rawF)
            [~, ~, spikeDepths] = ksDriftmap(rawF);
        end
        
        % load stimulus info
        data = load(fullfile(protocolFolder, db(iExp).subject, db(iExp).date, ...
            num2str(db(iExp).expOri), sprintf('%s_%d_%s_parameters.mat', ...
            db(iExp).date, db(iExp).expOri, db(iExp).subject)));
        parsGratings = data.parameters.Protocol;
        stimOnTL = readNPY(fullfile(alignDir, ...
            sprintf('mpep_%d_onsets_in_timeline_%d.npy', db(iExp).expOri, TLexp)));
        stimOffTL = readNPY(fullfile(alignDir, ...
            sprintf('mpep_%d_offsets_in_timeline_%d.npy', db(iExp).expOri, TLexp)));
        stimOn = applyCorrection(stimOnTL, bTLtoMaster);
        stimOff = applyCorrection(stimOffTL, bTLtoMaster);
        
        stimIDs = repmat((1:parsGratings.npfilestimuli)',parsGratings.nrepeats,1);
        [~,order] = sort(parsGratings.seqnums(:));
        stimSeq = stimIDs(order);
        stimDur = round(parsGratings.pars(strcmp(parsGratings.parnames, ...
            tags.stimDur{1}),:)' ./ 10 .* 60) ./ 60; % 60 Hz monitor frame rate
        stimOn_rel = round(parsGratings.pars(strcmp(parsGratings.parnames, ...
            tags.stimOnset{1}),:)' ./ 1000 .* 60) ./ 60;
        stimOff_rel = round(parsGratings.pars(strcmp(parsGratings.parnames, ...
            tags.stimOffset{1}),:)' ./ 1000 .* 60) ./ 60 - stimDur;
        stimOn_rel = stimOn_rel(stimSeq);
        stimOff_rel = stimOff_rel(stimSeq);
        stimOn = stimOn + stimOn_rel;
        stimOff = stimOff + stimOff_rel;
        
        stimDur = mean(stimOff - stimOn);
        tag = tags.contrast{ismember(tags.contrast, parsGratings.parnames)};
        blank = parsGratings.pars(strcmp(parsGratings.parnames,tag),:) == 0;
        tag = tags.ori{ismember(tags.ori, parsGratings.parnames)};
        directions = parsGratings.pars(strcmp(parsGratings.parnames,tag),:);
        directions(blank) = NaN;
        laserOn = parsGratings.pars(strcmp(parsGratings.parnames,tags.amplitude),:) > 0;
        repeats = setdiff(1:parsGratings.nrepeats, db(iExp).excludeReps);
        
        % get response for each cell
        units = sp(probe).cids;
        depths = NaN(length(units),1);
        for iCell = 1:length(units)
            depths(iCell) = nanmean(spikeDepths(sp(probe).clu == units(iCell)));
        end
        included = depths > db(iExp).V1depth(1) & depths < db(iExp).V1depth(2);
        units = units(included);
        tuning(iExp).cellIDs = units;
        tuning(iExp).depths = depths(included);
        tuning(iExp).class = classes(sp(probe).cgs(ismember(sp(probe).cids,units)));
        tuning(iExp).R2 = struct([]);
        tuning(iExp).isSuppressed = NaN(length(units),1);
        tuning(iExp).crossValExplVar = NaN(length(units),1);
        tuning(iExp).laserOn = laserOn;
        tuning(iExp).directions = directions;
        
        fprintf('    Cell (of %d):', length(units));
        conditions = repmat(laserOn'+1,1,length(repeats));
        for iCell = 1:length(units)
            fprintf(' %d', iCell)
            resp = NaN(length(directions), length(repeats));
            spInd = sp(probe).clu == units(iCell);
            for stim = 1:length(directions)
                stimInd = find(stimSeq == stim);
                stimInd = stimInd(repeats);
                binnedArray = timestampsToBinned(sp(probe).st(spInd), ...
                    stimOn(stimInd), binSizeGrating, [0 stimDur]);
                resp(stim,:) = sum(binnedArray,2) ./ stimDur;
            end
            tuning(iExp).cell(iCell).response = resp;
            baselineLaserOff = mean(reshape(resp(~laserOn & blank,:),[],1));
            respLaserOff = mean(resp(~laserOn & ~blank,:),2) - ...
                baselineLaserOff;
            [~,maxInd] = max(abs(respLaserOff));
            respSign = sign(respLaserOff(maxInd));
            tuning(iExp).isSuppressed(iCell) = -respSign;
            % if neuron is suppressed, multiply stim responses with -1
            if respSign < 0
                respFit = -resp;
            else
                respFit = resp;
            end
            
            curves = NaN(length(x),2);
            if strcmp(fixedPars, 'constant')
                % crossvalidate
                errors = models.crossvalidate( @models.fitConstant, ...
                    {resp'}, {directions, blank, conditions'}, 'stimSets');
                
                parameters = NaN(1, 2);
                predictions = NaN(size(resp));
                for c = 1:2
                    ind = conditions == c & ~blank';
                    parameters(c) = nanmean(resp(ind));
                    predictions(ind) = parameters(c);
                    curves(:,c) = parameters(c);
                end
            else
                % crossvalidate
                errors = models.crossvalidate( ...
                    @gratings.fitTuningCurveConditions_forCrossVal, ...
                    {respFit'}, {[directions(~blank)', ...
                    setdiff(1:length(directions),find(blank))'], ...
                    find(blank), conditions', fixedPars}, 'stimSets');
                
                [parameters, ~, predictions, R2] = ...
                    gratings.fitTuningCurveConditions(respFit, ...
                    [directions(~blank)', setdiff(1:length(directions),find(blank))'], ...
                    find(blank), conditions, fixedPars);
                if respSign < 0
                    parameters(2,:) = -parameters(2,:);
                    parameters(4,:) = -parameters(4,:);
                    predictions = -predictions;
                end
                for c = 1:2
                    curves(:,c) = gratings.orituneWrappedConditions(parameters(:,c),x);
                end
            end
            
            err = errors{1}(:,~blank)';
            respNoBlanks = resp(~blank,:);
            ind = ~isnan(err) & ~isnan(respNoBlanks);
            explVar = 1 - sum(err(ind).^2) / ...
                sum((respNoBlanks(ind)-mean(respNoBlanks(ind))).^2);
            tuning(iExp).crossValExplVar(iCell) = explVar;
            
            predNoBlanks = predictions(~blank,:);
            ind = ~isnan(respNoBlanks) & ~isnan(predNoBlanks);
            R2 = 1 - sum((respNoBlanks(ind) - predNoBlanks(ind)).^2) / ...
                sum((respNoBlanks(ind) - mean(respNoBlanks(ind))).^2);
            R2ToBlank = 1 - sum((respNoBlanks(ind) - predNoBlanks(ind)).^2) / ...
                sum((respNoBlanks(ind) - ...
                nanmean(reshape(resp(blank,:),[],1))).^2);
            tuning(iExp).R2(iCell).comparedToMean = R2;
            tuning(iExp).R2(iCell).comparedToBlankResp = R2ToBlank;
            for c = 1:size(parameters,2)
                tuning(iExp).cond(c).cell(iCell).parameters = parameters(:,c);
                tuning(iExp).cond(c).cell(iCell).curve = curves(:,c);
                tuning(iExp).cond(c).cell(iCell).blankResponses = ...
                    resp(conditions==c & repmat(blank',1,length(repeats)));
            end
        
            if doSave == 1
                save(fullfile(folderResults, fileNames{iPars}), 'tuning', 'x')
            end
        end
    end
end

%% Add to tuning structure which neurons are not tuned
data = load(fullfile(folderResults, 'tuning_prefDirSigmaDIFixed.mat'));
tuning = data.tuning;
x = data.x;
data = load(fullfile(folderResults, 'tuning_constantFit.mat'));
tunC = data.tuning;

for iExp = 1:length(tuning)
    tuned = zeros(length(tuning(iExp).cellIDs), 1);
    tuned(tuning(iExp).crossValExplVar > tunC(iExp).crossValExplVar) = 1;
    tuning(iExp).isTuned = tuned;
    ind = find(tuned == 0);
    for k = 1:length(ind)
        for c = 1:2
            par = tunC(iExp).cond(c).cell(ind(k)).parameters;
            tuning(iExp).cond(c).cell(ind(k)).parameters = par;
            tuning(iExp).cond(c).cell(ind(k)).curve = ...
                ones(length(x),1) .* par;
        end
        tuning(iExp).crossValExplVar(ind(k)) = ...
            tunC(iExp).crossValExplVar(ind(k));
        tuning(iExp).R2(ind(k)) = ...
            tunC(iExp).R2(ind(k));
    end
end
save(fullfile(folderResults, 'tuning_prefDirSigmaDIFixed_isTuned.mat'), ...
    'tuning', 'x');

%% Add line fits
data = load(fullfile(folderResults, 'tuning_prefDirSigmaDIFixed_isTuned.mat'));
tuning = data.tuning;
x = data.x;

for iExp = 1:length(tuning)
    if isempty(tuning(iExp).cellIDs)
        continue
    end
    for iCell = 1:length(tuning(iExp).cellIDs)
        dat1 = tuning(iExp).cond(1).cell(iCell); % x-axis: laser off
        if isempty(dat1.parameters)
            continue
        end
        dat2 = tuning(iExp).cond(2).cell(iCell); % y-axis: laser on
        if tuning(iExp).isTuned(iCell) == 1
            lineFit.intercept = dat2.parameters(4) - dat1.parameters(4) * ...
                dat2.parameters(2) / dat1.parameters(2);
            lineFit.slope = dat2.parameters(2) / dat1.parameters(2);
        else
            lineFit.intercept = dat2.parameters - dat1.parameters;
            lineFit.slope = NaN;
        end
        tuning(iExp).lineFit(iCell) = lineFit;
    end
end
save(fullfile(folderResults, 'tuning_prefDirSigmaDIFixed_lineFit.mat'), ...
    'tuning', 'x');

%% Plot slopes and intercepts of all neurons
data = load(fullfile(folderResults, 'tuning_prefDirSigmaDIFixed_lineFit.mat'));
tuning = data.tuning;

slopes = [];
intercepts = [];
mean_laserOn = [];
mean_laserOff = [];
depths = [];
for iExp = 1:length(tuning)
    slopes = [slopes; [tuning(iExp).lineFit.slope]'];
    intercepts = [intercepts; [tuning(iExp).lineFit.intercept]'];
    r = cat(3, tuning(iExp).cell.response); % [stim x rep x neuron]
    mean_laserOn = [mean_laserOn; ...
        squeeze(mean(mean(r(tuning(iExp).laserOn,:,:),2),1))];
    mean_laserOff = [mean_laserOff; ...
        squeeze(mean(mean(r(~tuning(iExp).laserOn,:,:),2),1))];
    V1 = db(iExp).V1depth;
    depths = [depths; (V1(2)-tuning(iExp).depths)/diff(V1)];
end
slopes(slopes<-5) = -5;
slopes(slopes>2) = 2;
intercepts(intercepts<-5) = -5;
intercepts(intercepts>5) = 5;
mean_ratios = mean_laserOn ./ mean_laserOff;
mean_ratios(mean_ratios>2) = 2;

tuned = ~isnan(slopes);
figure
plot(intercepts(tuned), slopes(tuned), 'k.')
xlabel('Intercept')
ylabel('Slope')
title('Tuned neurons')
figure
hist(intercepts(~tuned),-5:5)
xlabel('Intercept')
title('Untuned neurons')

figure('Position',[1570 550 300 550])
scatter(mean_ratios, depths)
hold on
plot([1 1],[0 1],'k')
xlim([-.2 2.2])
set(gca,'XTick',0:2,'XTickLabel',{'0','1','>=2'},'YDir','reverse')
xlabel('Response: laser on / laser off')
ylabel('Depth relative to V1 surface (scaled to max 1)')
