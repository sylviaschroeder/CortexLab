folderData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\electrophys\';
folderEye = '\\zserver.cortexlab.net\Data\EyeCamera';
% folderEye = 'J:\Eye';
folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\Electrophys_SC\nonVisualEffects\';

%% Define datasets
db_ephys_driftingGratings

%% Parameters
nonVisualSignal = 'pupil'; %'running' or 'pupil';

if strcmp(nonVisualSignal,'running')
    threshold = 100; %5; % 100 corresponds to 13 cm/s (speed = wheel * cmPerUnit; 
                         % cmPerUnit = 2*pi * radius / (4 * countsPerRotation)
    minDurPerTrial = .4;
    labels = {'no running','running'};
else
    threshPerc = 50;
    minDurPerTrial = .5;
    labels = {'small pupil, laser off','large pupil, laser off', ...
        'small pupil, laser on', 'large pupil, laser on'};
end

doPlot = 0;
doSave = 1;
runFitTwice = 1;

%% Compute tuning curves and linear fits for running and not running
% parameters: (1) pref. dir., (2) ampl. at pref. dir, (3) direction index,
%             (4) offset of tuning curve, (5) tuning width
parameterSets = {[], [1 3 5], 'constant', [1], [1 5], [1 2 3 5], [1 3 4 5]};
fileNames = {'tuning_nothingFixed.mat', ...
    'tuning_prefDirSigmaDIFixed.mat', ...
    'tuning_constantFit.mat', ...
    'tuning_prefDirFixed.mat', ...
    'tuning_prefDirSigmaFixed.mat', ...
    'tuning_prefDirAmpSigmaDIFixed.mat', ...
    'tuning_prefDirSigmaDIBaselineFixed.mat'};
if ~exist(fullfile(folderResults, nonVisualSignal), 'dir')
    mkdir(fullfile(folderResults, nonVisualSignal));
end
x = 0:359;
for iPars = 1:length(parameterSets)
    fprintf('Parameter set %d\n', iPars)
    fixedPars = parameterSets{iPars};
    tuning = struct([]);
    for iExp = 1:length(db)
        fprintf('  Dataset %d: %s %s exp.: %d\n', iExp, db(iExp).subject, ...
            db(iExp).date, db(iExp).exp);
        tuning(iExp).subject = db(iExp).subject;
        tuning(iExp).date = db(iExp).date;
        tuning(iExp).exp = db(iExp).exp;
        % load data (time, spikeCounts, cellIDs, dephts, stimMatrix,
        % directions, blanks, stimTimes, laserOn, laserTimesRelative,
        % nSpikes, sampleRate, spikeWidths, stimSequence, timelineToEphys,
        % waveforms)
        load(fullfile(folderData, db(iExp).subject, db(iExp).date, ...
            sprintf('%02d_data.mat', db(iExp).exp)));
        notBlanks = true(size(stimMatrix,1),1);
        notBlanks(blanks) = false;
        tuning(iExp).cellIDs = cellIDs;
        tuning(iExp).R2 = struct([]);
        tuning(iExp).isSuppressed = NaN(length(cellIDs),1);
        tuning(iExp).crossValExplVar = NaN(length(cellIDs),1);
        
        nonVisData = [];
        % load ball or pupil data
        switch nonVisualSignal
            case 'running'
                % load running
                load(fullfile(folderData, db(iExp).subject, db(iExp).date, ...
                    sprintf('%02d_running.mat', db(iExp).exp)));
                nonVisData = running;
                nonVisTime = time;
                ylab = 'Running speed (cm/s)';
            case 'pupil'
                if ~exist(fullfile(folderEye, db(iExp).subject, db(iExp).date, ...
                        sprintf('%02d_eye_processed.mat', db(iExp).exp)), 'file')
                    continue
                end
                data = load(fullfile(folderEye, db(iExp).subject, db(iExp).date, ...
                    sprintf('%02d_eye_processed.mat', db(iExp).exp)), 'results');
                pupil = data.results;
                nonVisData = nonVis.getPupilDiam(pupil);
                data = load(fullfile(folderEye, db(iExp).subject, db(iExp).date, ...
                    sprintf('%02d_eyeTime.mat', db(iExp).exp)));
                nonVisTime = data.eyeTime;
                threshold = prctile(nonVisData, threshPerc);
                ylab = 'Pupil diameter (a.u.)';
        end
        if doPlot == 1
            figure('Position',[1925 685 1915 420])
            plot(nonVisTime, nonVisData)
            hold on
            plot(nonVisTime([1 end]), [1 1] * threshold)
            xlabel('Time (s)')
            ylabel(ylab)
            title(sprintf('Dataset %d: %s %s exp.: %d\n', iExp, ...
                db(iExp).subject, db(iExp).date, ...
                db(iExp).exp),'Interpreter','none')
            xlim(nonVisTime([1 end]))
        end
        
        nanInds = isnan(nonVisData);
        nonVisual = interp1(nonVisTime(~nanInds), nonVisData(~nanInds), ...
            time, 'pchip');
        nonVisual = double(nonVisual > threshold);
        nanTimes = nonVisTime(isnan(nonVisData));
        nanTimes = hist(nanTimes, time)>0;
        nonVisual(nanTimes) = NaN;
        nonVisual = squeeze(ssLocal.getTracesPerStimulus(nonVisual, ...
            stimMatrix, [0 0])); % [stimuli x repetitions x time]
        % discard trials where nonvisual signal is NaN more than half of
        % the times
        ind = sum(isnan(nonVisual),3) > size(nonVisual,3)/2;
        nonVisual = nanmean(nonVisual,3); % [stimuli x repetition]
        nonVisual = double(nonVisual >= minDurPerTrial) + 1;
        nonVisual(ind) = NaN;
        % conditions: 1: small nonVis + no laser, 2: large nonVis + no laser,
        % 3: small nonVis + laser on, 4: large nonVis + laser on
        %     conditions = bsxfun(@plus, nonVisual, laserOn * 2);
        conditions = nonVisual;
%         tuning(iExp).nonVisual = nonVisual;
%         tuning(iExp).laserOn = laserOn;
%         tuning(iExp).directions = directions;
        for c = 1:4
            tuning(iExp).cond(c).name = labels{c};
        end
        
        responses = ssLocal.getTracesPerStimulus(spikeCounts, ...
            stimMatrix, [0 0]); % [neuron x stim x rep x time]
        stimDur = mean(stimTimes.offset - stimTimes.onset);
        responses = sum(responses,4) ./ stimDur; % in Hz, [neuron x stim x rep]
        
        fprintf('    Cell (of %d):', length(cellIDs));
        for iCell = 1:length(cellIDs)
            fprintf(' %d', iCell)
            resp = squeeze(responses(iCell,:,:)); % [stim x rep]
%             tuning(iExp).cell(iCell).response = resp;
            baselineLaserOff = mean(reshape(resp(~laserOn & ~notBlanks,:),[],1));
            respLaserOff = mean(resp(~laserOn & notBlanks,:),2) - ...
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
            
            curves = NaN(length(x),4);
            conds = conditions + laserOn .* 2;
            if strcmp(fixedPars, 'constant')
                % crossvalidate
                errors = models.crossvalidate( @models.fitConstant, ...
                    {resp'}, {directions, blanks, conds'}, 'stimSets');
                
                parameters = NaN(1, 4);
                predictions = NaN(size(resp));
                blankResponses = cell(1,4);
                for c = 1:4
                    ind = conds == c & notBlanks;
                    parameters(c) = nanmean(resp(ind));
                    predictions(ind) = parameters(c);
                    curves(:,c) = parameters(c);
                    blankResponses{c} = resp(conds==c & ~notBlanks);
                end
            else
                % crossvalidate
                errors = models.crossvalidate( ...
                    @gratings.fitTuningCurveConditionsAndLaser_forCrossVal, ...
                    {respFit'}, {directions, conditions', laserOn', ...
                    fixedPars}, 'stimSets');
                
                [parsLaserOff, parsLaserOn, predLaserOff, predLaserOn, ...
                    blanksLaserOff, blanksLaserOn] = ...
                    gratings.fitTuningCurveConditionsAndLaser(respFit, directions, ...
                    conditions, laserOn, fixedPars, runFitTwice);
                
                curves(:,1:2) = gratings.orituneWrappedConditions(parsLaserOff, ...
                    [x' x'], repmat([1 2], length(x), 1));
                curves(:,3:4) = gratings.orituneWrappedConditions(parsLaserOn, ...
                    [x' x'], repmat([1 2], length(x), 1));
                parameters = reshape([parsLaserOff, parsLaserOn]',5,[]);
                parameters(:,2) = sum(parameters(:,1:2),2);
                parameters(:,4) = sum(parameters(:,3:4),2);
                predictions = NaN(size(resp));
                predictions(~laserOn & notBlanks,:) = predLaserOff;
                predictions(laserOn & notBlanks,:) = predLaserOn;
                blankResponses = [blanksLaserOff, blanksLaserOn];
                if respSign < 0
                    parameters(2,:) = -parameters(2,:);
                    parameters(4,:) = -parameters(4,:);
                    predictions = -predictions;
                    blankResponses = cellfun(@times, blankResponses, ...
                        num2cell(-ones(size(blankResponses))), ...
                        'UniformOutput', false);
                    curves = -curves;
                end
            end
            
            err = errors{1}(:,notBlanks)';
            respNoBlanks = resp(notBlanks,:);
            ind = ~isnan(err) & ~isnan(respNoBlanks);
            explVar = 1 - sum(err(ind).^2) / ...
                sum((respNoBlanks(ind)-mean(respNoBlanks(ind))).^2);
            tuning(iExp).crossValExplVar(iCell) = explVar;
            
            predNoBlanks = predictions(notBlanks,:);
            ind = ~isnan(respNoBlanks) & ~isnan(predNoBlanks);
            R2 = 1 - sum((respNoBlanks(ind) - predNoBlanks(ind)).^2) / ...
                sum((respNoBlanks(ind) - mean(respNoBlanks(ind))).^2);
            R2ToBlank = 1 - sum((respNoBlanks(ind) - predNoBlanks(ind)).^2) / ...
                sum((respNoBlanks(ind) - ...
                nanmean(reshape(resp(blanks,:),[],1))).^2);
            tuning(iExp).R2(iCell).comparedToMean = R2;
            tuning(iExp).R2(iCell).comparedToBlankResp = R2ToBlank;
            for c = 1:size(parameters,2)
%                 tuning(iExp).cell(iCell).cond(c).parameters = parameters(:,c);
%                 tuning(iExp).cell(iCell).cond(c).curve = curves(:,c);
%                 tuning(iExp).cell(iCell).cond(c).blankResponses = ...
%                     blankResponses{c};
                tuning(iExp).cond(c).cell(iCell).parameters = parameters(:,c);
                tuning(iExp).cond(c).cell(iCell).curve = curves(:,c);
                tuning(iExp).cond(c).cell(iCell).directions = directions(notBlanks);
                r = NaN(size(resp));
                r(conds == c) = resp(conds == c);
                r(blanks,:) = [];
                tuning(iExp).cond(c).cell(iCell).responses = r;
                tuning(iExp).cond(c).cell(iCell).blankResponses = ...
                    blankResponses{c};
            end
        
            if doSave == 1
                save(fullfile(folderResults, nonVisualSignal, ...
                    fileNames{iPars}), 'tuning', 'x')
            end
        end
        fprintf('\n')
    end
end

%% Add depth of each neuron to tuning structure (TODO)

%% Compare tuning curve fits (which parameters fixed)
fileNames = {'tuning_nothingFixed.mat', ...
    'tuning_prefDirFixed.mat', ...
    'tuning_prefDirSigmaFixed.mat', ...
    'tuning_prefDirSigmaDIFixed.mat', ...
    'tuning_prefDirAmpSigmaDIFixed.mat', ...
    'tuning_prefDirSigmaDIBaselineFixed.mat', ...
    'tuning_constantFit.mat'};
groups = {'-', 'D_p', 'D_p, sig', 'D_p, sig, DI', 'D_p, R_p, sigma, DI', ...
    'D_p, sigma, offset, DI', 'const'};

explVar = cell(size(fileNames));
for iPars = 1:length(fileNames)
    data = load(fullfile(folderResults, nonVisualSignal, fileNames{iPars}));
    tuning = data.tuning;
    for iExp = 1:length(tuning)
        explVar{iPars} = [explVar{iPars}; ...
            tuning(iExp).crossValExplVar];
    end
end

explVar = cell2mat(explVar);
g = size(explVar,2);
bins = -.175:.05:1;
figure
for i = 1:g
    for j = 1:g
        subplot(g,g,(i-1)*g+j)
        if i == j
            hist(explVar(:,i), bins)
            xlim([-.2 1])
        else
            plot(explVar(:,j),explVar(:,i),'k.')
            axis([-.2 1 -.2 1])
            axis square
        end
        if i == 1
            title(groups{j})
        end
        if j == 1
            ylabel(groups{i})
        end
        if i == g && j == ceil(g/2)
            xlabel('Explained variance')
        end
    end
end
figure
hold on
for i = 1:g
    n = hist(explVar(:,i),bins);
    plot(bins,n,'LineWidth',2)
end
legend(groups{1:end})
xlabel('Explained variance')
ylabel('# Neurons')

anova1(explVar(:,1:end-1));

%% Add which neurons are not tuned to tuning structure

data = load(fullfile(folderResults, nonVisualSignal, ...
    'tuning_prefDirSigmaDIFixed.mat'));
tuning = data.tuning;
x = data.x;
data = load(fullfile(folderResults, nonVisualSignal, ...
    'tuning_constantFit.mat'));
tunC = data.tuning;

for iExp = 1:length(tuning)
    tuned = NaN(length(tuning(iExp).cellIDs), 1);
    tuned(~isnan(tuning(iExp).crossValExplVar)) = 0;
    tuned(tuning(iExp).crossValExplVar > tunC(iExp).crossValExplVar) = 1;
    tuning(iExp).isTuned = tuned;
    ind = find(tuned == 0);
    for k = 1:length(ind)
        for c = 1:4
            par = tunC(iExp).cond(c).cell(ind(k)).parameters;
            tuning(iExp).cond(c).cell(ind(k)).parameters = ...
                par;
            tuning(iExp).cond(c).cell(ind(k)).curve = ...
                ones(1,length(x)) .* par;
        end
        tuning(iExp).crossValExplVar(ind(k)) = ...
            tunC(iExp).crossValExplVar(ind(k));
        tuning(iExp).R2(ind(k)) = ...
            tunC(iExp).R2(ind(k));
    end
end
save(fullfile(folderResults, nonVisualSignal, ...
    'tuning_prefDirSigmaDIFixed_isTuned.mat'), 'tuning', 'x');

%% Add line fits
data = load(fullfile(folderResults, nonVisualSignal, ...
    'tuning_prefDirSigmaDIFixed_isTuned.mat'));
tuning = data.tuning;
x = data.x;

for iExp = 1:length(tuning)
    if isempty(tuning(iExp).cond)
        continue
    end
    for iCell = 1:length(tuning(iExp).cond(1).cell)
        dat1 = tuning(iExp).cond(1).cell(iCell);
        if isempty(dat1.parameters)
            continue
        end
        dat2 = tuning(iExp).cond(2).cell(iCell);
        if tuning(iExp).isTuned(iCell) == 1
            lineFit.intercept = dat2.parameters(4) - dat1.parameters(4) * ...
                dat2.parameters(2) / dat1.parameters(2);
            lineFit.slope = dat2.parameters(2) / dat1.parameters(2);
        else
            lineFit.intercept = dat2.parameters - dat1.parameters;
            lineFit.slope = NaN;
        end
        tuning(iExp).lineFit(iCell).laserOff = lineFit;
        
        dat1 = tuning(iExp).cond(3).cell(iCell);
        dat2 = tuning(iExp).cond(4).cell(iCell);
        clear lineFit
        if tuning(iExp).isTuned(iCell) == 1
            lineFit.intercept = dat2.parameters(4) - dat1.parameters(4) * ...
                dat2.parameters(2) / dat1.parameters(2);
            lineFit.slope = dat2.parameters(2) / dat1.parameters(2);
        else
            lineFit.intercept = dat2.parameters - dat1.parameters;
            lineFit.slope = NaN;
        end
        tuning(iExp).lineFit(iCell).laserOn = lineFit;
    end
end
save(fullfile(folderResults, nonVisualSignal, ...
    'tuning_prefDirSigmaDIFixed_lineFit.mat'), 'tuning', 'x');

%% Create tuning null distribution (shuffle conditions)
draws = 200;
fixedPars = [1 3 5];
data = load(fullfile(folderResults, nonVisualSignal, ...
    'tuning_prefDirSigmaDIFixed_isTuned.mat'));
tuning = data.tuning;
if ~exist(fullfile(folderResults, nonVisualSignal), 'dir')
    mkdir(fullfile(folderResults, nonVisualSignal));
end
for iExp = 1:length(tuning)
    fprintf('  Dataset %d: %s %s exp.: %d\n', iExp, tuning(iExp).subject, ...
        tuning(iExp).date, tuning(iExp).exp);
    if isempty(tuning(iExp).cond)
        continue
    end
    nullBehaviour(iExp).subject = tuning(iExp).subject;
    nullBehaviour(iExp).date = tuning(iExp).date;
    nullBehaviour(iExp).exp = tuning(iExp).exp;
    nullBehaviour(iExp).cellIDs = tuning(iExp).cellIDs;
    nullBehaviour(iExp).isSuppressed = tuning(iExp).isSuppressed;
    nullBehaviour(iExp).isTuned = tuning(iExp).isTuned;
    
    nullLaser(iExp).subject = tuning(iExp).subject;
    nullLaser(iExp).date = tuning(iExp).date;
    nullLaser(iExp).exp = tuning(iExp).exp;
    nullLaser(iExp).cellIDs = tuning(iExp).cellIDs;
    nullLaser(iExp).isSuppressed = tuning(iExp).isSuppressed;
    nullLaser(iExp).isTuned = tuning(iExp).isTuned;
    
    nonVis = tuning(iExp).nonVisual;
    laserOn = repmat(tuning(iExp).laserOn, 1, size(nonVis,2));
    directions = tuning(iExp).directions;
    fprintf('    Cell (of %d):', length(tuning(iExp).cellIDs));
    for iCell = 1:length(tuning(iExp).cell)
        fprintf(' %d', iCell)
        if isnan(tuning(iExp).isTuned(iCell))
            continue
        end
        resp = tuning(iExp).cell(iCell).response;
        np = 5;
        if tuning(iExp).isTuned(iCell) == 0
            np = 1;
        end
        parametersBehaviour = NaN(draws, 4, np);
        parametersLaser = NaN(draws, 4, np);
        if tuning(iExp).isTuned(iCell) == 0
            pb1 = NaN(draws,2);
            pb2 = NaN(draws,2);
            parfor k = 1:draws
                p = randperm(numel(nonVis));
                cond = reshape(nonVis(p), size(nonVis));
                for c = 1:2
                    ind = cond == c & ~laserOn;
                    pb1(k,c) = nanmean(resp(ind));
                    ind = cond == c & laserOn;
                    pb2(k,c) = nanmean(resp(ind));
                end
            end
            parametersBehaviour(:,1:2) = pb1;
            parametersBehaviour(:,3:4) = pb2;
            
            pl1 = NaN(draws,2);
            pl2 = NaN(draws,2);
            parfor k = 1:draws
                p = randperm(numel(laserOn));
                lo = reshape(laserOn(p), size(laserOn));
                for c = 1:2
                    ind = conditions == c & ~lo;
                    pl1(k,c) = nanmean(resp(ind));
                    ind = conditions == c & lo;
                    pl2(k,c) = nanmean(resp(ind));
                end
            end
            parametersLaser(:,1:2) = pl1;
            parametersLaser(:,3:4) = pl2;
        else
            if tuning(iExp).isSuppressed(iCell) > 0
                resp = -resp;
            end
            pb1 = NaN(draws,5);
            pb2 = NaN(draws,5);
            pb3 = NaN(draws,5);
            pb4 = NaN(draws,5);
%             parfor k = 1:draws
            for k = 1:draws
                p = randperm(numel(nonVis));
                cond = reshape(nonVis(p), size(nonVis));
                [parsOff, parsOn] = gratings.fitTuningCurveConditionsAndLaser( ...
                    resp, directions, cond, laserOn, fixedPars, runFitTwice);
                pb1(k,:) = parsOff(1:5);
                pb2(k,:) = parsOff(1:5)+parsOff(6:10);
                pb3(k,:) = parsOn(1:5);
                pb4(k,:) = parsOn(1:5)+parsOn(6:10);
            end
            parametersBehaviour(:,1,:) = pb1;
            parametersBehaviour(:,2,:) = pb2;
            parametersBehaviour(:,3,:) = pb3;
            parametersBehaviour(:,4,:) = pb4;
            
            pl1 = NaN(draws,5);
            pl2 = NaN(draws,5);
            pl3 = NaN(draws,5);
            pl4 = NaN(draws,5);
            parfor k = 1:draws
                p = randperm(numel(laserOn));
                lo = reshape(laserOn(p), size(laserOn));
                [parsOff, parsOn] = gratings.fitTuningCurveConditionsAndLaser( ...
                    resp, directions, nonVis, lo, fixedPars, runFitTwice);
                pl1(k,:) = parsOff(1:5);
                pl2(k,:) = parsOff(1:5)+parsOff(6:10);
                pl3(k,:) = parsOn(1:5);
                pl4(k,:) = parsOn(1:5)+parsOn(6:10);
            end
            parametersLaser(:,1,:) = pl1;
            parametersLaser(:,2,:) = pl2;
            parametersLaser(:,3,:) = pl3;
            parametersLaser(:,4,:) = pl4;
            
            if tuning(iExp).isSuppressed(iCell) > 0
                parametersBehaviour(:,:,2) = -parametersBehaviour(:,:,2);
                parametersBehaviour(:,:,4) = -parametersBehaviour(:,:,4);
                parametersLaser(:,:,2) = -parametersLaser(:,:,2);
                parametersLaser(:,:,4) = -parametersLaser(:,:,4);
            end
        end
        for c = 1:4
            nullBehaviour(iExp).cond(c).cell(iCell).parameters = squeeze(parametersBehaviour(:,c,:));
            nullLaser(iExp).cond(c).cell(iCell).parameters = squeeze(parametersLaser(:,c,:));
        end
        
        null = nullBehaviour;
        save(fullfile(folderResults, nonVisualSignal, ...
            'nullTunning_prefDirSigmaDIFixed_behaviour.mat'), 'null')
        
        null = nullLaser;
        save(fullfile(folderResults, nonVisualSignal, ...
            'nullTunning_prefDirSigmaDIFixed_laser.mat'), 'null')
    end
    fprintf('\n')
end

%% Plot tuning curves for each neuron
data = load(fullfile(folderResults, nonVisualSignal, ...
    'tuning_prefDirSigmaDIFixed_lineFit.mat'));
tuning = data.tuning;
x = data.x;
laserCond = [false, true];
conditions = {'laserOff', 'laserOn'};
behaviour = {'small pupil', 'large pupil'};
lin = {'-','-'};
cols = {'k', 'r', 'c', 'm'};
col2 = lines(1);
fPlots = fullfile(folderResults, ['plots_' nonVisualSignal], ...
    'tuningCurves');
if ~exist(fPlots, 'dir')
    mkdir(fPlots);
end
for iExp = 1:length(tuning)
    if isempty(tuning(iExp).cond)
        continue
    end
    fPlots = fullfile(folderResults, ['plots_' nonVisualSignal], ...
        'tuningCurves', [tuning(iExp).subject '_' tuning(iExp).date ...
        '_' num2str(tuning(iExp).exp)]);
    if ~exist(fPlots, 'dir')
        mkdir(fPlots);
    end
    for iCell = [2, 4, 7, 34, 35, 51, 95]%1:length(tuning(iExp).cellIDs)
        if isempty(tuning(iExp).cond(1).cell(iCell).parameters)
            continue
        end
        figure('Position', [3 715 1915 385])
        ax = [0 0];
        for l = 1:2
            subplot(1,4,l)
            hold on
            h = [0 0];
            for c = 1:2
                h(c) = plot(x, tuning(iExp).cond(c+(l-1)*2).cell(iCell).curve, lin{c}, ...
                    'Color', cols{c+(l-1)*2}, 'LineWidth',2);
                r = tuning(iExp).cell(iCell).response;
                r(tuning(iExp).laserOn~=laserCond(l) | tuning(iExp).nonVisual~=c) = NaN;
                r(end+1,:) = nanmean(r(tuning(iExp).directions==0,:),1);
                errorbar([tuning(iExp).directions;360], nanmean(r,2), ...
                    nanstd(r,0,2) ./ sqrt(sum(~isnan(r),2)), 'o', ...
                    'Color', cols{c+(l-1)*2}, 'MarkerFaceColor', cols{c+(l-1)*2})
                m = mean(tuning(iExp).cond(c+(l-1)*2).cell(iCell).blankResponses);
                s = std(tuning(iExp).cond(c+(l-1)*2).cell(iCell).blankResponses) / ...
                    length(tuning(iExp).cond(c+(l-1)*2).cell(iCell).blankResponses);
                fill(x([1 end end 1]), [[1 1].*(m+s), [1 1].*(m-s)], 'k', ...
                    'FaceColor', cols{c+(l-1)*2}, 'EdgeColor', 'none', ...
                    'FaceAlpha', 0.3)
                plot(x([1 end]), [m m], '--', 'Color', cols{c+(l-1)*2}, 'LineWidth', 2)
            end
            set(gca,'XTick',0:90:360)
            legend(h, tuning(iExp).conditionLabels((1:2)+(l-1)*2))
            legend('boxoff')
            title(sprintf('EV = %.2f, R2_{stimResp} = %.2f, R2_{tuning} = %.2f', ...
                tuning(iExp).crossValExplVar(iCell), ...
                tuning(iExp).R2(iCell).comparedToBlankResp, ...
                tuning(iExp).R2(iCell).comparedToMean))
            xlim([-10 370])
            xlabel('Direction (in degrees)')
            ylabel('Firing rate (spikes/s)')
            ax(l) = gca;
        end
        linkaxes(ax, 'y')
        
        subplot(1,4,3)
        hold on
        maxi = 0;
        for c = 1:2
            r = tuning(iExp).cell(iCell).response;
            r(tuning(iExp).nonVisual~=c) = NaN;
            r1 = r(tuning(iExp).laserOn==laserCond(1),:);
            r1 = nanmean(r1, 2);
            r2 = r(tuning(iExp).laserOn==laserCond(2),:);
            r2 = nanmean(r2, 2);
            plot(r1, r2, 'ok', 'MarkerFaceColor', cols{c}, ...
                'MarkerEdgeColor', 'none')
            maxi = max([maxi;r1;r2]);
        end
        plot([0 maxi],[0 maxi],'k:')
        if maxi > 0
            axis([0 maxi 0 maxi])
        end
        axis square
        xlabel(conditions{1})
        ylabel(conditions{2})
        legend(behaviour, 'location', 'NorthWest')
        legend('boxoff')
        
        subplot(1,4,4)
        hold on
        h = [0 0];
        r1 = cell(1,2);
        r2 = cell(1,2);
        maxi1 = [0 0];
        for l = 1:2
            r = tuning(iExp).cell(iCell).response;
            r(tuning(iExp).laserOn~=laserCond(l),:) = NaN;
            r1{l} = r;
            r1{l}(tuning(iExp).nonVisual~=1) = NaN;
            r1{l} = nanmean(r1{l},2);
            r2{l} = r;
            r2{l}(tuning(iExp).nonVisual~=2) = NaN;
            r2{l} = nanmean(r2{l},2);
            maxi1(l) = max(r1{l});
            r1{l} = r1{l} ./ maxi1(l);
            r2{l} = r2{l} ./ maxi1(l);
        end
        mini = min([cat(1,r1{:}); cat(1,r2{:})]);
        maxi = max([cat(1,r1{:}); cat(1,r2{:})]);
        if isnan(mini) || isinf(mini)
            close gcf
            continue
        end
        for l = 1:2
            plot(r1{l}, r2{l}, 'o', 'Color', cols{1+(l-1)*2}, ...
                'MarkerFaceColor', cols{1+(l-1)*2})
            intercept = tuning(iExp).lineFit(iCell) ...
                .(conditions{l}).intercept / maxi1(l);
            slope = tuning(iExp).lineFit(iCell).(conditions{l}).slope;
            if isnan(slope)
                slope = 1;
            end
%             if tuning(iExp).isSuppressed(iCell) == 1
%                 intercept = -intercept;
%             end
            h(l) = plot([mini maxi],[mini maxi].*slope+intercept, ...
                'Color',cols{1+(l-1)*2},'LineWidth',2);
        end
        plot([mini maxi],[mini maxi],'k:')
        axis([mini maxi mini maxi])
        axis square
        xlabel(behaviour{1})
        ylabel(behaviour{2})
        legend(h, conditions, 'location', 'NorthWest')
        legend('boxoff')
        set(gca,'box','off')
        title(['Unit ' num2str(iCell)])
%         title(sprintf('y = %.2f x + %.2f', slope, intercept))
          
%         fig = gcf;
%         fig.PaperPositionMode = 'auto';
%         print(fullfile(fPlots, sprintf('tuningCurve%03d.jpg', iCell)), ...
%             '-djpeg','-r0')
%         close gcf
    end
end

%% CONTINUE: (add depth for depth plot) Plot mean responses during laser off vs. laser on for population
laserConds = [false true false true]; % [off on off on]
pupilConds = [1 1 2 2]; % [small small large large]

data = load(fullfile(folderResults, nonVisualSignal, ...
    'tuning_prefDirSigmaDIFixed_isTuned.mat'));
tuning = data.tuning;

responses = [];

for iExp = 1:length(tuning)
    if isempty(tuning(iExp).cond)
        continue
    end
    for iCell = 1:length(tuning(iExp).cellIDs)
        resp = tuning(iExp).cell(iCell).response;
        means = zeros(1,4);
        for c = 1:4
            use = tuning(iExp).laserOn==laserConds(c) & ...
                tuning(iExp).nonVisual==pupilConds(c);
            r = resp;
            r(~use) = NaN;
            means(c) = nanmean(nanmean(r,2),1);
        end
        responses(end+1,:) = means;
    end
end

% Scatter: laser off vs laser on
cols = 'kr';
differences = cell(1,2);
h = [0 0 0 0];
figure
hold on
for b = 1:2
    diffs = responses{b} - responses{b+2}; % laser off - laser on
    nullDiffs = nullResponses{b+2} - nullResponses{b};
    confInts = prctile(nullDiffs, [2.5 97.5], 2);
    isSign = diffs < confInts(:,1) | diffs > confInts(:,2);
    
    h(2*b-1) = scatter(responses{b}(isSign), responses{b+2}(isSign), 30, cols(b), ...
        'filled', 'MarkerFaceAlpha', .3);
    h(2*b) = scatter(responses{b}(~isSign), responses{b+2}(~isSign), 10, cols(b), ...
        'filled', 'MarkerFaceAlpha', .3);
    
    differences{1} = [differences{1}; diffs(isSign)];
    differences{2} = [differences{2}; diffs(~isSign)];
end
maxi = max(cat(1,responses{:})) * 1.05;
plot([0 maxi], [0 maxi], 'k')
axis([0 maxi 0 maxi])
axis square
xlabel('laser off')
ylabel('laser on')
legend(h, 'small pupil, p < 0.05', 'small pupil, p \geq 0.05', ...
    'large pupil, p < 0.05', 'large pupil, p \geq 0.05', ...
    'Location', 'NorthWest')
legend('boxoff')
title(sprintf('n = %d', length(responses{1})*2))

% Histogram: difference of responses (laser off - laser on)
figure
% maxi = ceil(max(abs(cat(1, differences{:}))) * 2) / 2;
bins = -10+.25 : .5 : 10-.25;
n1 = hist(differences{1}, bins)';
n2 = hist(differences{2}, bins)';
bar(bins, [n1, n2], 'stacked')
colormap gray
xlabel('\Deltamean response (laser off - laser on)')
ylabel('# Neurons')
set(gca, 'box', 'off')
legend('p < 0.05', 'p \geq 0.05')
legend('boxoff')
m = median(cat(1,responses{3:4}) - cat(1,responses{1:2}));
p = signrank(cat(1,responses{3:4}), cat(1,responses{1:2}));
title(sprintf('Median diff: %.3f (signed rank test: p = %.4f)', m, p))

%% Plot difference index for minimum, maximum, ... of tuning curve for population
data = load(fullfile(folderResults, nonVisualSignal, ...
    'tuning_prefDirSigmaDIFixed_isTuned.mat'));
tuning = data.tuning;
degrees = data.x;
data = load(fullfile(folderResults, nonVisualSignal, ...
    'nullTuning_prefDirSigmaDIFixed_behaviour.mat'));
nullBeh = data.null;
data = load(fullfile(folderResults, nonVisualSignal, ...
    'nullTuning_prefDirSigmaDIFixed_laser.mat'));
nullLaser = data.null;

minis = cell(1,4);
maxis = cell(1,4);
stimMeans = cell(1,4);
grayScreen = cell(1,4);
nullBehMinis = cell(1,4);
nullBehMaxis = cell(1,4);
nullBehStimMeans = cell(1,4);
nullLaserMinis = cell(1,4);
nullLaserMaxis = cell(1,4);
nullLaserStimMeans = cell(1,4);
isSuppr = [];

for iExp = 1:length(tuning)
    fprintf('Dataset %d: %s %s exp.: %d\n', iExp, tuning(iExp).subject, ...
        tuning(iExp).date, tuning(iExp).exp);
    if isempty(tuning(iExp).cond)
        continue
    end
    for iCell = 1:length(tuning(iExp).cellIDs)
        if isnan(tuning(iExp).isTuned(iCell))
            continue
        end
        for c = 1:4
            pars = tuning(iExp).cond(c).cell(iCell).parameters;
            curve = tuning(iExp).cond(c).cell(iCell).curve;
            if length(pars) == 1
                m = pars;
                mi = NaN;
                ma = NaN;
            else
                m = mean(curve);
                oris = mod(pars(1) + [0 90 180], 360);
                sr = gratings.orituneWrappedConditions(pars, oris);
                ma = sr(1);
                if sr(1)-sr(2)>0
                    [mi, ind] = min(sr(2:3));
                else % suppressed
                    [mi, ind] = max(sr(2:3));
                end
                ind = ind+1;
            end
            minis{c}(end+1,1) = mi;
            maxis{c}(end+1,1) = ma;
            stimMeans{c}(end+1,1) = m;
            grayScreen{c}(end+1,1) = mean(tuning(iExp).cond(c) ...
                .cell(iCell).blankResponses);
            
            pars = nullBeh(iExp).cond(c).cell(iCell).parameters;
            if size(pars,2) == 1
                m = pars;
                mi = NaN(size(m));
                ma = NaN(size(m));
            else
                sr = NaN(size(pars,1), 3);
                curves = NaN(size(pars,1), length(degrees));
                for p = 1:size(pars,1)
                    oris = mod(pars(p,1) + [0 90 180], 360);
                    sr(p,:) = gratings.orituneWrappedConditions(pars(p,:), oris);
                    curves(p,:) = gratings.orituneWrappedConditions(pars(p,:), degrees);
                end
                m = mean(curves,2);
                ma = sr(:,1);
                mi = sr(:,ind);
            end
            nullBehMinis{c}(end+1,:) = mi;
            nullBehMaxis{c}(end+1,:) = ma;
            nullBehStimMeans{c}(end+1,:) = m;
            
            pars = nullLaser(iExp).cond(c).cell(iCell).parameters;
            if size(pars,2) == 1
                m = pars;
                mi = NaN(size(m));
                ma = NaN(size(m));
            else
                sr = NaN(size(pars,1), 3);
                curves = NaN(size(pars,1), length(degrees));
                for p = 1:size(pars,1)
                    oris = mod(pars(p,1) + [0 90 180], 360);
                    sr(p,:) = gratings.orituneWrappedConditions(pars(p,:), oris);
                    curves(p,:) = gratings.orituneWrappedConditions(pars(p,:), degrees);
                end
                m = mean(curves,2);
                ma = sr(:,1);
                mi = sr(:,ind);
            end
            nullLaserMinis{c}(end+1,:) = mi;
            nullLaserMaxis{c}(end+1,:) = ma;
            nullLaserStimMeans{c}(end+1,:) = m;
        end
    end
    isSuppr = [isSuppr; tuning(iExp).isSuppressed];
end

mods = cell(1,4);
nullBehMods = cell(1,4);
nullLaserMods = cell(1,4);
j = isSuppr==1;
for c = 1:4
    mods{c} = abs(maxis{c} - minis{c});
    nullBehMods{c} = abs(nullBehMaxis{c} - nullBehMinis{c});
    nullLaserMods{c} = abs(nullLaserMaxis{c} - nullLaserMinis{c});
    
    m = minis{c}(j);
    minis{c}(j) = maxis{c}(j);
    maxis{c}(j) = m;    
end

% Plot scatters: diff. index of various measures against each other
measures = {grayScreen, mods, stimMeans};
labels = {'Gray screen', 'Tuning modulation', 'Mean Response'};
modFuns = @(a,b) (b-a)./(abs(a)+abs(b));
funLabels = 'Normalised difference';
ticks = -1:.5:1;

for m1 = 1:length(measures)
    d1 = modFuns(measures{m1}{1}, measures{m1}{2});
    for m2 = m1+1:length(measures)
        d2 = modFuns(measures{m2}{1}, measures{m2}{2});
        
        figure
        plot(d1, d2, 'k.')
        xlabel(labels{m1})
        ylabel(labels{m2})
        title('Difference: large vs. small pupil')
        axis square
        axis([-1 1 -1 1])
        set(gca, 'box', 'off', 'XTick', ticks, 'YTick', ticks)
    end
end





measures = {minis, maxis, mods, stimMeans};
nullBehMeasures = {nullBehMinis, nullBehMaxis, nullBehMods, nullBehStimMeans};
nullLaserMeasures = {nullLaserMinis, nullLaserMaxis, nullLaserMods, nullLaserStimMeans};
labels = {'Minimum', 'Maximum', 'Tuning modulation', 'Mean'};
modFuns = @(a,b) (b-a)./(abs(a)+abs(b));
funLabels = 'Diff-index (large-small)/(large+small)';
binSizes = 0.05;
mini = -.975;
maxi = .975;

% Plot scatters: diff. index - laser off vs. on
for m = 1:length(measures)
    diffOff = modFuns(measures{m}{1}, measures{m}{2});
    diffOn = modFuns(measures{m}{3}, measures{m}{4});
    pseudoBehOff = modFuns(nullBehMeasures{m}{1}, nullBehMeasures{m}{2});
    pseudoBehOn = modFuns(nullBehMeasures{m}{3}, nullBehMeasures{m}{4});
    pseudoLaserOff = modFuns(nullLaserMeasures{m}{1}, nullLaserMeasures{m}{2});
    pseudoLaserOn = modFuns(nullLaserMeasures{m}{3}, nullLaserMeasures{m}{4});
    
    confIntBehOff = prctile(pseudoBehOff, [2.5 97.5], 2);
    signBeh = diffOff < confIntBehOff(:,1) | diffOff > confIntBehOff(:,2);
    
    diffs = diffOff - diffOn;
    pseudoDiffsLaser = pseudoLaserOff - pseudoLaserOn;
    confIntLaser = prctile(pseudoDiffsLaser, [2.5 97.5], 2);
    signLaser = diffs < confIntLaser(:,1) | diffs > confIntLaser(:,2);
    
    n = [0 0 0 0];
    h = [0 0 0 0];
    figure('Position', [635 685 660 420])
    hold on
    fill([-1 1 1 -1],[-1 1 -1 1], 'k', 'EdgeColor', 'none', ...
        'FaceColor', 'k', 'FaceAlpha', .1)
    % (1) significant for laser (large dots)
    % (1a) not significant for pupil when laser off (black dots)
    ind = signLaser & ~signBeh;
    n(1) = sum(ind);
    h(1) = plot(diffOff(ind), diffOn(ind), 'o', 'MarkerSize', 8, ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w');
    % (1b) significant for pupil when laser off (red dots)
    ind = signLaser & signBeh;
    n(2) = sum(ind);
    h(2) = plot(diffOff(ind), diffOn(ind), 'o', 'MarkerSize', 8, ...
        'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'w');
    % (2) not significant for laser (small dots)
    % (2a) not significant for pupil when laser off (black dots)
    ind = ~signLaser & ~signBeh;
    n(3) = sum(ind);
    h(3) = plot(diffOff(ind), diffOn(ind), 'o', 'MarkerSize', 5, ...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k');
    % (2b) significant for pupil when laser off (red dots)
    ind = ~signLaser & signBeh;
    n(4) = sum(ind);
    h(4) = plot(diffOff(ind), diffOn(ind), 'o', 'MarkerSize', 5, ...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r');
    
    plot([-1 1], [-1 1], 'k:')
    
    legend(h, ['p_{laser}<.05, p_{pupil}\geq.05 (n=' num2str(n(1)) ')'], ...
        ['p_{laser}<.05, p_{pupil}<.05 (n=' num2str(n(2)) ')'], ...
        ['p_{laser}\geq.05, p_{pupil}\geq.05 (n=' num2str(n(3)) ')'], ...
        ['p_{laser}\geq.05, p_{pupil}<.05 (n=' num2str(n(4)) ')'], ...
        'Location', 'northeastoutside')
    axis([-1 1 -1 1])
    axis square
    xlabel('laser off')
    ylabel('laser on')
    title(['\Delta' sprintf('%s: off - on = %.3f (p = %.3f)', labels{m}, ...
        nanmedian(diffs), signrank(diffOff,diffOn))])
end

%% Plot delta-alpha (amplitude) and -beta (offset) for population
data = load(fullfile(folderResults, nonVisualSignal, ...
    'tuning_prefDirSigmaDIFixed_isTuned.mat'));
tuning = data.tuning;
data = load(fullfile(folderResults, nonVisualSignal, ...
    'nullTuning_prefDirSigmaDIFixed_behaviour.mat'));
nullBeh = data.null;
data = load(fullfile(folderResults, nonVisualSignal, ...
    'nullTuning_prefDirSigmaDIFixed_laser.mat'));
nullLaser = data.null;

% (1) Tuned neurons
deltaAlpha_off = []; % amplitude differences (large - small pupil), laser off
deltaBeta_off = []; % offset differences (large - small pupil), laser off
deltaAlpha_on = []; % amplitude differences (large - small pupil), laser on
deltaBeta_on = []; % offset differences (large - small pupil), laser on
nullDeltaAlpha_off = []; % behavioural tags randomised
nullDeltaBeta_off = [];
nullAlphaLaserOff = []; % laser tags randomised
nullAlphaLaserOn = [];
nullBetaLaserOff = [];
nullBetaLaserOn = [];
isSuppressed =[];
for iExp = 1:length(tuning)
    if isempty(tuning(iExp).cond)
        continue
    end
    tuned = tuning(iExp).isTuned == 1;
    pars1 = cat(2, tuning(iExp).cond(1).cell(tuned).parameters);
    pars2 = cat(2, tuning(iExp).cond(2).cell(tuned).parameters);
    pars3 = cat(2, tuning(iExp).cond(3).cell(tuned).parameters);
    pars4 = cat(2, tuning(iExp).cond(4).cell(tuned).parameters);
    nullPars1 = cat(3, nullBeh(iExp).cond(1).cell(tuned).parameters);
    nullPars2 = cat(3, nullBeh(iExp).cond(2).cell(tuned).parameters);
    nullLaser1 = cat(3, nullLaser(iExp).cond(1).cell(tuned).parameters);
    nullLaser2 = cat(3, nullLaser(iExp).cond(2).cell(tuned).parameters);
    nullLaser3 = cat(3, nullLaser(iExp).cond(3).cell(tuned).parameters);
    nullLaser4 = cat(3, nullLaser(iExp).cond(4).cell(tuned).parameters);
    
    normPars1 = abs(sum(pars1([2 4],:), 1))';
    normPars3 = abs(sum(pars3([2 4],:), 1))';
    deltaAlpha_off = [deltaAlpha_off; (pars2(2,:) - pars1(2,:))' ./ normPars1];
    deltaBeta_off = [deltaBeta_off; (pars2(4,:) - pars1(4,:))' ./ normPars1];
    deltaAlpha_on = [deltaAlpha_on; (pars4(2,:) - pars3(2,:))' ./ normPars3];
    deltaBeta_on = [deltaBeta_on; (pars4(4,:) - pars3(4,:))' ./ normPars3];
    
    nullDeltaAlpha_off = [nullDeltaAlpha_off; bsxfun(@rdivide, ...
        squeeze(nullPars2(:,2,:) - nullPars1(:,2,:))', normPars1)];
    nullDeltaBeta_off = [nullDeltaBeta_off; bsxfun(@rdivide, ...
        squeeze(nullPars2(:,4,:) - nullPars1(:,4,:))', normPars1)];
    
    nullAlphaLaserOff = [nullAlphaLaserOff; bsxfun(@rdivide, ...
        squeeze(nullLaser2(:,2,:) - nullLaser1(:,2,:))', normPars1)];
    nullBetaLaserOff = [nullBetaLaserOff; bsxfun(@rdivide, ...
        squeeze(nullLaser2(:,4,:) - nullLaser1(:,4,:))', normPars1)];
    nullAlphaLaserOn = [nullAlphaLaserOn; bsxfun(@rdivide, ...
        squeeze(nullLaser4(:,2,:) - nullLaser3(:,2,:))', normPars3)];
    nullBetaLaserOn = [nullBetaLaserOn; bsxfun(@rdivide, ...
        squeeze(nullLaser4(:,4,:) - nullLaser3(:,4,:))', normPars3)];
    
    isSuppressed = [isSuppressed; tuning(iExp).isSuppressed(tuned)];
end
ind = isSuppressed == 1;
deltaAlpha_off(ind) = -deltaAlpha_off(ind);
deltaBeta_off(ind) = -deltaBeta_off(ind);
deltaAlpha_on(ind) = -deltaAlpha_on(ind);
deltaBeta_on(ind) = -deltaBeta_on(ind);
nullDeltaAlpha_off(ind,:) = -nullDeltaAlpha_off(ind,:);
nullDeltaBeta_off(ind,:) = -nullDeltaBeta_off(ind,:);
nullAlphaLaserOff = -nullAlphaLaserOff;
nullAlphaLaserOn = -nullAlphaLaserOn;
nullBetaLaserOff = -nullBetaLaserOff;
nullBetaLaserOn = -nullBetaLaserOn;

confInt = prctile(nullDeltaAlpha_off, [2.5 97.5], 2);
signAlpha = deltaAlpha_off < confInt(:,1) | deltaAlpha_off > confInt(:,2);
confInt = prctile(nullDeltaBeta_off, [2.5 97.5], 2);
signBeta = deltaBeta_off < confInt(:,1) | deltaBeta_off > confInt(:,2);

% Scatterplots
% deltaAlpha_off vs deltaAlpha_on
figure
hold on
scatter(deltaAlpha_off(~signAlpha), deltaAlpha_on(~signAlpha), ...
    30, [.7 .7 .7], 'filled', 'MarkerEdgeColor', 'none')
scatter(deltaAlpha_off(signAlpha), deltaAlpha_on(signAlpha), ...
    30, 'k', 'filled', 'MarkerEdgeColor', 'none')
plot([-3 22], [-3 22], 'k:')
alpha(.5)
legend('p \geq 0.05', 'p< 0.05', 'Location', 'NorthWest')
xlim([-3 22])
ylim([-3 22])
axis square
xlabel('laser off')
ylabel('laser on')
title('\Delta\alpha')

figure
hold on
scatter(deltaBeta_off(~signBeta), deltaBeta_on(~signBeta), ...
    30, [.7 .7 .7], 'filled', 'MarkerEdgeColor', 'none')
scatter(deltaBeta_off(signBeta), deltaBeta_on(signBeta), ...
    30, 'k', 'filled', 'MarkerEdgeColor', 'none')
plot([-10 5], [-10 5], 'k:')
alpha(.5)
legend('p \geq 0.05', 'p< 0.05', 'Location', 'NorthWest')
xlim([-10 5])
ylim([-10 5])
axis square
xlabel('laser off')
ylabel('laser on')
title('\Delta\beta')

%% Plot slopes and intercepts (population)
minR2 = 0.3;
conditions = {'laserOff','laserOn'};

data = load(fullfile(folderResults, nonVisualSignal, ...
    'tuning_prefDirSigmaDIFixed_lineFit.mat'));
% data = load(fullfile(folderResults, nonVisualSignal, ...
%     'tuning_prefDirSigmaDIFixed_bootstrap.mat'));
tuning = data.tuning;

% plot intercepts and slopes
% normalize all parameters so that max. response at small pupil/no running
% is 1
slopes = [];
intercepts = [];
% prctlSlopes = [];
% prctlIntercepts = [];
maxRespDiffs = [];
R2 = [];
isSuppr = [];
for iExp = 1:length(tuning)
    data = tuning(iExp);
    if isempty(data.cond)
        continue
    end
    neurons = find(~cellfun(@isempty, {data.cond(1).cell.parameters}));
    for j = 1:length(neurons)
        iCell = neurons(j);
        sl = [NaN NaN];
        intc = [NaN NaN];
        rDiffs = [NaN NaN];
        isSuppr(end+1,1) = data.isSuppressed(iCell);
        R2(end+1,1) = data.R2(iCell).comparedToMean;
        if R2(end) < minR2
            slopes = [slopes; sl];
            intercepts = [intercepts; intc];
            maxRespDiffs = [maxRespDiffs; rDiffs];
            continue
        end
        for l = 1:2
            maxR = max(nanmean(data.cond(1+(l-1)*2).cell(iCell).responses, 2));
            sl(l) = data.lineFit(iCell).(conditions{l}).slope;
            intc(l) = data.lineFit(iCell).(conditions{l}).intercept / maxR;
            if isSuppr(end) == 1
                minR = min(nanmean(data.cond(1+(l-1)*2).cell(iCell).responses, 2));
                rDiffs(l) = -(minR * sl(l) + intc(l) - minR);
            else
                rDiffs(l) = sl(l) + intc(l) - 1;
            end
        end
        slopes = [slopes; sl];
        intercepts = [intercepts; intc];
        maxRespDiffs = [maxRespDiffs; rDiffs];
    end
end

% Plot max. response differences
bins = -1:.1:3;
figure
hist(maxRespDiffs(:,1), bins);
h = findobj(gca,'Type','patch');
h.FaceColor = 'k';
set(gca, 'box', 'off')
xlim(bins([1 end]))
xlabel('Max. resp. diff.')
title('Laser off')
ax1 = gca;
figure
hist(maxRespDiffs(:,2), bins);
h = findobj(gca,'Type','patch');
h.FaceColor = 'k';
set(gca, 'box', 'off')
xlim(bins([1 end]))
xlabel('Max. resp. diff.')
title('Laser on')
ax2 = gca;
linkaxes([ax1 ax2], 'y')
figure
mrd = maxRespDiffs;
mrd(mrd>bins(end)) = bins(end);
mrd(mrd<bins(1)) = bins(1);
plot(mrd(:,1), mrd(:,2), 'k.')
hold on
plot(bins([1 end]), bins([1 end]), ':k')
axis square
set(gca, 'box', 'off')
axis([bins([1 end]), bins([1 end])])
xlabel('Max. resp. diff. - laser off')
ylabel('Max. resp. diff. - laser on')

%% Plot characteristics of each neuron against each other (ggobi)
data = load(fullfile(folderResults, nonVisualSignal, ...
    'tuning_prefDirSigmaDIFixed_lineFit.mat'));
tuning = data.tuning;
x = data.x;
baseline_small = [];
baseline_large = [];
meanStim_small = [];
meanStim_large = [];
maxStim_small = [];
maxStim_large = [];
modulation_small = [];
modulation_large = [];
prefStim_small = [];
prefStim_large = [];
suppressed_mean = [];
suppressed_max = [];
for iExp = 1:length(tuning)
    if isempty(tuning(iExp).cond)
        continue
    end
    for iCell = 1:length(tuning(iExp).cellIDs)
        if isempty(tuning(iExp).cond(1).cell(iCell).parameters)
            continue
        end
        baseline_small(end+1,1) = mean(tuning(iExp).cond(1).cell(iCell).blankResponses);
        baseline_large(end+1,1) = mean(tuning(iExp).cond(2).cell(iCell).blankResponses);
        meanResp_small = nanmean(tuning(iExp).cond(1).cell(iCell).responses, 2);
        meanResp_large = nanmean(tuning(iExp).cond(2).cell(iCell).responses, 2);
        meanStim_small(end+1,1) = nanmean(meanResp_small);
        meanStim_large(end+1,1) = nanmean(meanResp_large);
        maxStim_small(end+1,1) = max(meanResp_small);
        maxStim_large(end+1,1) = max(meanResp_large);
        [~,oriInd] = min(abs(x-tuning(iExp).cond(1).cell(iCell).parameters(1)));
        modulation_small(end+1,1) = abs(diff(tuning(iExp).cond(1).cell(iCell) ...
            .curve([mod(oriInd+179,360)+1 oriInd])));
        modulation_large(end+1,1) = abs(diff(tuning(iExp).cond(2).cell(iCell) ...
            .curve([mod(oriInd+179,360)+1 oriInd])));
        prefStim_small(end+1,1) =  sum(tuning(iExp).cond(1).cell(iCell) ...
            .parameters([2 4]));
        prefStim_large(end+1,1) = sum(tuning(iExp).cond(2).cell(iCell) ...
            .parameters([2 4]));
        suppressed_mean(end+1,1) = (meanStim_small(end)-baseline_small(end));
        [~,ind] = max(abs(meanResp_small - baseline_small(end)));
        suppressed_max(end+1,1) = (meanResp_small(ind) - baseline_small(end));
    end
end
prefStim_small(prefStim_small<0) = 0;
prefStim_large(prefStim_large<0) = 0;

ratio_baseline = log10((baseline_large+.5) ./ (baseline_small+.5));
ratio_modulation = log10((modulation_large+.5) ./ (modulation_small+.5));
ratio_meanStim = log10((meanStim_large+.5) ./ (meanStim_small+.5));
ratio_prefStim = log10((prefStim_large+.5) ./ (prefStim_small+.5));

diff_baseline = baseline_large - baseline_small;
diff_modulation = modulation_large - modulation_small;
diff_meanStim = meanStim_large - meanStim_small;
diff_prefStim = prefStim_large - prefStim_small;

mean_baseline = mean([baseline_large baseline_small],2);
mean_modulation = mean([modulation_large modulation_small],2);
mean_meanStim = mean([meanStim_large meanStim_small],2);
mean_prefStim = mean([prefStim_large prefStim_small],2);

n=6;
figure('Position',[1675 800 2385 315])
subplot(1,n,1)
plot(ratio_baseline, ratio_modulation, 'k.')
xlabel('Baseline (large / small)')
ylabel('Tuning modulation (large / small)')
title('Log ratios')
subplot(1,n,2)
plot(ratio_baseline, ratio_meanStim, 'k.')
xlabel('Baseline (large / small)')
ylabel('Mean stim resp (large / small)')
subplot(1,n,3)
plot(ratio_baseline, ratio_prefStim, 'k.')
xlabel('Baseline (large / small)')
ylabel('Pref stim resp (large / small)')
subplot(1,n,4)
plot(ratio_modulation, ratio_prefStim, 'k.')
xlabel('Tuning modulation (large / small)')
ylabel('Pref stim resp (large / small)')
subplot(1,n,5)
plot(ratio_meanStim, ratio_prefStim, 'k.')
xlabel('Mean stim resp (large / small)')
ylabel('Pref stim resp (large / small)')
subplot(1,n,6)
plot(ratio_meanStim, ratio_modulation, 'k.')
xlabel('Mean stim resp (large / small)')
ylabel('Tuning modulation (large / small)')

n=6;
figure('Position',[1675 410 2385 315])
subplot(1,n,1)
plot(diff_baseline, diff_modulation, 'k.')
xlabel('Baseline (large - small)')
ylabel('Tuning modulation (large - small)')
title('Differences')
subplot(1,n,2)
plot(diff_baseline, diff_meanStim, 'k.')
xlabel('Baseline (large - small)')
ylabel('Mean stim resp (large - small)')
subplot(1,n,3)
plot(diff_baseline, diff_prefStim, 'k.')
xlabel('Baseline (large - small)')
ylabel('Pref stim resp (large - small)')
subplot(1,n,4)
plot(diff_modulation, diff_prefStim, 'k.')
xlabel('Tuning modulation (large - small)')
ylabel('Pref stim resp (large - small)')
subplot(1,n,5)
plot(diff_meanStim, diff_prefStim, 'k.')
xlabel('Mean stim resp (large - small)')
ylabel('Pref stim resp (large - small)')
subplot(1,n,6)
plot(diff_meanStim, diff_modulation, 'k.')
xlabel('Mean stim resp (large - small)')
ylabel('Tuning modulation (large - small)')

n=6;
figure('Position',[1675 20 2385 315])
subplot(1,n,1)
plot(diff_baseline./mean_baseline, diff_modulation./mean_modulation, 'k.')
xlabel('Baseline (large - small, norm.)')
ylabel('Tuning modulation (large - small, norm.)')
title('Differences, divided by means')
subplot(1,n,2)
plot(diff_baseline./mean_baseline, diff_meanStim./mean_meanStim, 'k.')
xlabel('Baseline (large - small, norm.)')
ylabel('Mean stim resp (large - small, norm.)')
subplot(1,n,3)
plot(diff_baseline./mean_baseline, diff_prefStim./mean_prefStim, 'k.')
xlabel('Baseline (large - small, norm.)')
ylabel('Pref stim resp (large - small, norm.)')
subplot(1,n,4)
plot(diff_modulation./mean_modulation, diff_prefStim./mean_prefStim, 'k.')
xlabel('Tuning modulation (large - small, norm.)')
ylabel('Pref stim resp (large - small, norm.)')
subplot(1,n,5)
plot(diff_meanStim./mean_meanStim, diff_prefStim./mean_prefStim, 'k.')
xlabel('Mean stim resp (large - small, norm.)')
ylabel('Pref stim resp (large - small, norm.)')
subplot(1,n,6)
plot(diff_meanStim./mean_meanStim, diff_modulation./mean_modulation, 'k.')
xlabel('Mean stim resp (large - small, norm.)')
ylabel('Tuning modulation (large - small, norm.)')

n=4;
bins = -1.95:.1:2;
figure('Position',[2 800 1915 315])
subplot(1,n,1)
hist(diff_baseline./mean_baseline, bins)
xlabel('Baseline (large - small, norm.)')
subplot(1,n,2)
hist(diff_meanStim./mean_meanStim, bins)
xlabel('Mean stim resp (large - small, norm.)')
subplot(1,n,3)
hist(diff_modulation./mean_modulation, bins)
xlabel('Tuning modulation (large - small, norm.)')
subplot(1,n,4)
hist(diff_prefStim./mean_prefStim, bins)
xlabel('Pref stim resp (large - small, norm.)')

n=4;
figure('Position',[2 410 1915 315])
subplot(1,n,1)
plot(diff_baseline, mean_baseline, 'k.')
xlabel('Baseline (large - small)')
ylabel('Baseline (mean)')
title('Differences vs. means')
subplot(1,n,2)
plot(diff_meanStim, mean_meanStim, 'k.')
xlabel('Mean stim resp (large - small)')
ylabel('Mean stim resp (mean)')
subplot(1,n,3)
plot(diff_modulation, mean_modulation, 'k.')
xlabel('Tuning modulation (large - small)')
ylabel('Tuning modulation (mean)')
subplot(1,n,4)
plot(diff_prefStim, mean_prefStim, 'k.')
xlabel('Pref stim resp (large - small)')
ylabel('Pref stim resp (mean)')

n=4;
figure('Position',[2 20 1915 315])
subplot(1,n,1)
plot(suppressed_max, diff_baseline./mean_baseline, 'k.')
xlabel('Suppressed by contrast (max-baseline)')
ylabel('Baseline (large - small, norm.)')
title('Features of suppressed-by-contrast cells')
subplot(1,n,2)
plot(suppressed_max, diff_modulation./mean_modulation, 'k.')
xlabel('Suppressed by contrast (max-baseline)')
ylabel('Tuning modulation (large - small, norm.)')
subplot(1,n,3)
plot(suppressed_max, diff_meanStim./mean_meanStim, 'k.')
xlabel('Suppressed by contrast (max-baseline)')
ylabel('Mean stim resp (large - small, norm.)')
subplot(1,n,4)
plot(suppressed_max, diff_prefStim./mean_prefStim, 'k.')
xlabel('Suppressed by contrast (max-baseline)')
ylabel('Pref stim resp (large - small, norm.)')

% input(1).baseline_small = baseline_small;
% input(1).baseline_large = baseline_large;
% input(1).meanStim_small = meanStim_small;
% input(1).meanStim_large = meanStim_large;
% input(1).maxStim_small = maxStim_small;
% input(1).maxStim_large = maxStim_large;
% input(1).modulation_small = modulation_small;
% input(1).modulation_large = modulation_large;
% input(1).baseline_ratio = baseline_ratio;
% input(1).prefStim_ratio = prefStim_ratio;
% input(1).suppressed_mean = suppressed_mean;
% input(1).suppressed_max = suppressed_max;
% ggobi(input);
        
%% Plot tuning curves (OLD, for SfN2016)
x = 0:360;
cols = 'kr';
lin = {'-','--'};
dirs = data1.directions(:,1);
condBlank = conditions(data1.blanks,:);
cols2 = 'kc';
for iCell = 1:size(responses,1)
    resp = squeeze(responses(iCell,:,:));
    br = squeeze(blankResponses(iCell,:,:));
    figure('position',[330 655 1530 420])
    subplot(1,3,1)
    hold on
    subplot(1,3,2)
    hold on
    h = zeros(1,4);
    maxi = 0;
    mini = 0;
    for c = 1:4
        col = mod(c,2);
        if col == 0
            col = 2;
        end
        col = cols(col);
        subplot(1,3,ceil(c/2))
        r = NaN(size(resp));
        r(condStim==c) = resp(condStim==c);
        m = nanmean(r,2);
        s = nanstd(r,0,2) ./ sqrt(sum(condStim==c,2));
        errorbar(dirs, m, s, ['o' col]);
        y = orituneWrapped(squeeze(parameters(iCell,c,:)), x);
        h(c) = plot(x,y,[col,lin{ceil(c/2)}],'LineWidth',2);
        if max(m+s) > maxi
            maxi = max(m+s);
        end
        if min(m-s) < mini
            mini = min(m-s);
        end
        r = NaN(size(br));
        r(condBlank==c) = br(condBlank==c);
        m = nanmean(r(:));
        s = nanstd(r(:)) / sqrt(sum(condBlank(:)==c));
        fill(x([1 end end 1]), [[1 1].*(m+s),[1 1].*(m-s)], 'k', ...
            'FaceColor',col,'EdgeColor','none','FaceAlpha',.3)
        plot(x([1 end]),[m m],[col lin{ceil(c/2)}],'LineWidth',1)
        if max(m+s) > maxi
            maxi = max(m+s);
        end
        if min(m-s) < mini
            mini = min(m-s);
        end
    end
    rng = maxi-mini;
    maxi = maxi + 0.05*rng;
    mini = mini - 0.05*rng;
    for k = 1:2
        subplot(1,3,k)
        xlim([-5 365])
        ylim([mini maxi])
        xlabel('Direction (degrees)')
        ylabel('Firing rate (Hz)')
    end
    subplot(1,3,1)
    legend(h(1:2), {'not running','running'})
    title('V1 intact')
    subplot(1,3,2)
    legend(h(3:4), {'not running','running'})
    title('V1 inactivated')
    
    % plot not running vs running
    subplot(1,3,3)
    hold on
    h = zeros(1,2);
    for c = 1:2
        c1 = c*2-1;
        c2 = c*2;
        r1 = NaN(size(resp));
        r1(condStim==c1) = resp(condStim==c1);
        r2 = NaN(size(resp));
        r2(condStim==c2) = resp(condStim==c2);
        ind = ~all(isnan([r1,r2]),2);
        resp1 = nanmean(r1(ind,:),2);
        resp2 = nanmean(r2(ind,:),2);
        maxi = max([resp1; resp2]);
        resp1 = resp1 ./ maxi;
        resp2 = resp2 ./ maxi;
        h(c) = plot(resp1, resp2, ['o' cols2(c)]);
        mini = min([0; resp1; resp2]);
        maxi = 1;
        rng = maxi-mini;
        mini = mini-.05*rng;
        maxi = maxi+.05*rng;
        if ~isnan(intercepts(iCell,c))
            plot([mini maxi],[mini maxi].*slopes(iCell,c)+intercepts(iCell,c), ...
                'Color',cols2(c),'LineWidth',2)
        end
    end
    plot([mini maxi],[mini maxi],'k:')
    axis([mini maxi mini maxi])
    axis square
    xlabel('not running')
    ylabel('running')
    legend(h,{['V1 intact (' modulations{iCell,1} ')'], ...
        ['V1 inactivated (' modulations{iCell,2} ')']})
    
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(fullfile('C:\RESULTS\Electrophys_SC\nonVisualEffects\plots_running\tuningCurves', ...
        sprintf('tuningCurve%03d.jpg', iCell)), '-djpeg','-r0')
    close gcf
end

%% OLD (fitting tuning curves for different running speeds)
% running = double(data1.running > runningThresh);
% running = ssLocal.getTracesPerStimulus(running, data1.stimMatrix, [1 0]); % [stimuli x repetitions]
% running = squeeze(mean(running, 4));
% running = double(running >= minDurPerTrial) + 1;
% responses = ssLocal.getTracesPerStimulus(data1.spikeRates, data1.stimMatrix, [0 0]);
% fs = 1/median(diff(data1.time));
% stimDur = mean(sum(data1.stimMatrix,2) ./ size(responses,3)) / fs;
% responses = sum(responses,4) ./ stimDur; % in Hz, [neuron x stim x rep]
% blankResponses = responses(:,data1.blanks,:);
% responses = responses(:,data1.directions(:,2),:);
% % conditions: 1: not running + no laser, 2: running + no laser, 3: not
% % running + laser on, 4: running + laser on
% conditions = running + data1.laserOn * 2;
% condStim = conditions(data1.directions(:,2),:);
% stim = repmat(data1.directions(:,1),1,size(responses,3));
% parameters = zeros(size(responses,1), 4, 5);
% intercepts = NaN(size(responses,1),2);
% slopes = NaN(size(responses,1),2);
% maxRespDiffs = NaN(size(responses,1),2);
% modulations = cell(size(responses,1),2);
% for iCell = 1:size(responses,1)
%     fprintf('%d ',iCell)
%     resp = squeeze(responses(iCell,:,:));
%     if sum(resp(:)>0) < 2
%         continue
%     end
%     
%     parsInit = fitoriWrapped(stim(:), resp(:));
%     pars = NaN(4, length(parsInit));
%     for cond = 1:4
%         if all(resp(condStim==cond) == 0)
%             continue
%         end
%         pars(cond,:) = fitoriWrapped(stim(condStim==cond), ...
%             resp(condStim==cond), [], [parsInit(1) NaN NaN NaN parsInit(5)]);
%     end
%     parameters(iCell,:,:) = pars;
%     
%     for c = 1:2
%         c1 = c*2-1;
%         c2 = c*2;
%         r1 = NaN(size(resp));
%         r1(condStim==c1) = resp(condStim==c1);
%         r2 = NaN(size(resp));
%         r2(condStim==c2) = resp(condStim==c2);
%         r = nansum(cat(3,r1,r2),3);
%         ind = ~all(isnan([r1,r2]),2);
%         s = stim(ind,:);
%         run = condStim(ind,:);
%         r = r(ind,:);
%         r1 = r1(ind,:);
%         r2 = r2(ind,:);
%         pVals = anovan(r(:),[s(:),run(:)],'model','interaction','display','off');
%         
%         resp1 = nanmean(r1,2);
%         resp2 = nanmean(r2,2);
%         maxi = max([resp1; resp2]);
%         resp1 = resp1 ./ maxi;
%         resp2 = resp2 ./ maxi;
%         if all(pVals >= 0.05) % neuron not modulated by stimulus or nonvisual signal
%             modulations{iCell,c} = 'none';
%         elseif pVals(2) < 0.05 && pVals(3) >= 0.05 % purely additive
%             [y,gof] = fit(resp1, resp2, @(a,x) x+a, ...
%                 'StartPoint', median(resp1)-median(resp2));
%             intercepts(iCell,c) = y.a;
%             slopes(iCell,c) = 1;
%             maxRespDiffs(iCell,c) = y.a;
%             modulations{iCell,c} = 'additive';
%         elseif pVals(2) >= 0.05 && pVals(3) < 0.05 % purely multiplicative
%             meanResp = nanmean([resp1; resp2]);
%             resp1 = resp1 - meanResp; % subtract mean response so that no intercept is necessary when fitting
%             resp2 = resp2 - meanResp;
%             [y,gof] = fit(resp1, resp2, @(a,x) x.*a, ...
%                 'StartPoint', nanmean(resp1./resp2));
%             if gof.adjrsquare < 0.3 % check whether the interaction is multiplicative
%                 modulations{iCell,c} = 'none';
%             else
%                 intercept = meanResp - y.a * meanResp; % shift fitted line by meanResp on x and y axis
%                 % calculate response difference;
%                 intercepts(iCell,c) = intercept;
%                 slopes(iCell,c) = y.a;
%                 % calculate response difference
%                 if abs(y.a) <= 1
%                     % if the absolute slope is smaller one, then the response
%                     % to the most driving stimulus during small pupil is 1
%                     respDiff = y.a+intercept - 1;
%                 else
%                     % if the absolute slope is larger one, then the response
%                     % to the most driving stimulus during large pupil is 1
%                     respDiff = 1 - (1-intercept)/y.a;
%                 end
%                 maxRespDiffs(iCell,c) = respDiff;
%                 modulations{iCell,c} = 'multiplicative';
%             end
%         elseif pVals(2) < 0.05 && pVals(3) < 0.05 % mixed additive and muliplicative
%             y = fitlm(resp1,resp2);
%             coefficients = y.Coefficients.Estimate;
%             intercepts(iCell,c) = coefficients(1);
%             slopes(iCell,c) = coefficients(2);
%             if abs(coefficients(2)) <= 1
%                 respDiff = coefficients(2)+coefficients(1) - 1;
%             else
%                 respDiff = 1 - (1-coefficients(1))/coefficients(2);
%             end
%             maxRespDiffs(iCell,c) = respDiff;
%             modulations{iCell,c} = 'mixed';
%         else % only modulated by stimulus
%             intercepts(iCell,c) = 0;
%             slopes(iCell,c) = 1;
%             maxRespDiffs(iCell,c) = 0;
%             modulations{iCell,c} = 'stimOnly';
%         end
%     end
% end
% fprintf('\n')

%% Plot max response differences between running vs not running (OLD, for SfN2016)
maxi = 1.6;
col = lines(1);
figure
plot(maxRespDiffs(:,1),maxRespDiffs(:,2),'o','color',col)
hold on
plot([-maxi maxi],[-maxi maxi],'k:')
plot([-maxi maxi],[0 0],'k:')
plot([0 0],[-maxi maxi],'k:')
ind=~any(isnan(maxRespDiffs),2);
p = signtest(maxRespDiffs(ind,1),maxRespDiffs(ind,2));
m = median(maxRespDiffs(ind,1)-maxRespDiffs(ind,2));
plot([-maxi maxi],[-maxi maxi]-m, 'color',col)
axis([-maxi maxi -maxi maxi])
axis square
set(gca,'box','off')
title(sprintf('Response diff. at pref. stim.: no running - running (n = %d, p = %.3f)',sum(ind), p))
xlabel('when V1 intact')
ylabel('when V1 inactivated')