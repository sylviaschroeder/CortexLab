folderROIData = 'C:\DATA\InfoStructs';
folderResults = 'C:\RESULTS\nonVisualEffects\modelGratingResp\';

%% Parameters
nonVisualSignal = 'pupil'; % 'running' or 'pupil'
minR2 = -Inf;
minEV = 0;

% smoothing (low-pass filter) of non-visual signal before
smoothing = 3; %in sec
filtPoly = 3;

doSave = 1;

data = load(fullfile(folderResults, 'kernelFit', 'results.mat'));
results = data.results;
nonVisPerTrial = struct('plane', cell(length(results),1));

%% Cross-validate several models to fit neural responses
% load results
for k=1:length(results)
    fprintf('Dataset %d: %s %s exp.: %d\n', k, results(k).subject, ...
        results(k).date, results(k).exp);
    folder = [folderROIData filesep results(k).subject filesep ...
        results(k).date filesep num2str(results(k).exp)];
    fileStart = [results(k).date '_' num2str(results(k).exp) '_' ...
        results(k).subject];
    file = [fileStart '_2P_plane%03d_ROI.mat'];
    for iPlane = 1:length(results(k).planes)
        display(['  Plane ' num2str(results(k).planes(iPlane))])
        data=load(fullfile(folder, sprintf(file,results(k).planes(iPlane))));
        meta = data.meta;
        
        if iPlane == 1
            nonVisData = [];
            % load ball or pupil data
            switch nonVisualSignal
                case 'running'
                    ballData = nonVis.getRunningSpeed(meta);
                    if ~isempty(ballData)
                        filtWindow = ceil(smoothing / median(diff(ballData.t)));
                        if mod(filtWindow,2) == 0
                            filtWindow = filtWindow-1;
                        end
                        total = ballData.total ./ median(diff(ballData.t)) ./ 53;
                        nonVisData = sgolayfilt(total, filtPoly, filtWindow);
                        nonVisTime = ballData.t;
                    end
                case 'pupil'
                    [pupilData, nonVisTime] = nonVis.loadPupilData(meta);
                    if ~isempty(pupilData)
                        nonVisTime(length(pupilData.x)+1:end) = [];
                        nonVisData = nonVis.getPupilDiam(pupilData);
                    end
            end
        end
        
        [stimTimes, stimSeq, stimMatrix, frameTimes, samplingRate] = ...
            ssLocal.getStimulusResponseInfo(meta);
        stimDuration = mean(stimTimes.offset - stimTimes.onset);
        [~,blankStim] = gratings.getOrientations(stimSeq);
        
        % get nonvisual data (average across 1 s before stim onset and stim
        % offset) (code changed 18.10.2016)
        offsets = [1 0];
        if ~isempty(nonVisData)
            nonVisInt = interp1(nonVisTime, nonVisData, frameTimes, 'pchip');
            nonVisual = ssLocal.getTracesPerStimulus(nonVisInt(:), ...
                stimMatrix, offsets);
            nonVisual(isnan(nonVisual)) = 0;
            nonVisual = permute(nonVisual, [4 3 2 1]);
        end
        nonVisTimeAvg = squeeze(mean(nonVisual,1));
        nonvisCentered = nonVisTimeAvg - mean(nonVisTimeAvg(:));
        nonvisCentered(:,blankStim) = [];
        nonVisPerTrial(k).plane(iPlane).centered = nonvisCentered;
        nonVisPerTrial(k).plane(iPlane).continuous = nonVisInt;
        nonVisPerTrial(k).plane(iPlane).trialsWithBlanks = nonVisTimeAvg;
        [results(k).plane(iPlane).crossVal(1:length( ...
            results(k).plane(iPlane).cellIDs)).explVars] = deal([]);
        for iCell = 1:length(results(k).plane(iPlane).cellIDs)
            nChars = fprintf('    Cell %d of %d', iCell, ...
                length(results(k).plane(iPlane).cellIDs));
            modelPars = results(k).plane(iPlane).modelPars(iCell);
            
            if isempty(modelPars.alphaEachTrial)
                fprintf(repmat('\b', 1, nChars));
                continue
            end
            
            % center variables (subtract mean)
            baselineCentered = modelPars.baselineEachTrial - ...
                mean(modelPars.baselineEachTrial(:));
            baselineCentered(:,blankStim) = [];
            
            % run crossvalidation
            resp = modelPars.alphaEachTrial;
            resp(:,blankStim) = [];
            
% (A) Crossvalidate models, only get back explained variance
            % (1) Use multiple linear models and models with interaction
            % between variables
            crossVals = models.runCrossValsSingleCell( ...
                resp, nonvisCentered, baselineCentered);          
            results(k).plane(iPlane).crossVal(iCell).explVars = ...
                crossVals.explVars;
            results(k).plane(iPlane).crossVal(iCell).formulas = ...
                crossVals.formulas;
            % (2) Use multiplicative model
%             error = models.crossvalidate(@models.testMultiplicativeModel, ...
%                 {resp}, {nonvisCentered});
%             explVars = 1 - sum(error{1}(:).^2) / ...
%                 sum((resp(:)-mean(resp(:))).^2);
%             results(k).plane(iPlane).crossVal(iCell).explVars = explVars;
%             results(k).plane(iPlane).crossVal(iCell).formulas = ...
%                 sprintf('y ~ 1 + stim*%s (through [x0 y0])',nonVisualSignal);

% (B) only fit one model and get parameters
%             stim = reshape(repmat(1:size(nonvisCentered,2), ...
%                 size(nonvisCentered,1), 1), [], 1);
%             nonvis = nonvisCentered(:);
%             base = baselineCentered(:);
%             resp = resp(:);
%             tbl = table(nonvis, base, stim, resp);
% %             mdl = fitlm(tbl, 'resp ~ 1 + stim + nonvis + base', ...
% %                 'CategoricalVars', 3);
%             mdl = fitlm(tbl, 'resp ~ 1 + stim + nonvis', ...
%                 'CategoricalVars', 3);
%             parameters = models.getModelParameters(mdl);
%             modelPars.parameters = parameters;
            
            
%             results(k).plane(iPlane).modelPars(iCell).parameters = ...
%                 modelPars.parameters;

% END (A) and (B)
            fprintf(repmat('\b', 1, nChars));
        end
        if doSave == 1
            if ~exist(fullfile(folderResults, nonVisualSignal), 'dir')
                mkdir(fullfile(folderResults, nonVisualSignal));
            end
%             save(fullfile(folderResults, nonVisualSignal, ...
%                 'linFit_nonvis_baseline.mat'), 'results');
%             save(fullfile(folderResults, nonVisualSignal, ...
%                 'linFit_nonvis_baseline_fullAdditiveModelOnly.mat'), 'results');
%             save(fullfile(folderResults, nonVisualSignal, ...
%                 'linFit_nonvis_baseline_stimNonvisModelOnly.mat'), 'results');
%             save(fullfile(folderResults, nonVisualSignal, ...
%                 'nonVis.mat'), 'nonVisPerTrial');
            save(fullfile(folderResults, nonVisualSignal, ...
                'multplic_nonvis.mat'), 'results');
        end
    end
end

%% Cross-validate models including running and pupil
labels = {'running','pupil'};
% load results
for k=1:length(results)
    fprintf('Dataset %d: %s %s exp.: %d\n', k, results(k).subject, ...
        results(k).date, results(k).exp);
    folder = [folderROIData filesep results(k).subject filesep ...
        results(k).date filesep num2str(results(k).exp)];
    fileStart = [results(k).date '_' num2str(results(k).exp) '_' ...
        results(k).subject];
    file = [fileStart '_2P_plane%03d_ROI.mat'];
    for iPlane = 1:length(results(k).planes)
        display(['  Plane ' num2str(results(k).planes(iPlane))])
        data=load(fullfile(folder, sprintf(file,results(k).planes(iPlane))));
        meta = data.meta;
        
        if iPlane == 1
            % load ball and pupil data
            nonVisData = cell(1,2);
            nonVisTime = cell(1,2);
            ballData = nonVis.getRunningSpeed(meta);
            if ~isempty(ballData)
                filtWindow = ceil(smoothing / median(diff(ballData.t)));
                if mod(filtWindow,2) == 0
                    filtWindow = filtWindow-1;
                end
                total = ballData.total ./ median(diff(ballData.t)) ./ 53;
                nonVisData{1} = sgolayfilt(total, filtPoly, filtWindow);
                nonVisTime{1} = ballData.t;
            end
            [pupilData, nonVisTime{2}] = nonVis.loadPupilData(meta);
            if ~isempty(pupilData)
                nonVisTime{2}(length(pupilData.x)+1:end) = [];
                nonVisData{2} = nonVis.getPupilDiam(pupilData);
            end
        end
        
        [stimTimes, stimSeq, stimMatrix, frameTimes, samplingRate] = ...
            ssLocal.getStimulusResponseInfo(meta);
        stimDuration = mean(stimTimes.offset - stimTimes.onset);
        [~,blankStim] = gratings.getOrientations(stimSeq);
        
        offsets = [1 0];
        nonvisZscored = cell(1,2);
        for nv = 1:2
            nonVisInt = interp1(nonVisTime{nv}, nonVisData{nv}, frameTimes, 'pchip');
            nonVisual = ssLocal.getTracesPerStimulus(nonVisInt(:), ...
                stimMatrix, offsets);
%             nonVisual(isnan(nonVisual)) = 0;
            nonVisual = permute(nonVisual, [4 3 2 1]);
            nonVisTimeAvg = squeeze(nanmean(nonVisual,1));
            nonvisZscored{nv} = (nonVisTimeAvg - mean(nonVisTimeAvg(:))) ./ ...
                std(nonVisTimeAvg(:));
            nonvisZscored{nv}(:,blankStim) = [];
            nonVisPerTrial(k).plane(iPlane).(labels{nv}).zscored = nonvisZscored{nv};
            nonVisPerTrial(k).plane(iPlane).(labels{nv}).continuous = nonVisInt;
            nonVisPerTrial(k).plane(iPlane).(labels{nv}).trialsWithBlanks = nonVisTimeAvg;
        end
        [results(k).plane(iPlane).crossVal(1:length( ...
            results(k).plane(iPlane).cellIDs)).explVars] = deal([]);
        for iCell = 1:length(results(k).plane(iPlane).cellIDs)
            nChars = fprintf('    Cell %d of %d', iCell, ...
                length(results(k).plane(iPlane).cellIDs));
            modelPars = results(k).plane(iPlane).modelPars(iCell);
            
            if isempty(modelPars.alphaEachTrial)
                fprintf(repmat('\b', 1, nChars));
                continue
            end
            
            % run crossvalidation
            resp = modelPars.alphaEachTrial;
            resp(:,blankStim) = [];
            
% (A) fit all possible models
            crossVals = models.runCrossValsSingleCell( ...
                resp, nonvisZscored{1}, nonvisZscored{2});
            f = strrep(crossVals.formulas,'nonvis','running');
            f = strrep(f,'base','pupil');
            crossVals.formulas = f;
            
            % get best fitting parameters
%             [explVar, best] = max(crossVals.explVars);
%             if explVar <= minEV
%                 modelPars.stim = [];
%                 modelPars.nonvis = [];
%                 modelPars.baseline = [];
%                 modelPars.interact = [];
%                 modelPars.parameters = [];
%             else
%                 stim = reshape(repmat(1:size(nonvisCentered,2), ...
%                     size(nonvisCentered,1), 1), [], 1);
%                 nonvis = nonvisCentered(:);
%                 base = baselineCentered(:);
%                 resp = resp(:);
%                 tbl = table(nonvis, base, stim, resp);
%                 mdl = fitlm(tbl, crossVals.formulas{best}, 'CategoricalVars', 3);
%                 parameters = models.getModelParameters(mdl);
%                 modelPars.parameters = parameters;
%             end
            results(k).plane(iPlane).crossVal(iCell).explVars = ...
                crossVals.explVars;
            results(k).plane(iPlane).crossVal(iCell).formulas = ...
                crossVals.formulas;
% (B) only fit one specific model and get parameters
%             stim = reshape(repmat(1:size(resp,2), ...
%                 size(resp,1), 1), [], 1);
%             running = nonvisZscored{1}(:);
%             pupil = nonvisZscored{2}(:);
%             resp = resp(:);
%             tbl = table(running, pupil, stim, resp);
% %             mdl = fitlm(tbl, 'resp ~ 1 + stim*(running+pupil)', ...
% %                 'CategoricalVars', 3);
%             mdl = fitlm(tbl, 'resp ~ 1 + stim+running+pupil', ...
%                 'CategoricalVars', 3);
%             parameters = models.getModelParameters(mdl);
%             modelPars.parameters = parameters;
%             results(k).plane(iPlane).modelPars(iCell).parameters = ...
%                 modelPars.parameters;
% END (A) and (B) ---------------------------------------------------------
            fprintf(repmat('\b', 1, nChars));
        end
        if doSave == 1
            if ~exist(fullfile(folderResults, nonVisualSignal), 'dir')
                mkdir(fullfile(folderResults, nonVisualSignal));
            end
%             save(fullfile(folderResults, nonVisualSignal, ...
%                 'linFit_nonvis_baseline.mat'), 'results');
%             save(fullfile(folderResults, nonVisualSignal, ...
%                 'linFit_nonvis_baseline_fullAdditiveModelOnly.mat'), 'results');
%             save(fullfile(folderResults, nonVisualSignal, ...
%                 'linFit_nonvis_baseline_stimNonvisModelOnly.mat'), 'results');
%             save(fullfile(folderResults, ...
%                 'linFit_running_pupil.mat'), 'results');
            save(fullfile(folderResults, ...
                'linFit_running_pupil_interactiveOnly_parameters.mat'), 'results');
%             save(fullfile(folderResults, ...
%                 'nonVis_runningAndPupil.mat'), 'nonVisPerTrial');
        end
    end
end

%% Check correlation between baseline and nonvisual signal (continuous signals)
for iExp = 1:length(results)
    fprintf('Dataset %d: %s %s exp.: %d\n', iExp, results(iExp).subject, ...
        results(iExp).date, results(iExp).exp);
    folder = [folderROIData filesep results(iExp).subject filesep ...
        results(iExp).date filesep num2str(results(iExp).exp)];
    fileStart = [results(iExp).date '_' num2str(results(iExp).exp) '_' ...
        results(iExp).subject];
    file = [fileStart '_2P_plane%03d_ROI.mat'];
    for iPlane = 1:length(results(iExp).plane)
        % load meta
        data=load(fullfile(folder, sprintf(file,results(iExp).planes(iPlane))));
        meta = data.meta;
        frameTimes = ppbox.getFrameTimes(meta);
        if iPlane == 1
            nonVisData = [];
            % load ball or pupil data
            switch nonVisualSignal
                case 'running'
                    ballData = nonVis.getRunningSpeed(meta);
                    if ~isempty(ballData)
                        filtWindow = ceil(smoothing / median(diff(ballData.t)));
                        if mod(filtWindow,2) == 0
                            filtWindow = filtWindow-1;
                        end
                        nonVisData = sgolayfilt(ballData.total, filtPoly, filtWindow);
                        nonVisTime = ballData.t;
                    end
                case 'pupil'
                    [pupilData, nonVisTime] = nonVis.loadPupilData(meta);
                    if ~isempty(pupilData)
                        nonVisTime(length(pupilData.x)+1:end) = [];
                        nonVisData = nonVis.getPupilDiam(pupilData);
                    end
            end
        end
        nonVisInt = interp1(nonVisTime, nonVisData, frameTimes, 'pchip')';
        nonvisZScored = (nonVisInt-nanmean(nonVisInt))./nanstd(nonVisInt);
        fPlots = fullfile(folderResults, ['plots_' nonVisualSignal], ...
            'corrWithFittedBaseline', ...
            [results(iExp).subject '_' results(iExp).date '_' ...
            num2str(results(iExp).exp)], num2str(results(iExp).planes(iPlane)));
        if ~exist(fPlots, 'dir')
            mkdir(fPlots);
        end
        for iCell = 1:length(results(iExp).plane(iPlane).cellIDs)
            baseline = results(iExp).plane(iPlane).modelPars(iCell).baseline;
            
            if isempty(baseline)
                j = results(iExp).plane(iPlane).cellIDs(iCell);
                baseline = (meta.Fcorr(:,j) - meta.F0(:,j)) ./ ...
                    max(1, mean(meta.F0(:,j)));
                baselineZScored = (baseline-mean(baseline))./std(baseline);
                labels = {'Nonvisual signal','Activity','F0'};
            else
                baselineZScored = (baseline-mean(baseline))./std(baseline);
                labels = {'Nonvisual signal','Baseline','F0'};
            end
            F0 = meta.F0(:,results(iExp).plane(iPlane).cellIDs(iCell));
            F0 = (F0-mean(F0))./std(F0);
            figure('position',[2 678 1918 420])
            subplot(1,2,1)
            plot(frameTimes,[nonvisZScored, baselineZScored-4, F0-8])
            legend(labels)
            axis tight
            xlabel('Time (in s)')
            title(sprintf('Neuron %d, %s %s exp %d, plane %d', ...
                results(iExp).plane(iPlane).cellIDs(iCell), ...
                results(iExp).subject, results(iExp).date, ...
                results(iExp).exp, results(iExp).planes(iPlane)), ...
                'interpreter','none')
            
            subplot(1,4,3)
            plot(nonVisInt, baseline, 'k.')
            xlabel('Nonvisual signal')
            ylabel('Baseline')
            subplot(1,4,4)
            [r,lags] = xcorr(nonvisZScored,baselineZScored,2000,'unbiased');
            l = zeros(size(lags));
            ind = abs(lags)>0;
            l(ind) = frameTimes(abs(lags(ind))) .* sign(lags(ind));
            plot(l, r)
            hold on
            maxi=max(r);
            mini=min(r);
            rn = maxi-mini;
            maxi=maxi+.1*rn;
            mini=mini-.1*rn;
            plot(l([1 end]),[0 0], 'k:')
            plot([0 0],[mini maxi], 'k:')
            axis tight
            xlabel('Lag (baseline+x vs. nonvis)')
            ylabel('Corr.')
            
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(fullfile(fPlots, sprintf('baselineCorr%03d.jpg', iCell)), ...
                '-djpeg','-r0')
%             saveas(gcf, fullfile(fPlots, sprintf('baselineCorr%03d.jpg', iCell)));
            close gcf
        end
    end
end

%% Compare correlation of baseline and nonvisual signal to their linear coefficients in prediction neural response
data = load(fullfile(folderResults, nonVisualSignal, 'nonVis.mat'));
nonVisPerTrial = data.nonVisPerTrial;
data = load(fullfile(folderResults, nonVisualSignal, ...
    'linFit_nonvis_baseline_fullAdditiveModelOnly.mat'));
results = data.results;
rhosBaseline = [];
pValuesBaseline = [];
rhosResponse = [];
pValuesResponse = [];
nonVisual = [];
baseline = [];
for iExp = 1:length(results)
    fprintf('Dataset %d: %s %s exp.: %d\n', iExp, results(iExp).subject, ...
        results(iExp).date, results(iExp).exp);
    folder = [folderROIData filesep results(iExp).subject filesep ...
        results(iExp).date filesep num2str(results(iExp).exp)];
    fileStart = [results(iExp).date '_' num2str(results(iExp).exp) '_' ...
        results(iExp).subject];
    file = [fileStart '_2P_plane%03d_ROI.mat'];
    for iPlane = 1:length(results(iExp).plane)
        data=load(fullfile(folder, sprintf(file,results(iExp).planes(iPlane))));
        meta = data.meta;
        stimSeq = ppbox.getStimSequence(meta);
        [~,blankStim] = gratings.getOrientations(stimSeq);
        nv = nonVisPerTrial(iExp).plane(iPlane).centered;
        stimEnd = round(results(iExp).plane(iPlane).stimDuration / ...
            median(diff(results(iExp).plane(iPlane).kernelTime)));
        for iCell = 1:length(results(iExp).plane(iPlane).cellIDs)
            if isempty(results(iExp).plane(iPlane).modelPars(iCell).baselineEachTrial)
                continue
            end
            kernelSign = sign(sum(results(iExp).plane(iPlane)...
                .modelPars(iCell).kernel(1:stimEnd)));
            
            bl = results(iExp).plane(iPlane).modelPars(iCell).baselineEachTrial;
            bl(:,blankStim) = [];
            [r,p] = corr(bl(:),nv(:));
            rhosBaseline = [rhosBaseline;r];
            pValuesBaseline = [pValuesBaseline;p];
            resp = results(iExp).plane(iPlane).modelPars(iCell).alphaEachTrial;
            resp(:,blankStim) = [];
            resp = resp * kernelSign;
            [r,p] = corr(resp(:),nv(:));
            rhosResponse = [rhosResponse; r];
            pValuesResponse = [pValuesResponse; p];
            nonVisual = [nonVisual; results(iExp).plane(iPlane) ...
                .modelPars(iCell).parameters.nonvisTotal] * kernelSign;
            baseline = [baseline; results(iExp).plane(iPlane) ...
                .modelPars(iCell).parameters.baseTotal] * kernelSign;
        end
    end
end

figure
hist(rhosBaseline,-.5:.1:.5)
xlabel('Corr. of baseline and non visual')
ylabel('#Neurons')
figure
hist(rhosResponse,-.5:.1:.5)
xlabel('Corr. of visual response and non visual')
ylabel('#Neurons')
figure
scatter(rhosBaseline, rhosResponse)
xlabel('Corr. of baseline and non visual')
ylabel('Corr. of visual response and non visual')

figure
scatter(rhosBaseline, nonVisual)
xlabel('Corr. of baseline and non visual')
ylabel('Coeff. of non visual in linear model')
figure
scatter(rhosBaseline, baseline)
xlabel('Corr. of baseline and non visual')
ylabel('Coeff. of baseline in linear model')


%% Plot which models predict the cell responses best
data = load(fullfile(folderResults, ...
    'linFit_running_pupil_additiveOnly_parameters.mat'));
results = data.results;
data = load(fullfile(folderResults, 'running', 'multplic_nonvis.mat'));
resultsMultRunning = data.results;
data = load(fullfile(folderResults, 'pupil', 'multplic_nonvis.mat'));
resultsMultPupil = data.results;
influences = [];
mods = [results(1).plane(1).crossVal(1).formulas, ...
    resultsMultPupil(1).plane(1).crossVal(1).formulas, ...
    resultsMultRunning(1).plane(1).crossVal(1).formulas, ...
    {'none', 'no stim resp'}];
gadGroups = [-1 0 1];
explainedVariance = [];
gadAll = [];
for iExp = 1:length(results)
    fprintf('Dataset %d: %s %s exp.: %d\n', iExp, results(iExp).subject, ...
        results(iExp).date, results(iExp).exp);
    folder = [folderROIData filesep results(iExp).subject filesep ...
        results(iExp).date filesep num2str(results(iExp).exp)];
    fileStart = [results(iExp).date '_' num2str(results(iExp).exp) '_' ...
        results(iExp).subject];
    file = [fileStart '_2P_plane%03d_ROI.mat'];
    inflPl = zeros(length(mods),3);
    for iPlane = 1:length(results(iExp).plane)
        % load meta
        info=load(fullfile(folder, sprintf(file,results(iExp).planes(iPlane))));
        meta = info.meta;
        gad = meta.ROI.isGad;
        isGad = gad(results(iExp).plane(iPlane).cellIDs);
        gadAll = [gadAll; gad];
        
        allExplVars = NaN(length(results(iExp).plane(iPlane).crossVal), ...
            length(results(1).plane(1).crossVal(1).formulas));
        explVarMult = NaN(length(results(iExp).plane(iPlane).crossVal),2);
        for iCell = 1:length(results(iExp).plane(iPlane).crossVal)
            ev = results(iExp).plane(iPlane).crossVal(iCell).explVars;
            if ~isempty(ev)
                allExplVars(iCell,:) = ev;
                explVarMult(iCell,1) = resultsMultPupil(iExp) ...
                    .plane(iPlane).crossVal(iCell).explVars;
                explVarMult(iCell,2) = resultsMultRunning(iExp) ...
                    .plane(iPlane).crossVal(iCell).explVars;
            end
        end
        explainedVariance = [explainedVariance; [allExplVars, explVarMult]];
        [maxExplVars, inds] = max([allExplVars, explVarMult], [], 2);
        % detect neurons for which there was no stimulus response
        % (fitting kernel+baseline was worse than just fitting baseline)
        noStimInds = isnan(maxExplVars);
        inds(noStimInds) = 0;
        infl = zeros(length(mods),3);
        for g = 1:3
            for m = 1:length(results(1).plane(1).crossVal(1).formulas)+1
                infl(m,g) = sum(inds==m & maxExplVars>=minEV & isGad==gadGroups(g));
            end
            infl(end-1,g) = sum((inds>0 & maxExplVars<minEV) & ...
                isGad==gadGroups(g));
            infl(end,g) = sum(inds==0 & isGad==gadGroups(g));
        end
        inflPl = inflPl + infl;
    end
    influences = cat(3, influences, inflPl);
end
reorder = [1:3 6 5 14 13 4 8 7 9 12 10 11 15:16];
influences = influences(reorder,:,:);
mods = mods(reorder);
mods = strrep(mods,'resp','y');
mods{end-1} = 'expl. var. < 0';
mods{end} = 'no stim response';
reorder = [1:3 6 5 14 13 4 8 7 9 12 10 11];
explainedVariance = explainedVariance(:,reorder);


bars = sum(influences, 3);
figure
barh(bsxfun(@rdivide, bars, sum(bars,1)) .* 100)
set(gca, 'YTick', 1:size(bars,1), 'YTickLabel', mods, 'YDir','reverse')
xlabel('Neurons (%)')
title('Effect of stimulus, running and pupil')
figure
barh(sum(bars,2))
set(gca, 'YTick', 1:size(bars,1), 'YTickLabel', mods, 'YDir','reverse', 'box','off')
xlabel('# Neurons')
title('Best model for each neuron')
% plot numbers of exc, NA, inh neurons
figure
bar(sum(bars,1))
set(gca, 'XTick', 1:3, 'XTickLabel', {'exc', 'NA', 'inh'})
ylabel('# Neurons')
title('Cell types')
bars_old = bars;

% without 'bad kernel fit' model
bars2 = sum(influences(1:end-1,:,:), 3);
figure
bar(bsxfun(@rdivide, bars2, sum(bars2,1)) .* 100)
set(gca, 'XTick', 1:size(bars2,1), 'XTickLabel', mods(1:end-1), ...
    'XTickLabelRotation', 90)
ylabel('Neurons (%)')
title(sprintf('Effect of stimulus, %s and baseline', nonVisualSignal))
figure
bar(sum(bars2,2))
set(gca, 'XTick', 1:size(bars2,1), 'XTickLabel', mods(1:end-1), ...
    'XTickLabelRotation', 90)
ylabel('# Neurons')
title(sprintf('Effect of stimulus, %s and baseline', nonVisualSignal))
% plot numbers of exc, NA, inh neurons
figure
bar(sum(bars2,1))
set(gca, 'XTick', 1:3, 'XTickLabel', {'exc', 'NA', 'inh'})
ylabel('# Neurons')
title('Cell types')


% Compare R^2 values across models
R2 = cell(length(results),1);
% diffsFullModel = cell(length(results),1);
for iExp = 1:length(results)
    for iPlane = 1:length(results(iExp).plane)
        allExplVars = NaN(length(results(iExp).plane(iPlane).crossVal), ...
            length(results(1).plane(1).crossVal(1).formulas));
        explVarMult = NaN(length(results(iExp).plane(iPlane).crossVal),2);
        for iCell = 1:length(results(iExp).plane(iPlane).crossVal)
            ev = results(iExp).plane(iPlane).crossVal(iCell).explVars;
            if ~isempty(ev)
                allExplVars(iCell,:) = ev;
                explVarMult(iCell,1) = resultsMultPupil(iExp) ...
                    .plane(iPlane).crossVal(iCell).explVars;
                explVarMult(iCell,2) = resultsMultRunning(iExp) ...
                    .plane(iPlane).crossVal(iCell).explVars;
            end
        end
        R2{iExp,iPlane} = [allExplVars, explVarMult];
%         d = [];
%         for k = setdiff(1:12,4)
%             d = [d, diff(allExplVars(:,[k 4]),1,2)];
%         end
%         diffsFullModel{iExp,iPlane} = d;
    end
end
allR2 = cat(1,R2{:});
allR2 = allR2(:,reorder);
ind = all(allR2 < minEV,2);
allR2pos = allR2;
allR2pos(ind,:) = NaN;
figure
% boxplot(allR2)
boxplot(allR2pos,'orientation','horizontal')
set(gca, 'YTick', 1:size(allR2,2), 'YTickLabel', mods, 'YDir','reverse', ...
    'XGrid','on')
xlabel('R^2','interpreter','tex')
xlim([-1 1])
title('R^2 for all models and all neurons')

% check difference of R2 between overall best model
% (1+stim+nonvis+baseline) to other models for each neuron
alldiffs = [];
for k = setdiff(1:12,4)
    alldiffs = [alldiffs, diff(allR2(:,[k 4]),1,2)];
end
figure,boxplot(alldiffs)
ylim([-1 1.5])
title('Diff. in R^2 (resp~1+stim+nonvis+basel minus ...)')
ylabel('R^2')
set(gca, 'XTick', 1:size(alldiffs,2), 'XTickLabel', mods(setdiff(1:12,4)), ...
    'XTickLabelRotation', 90)

% scatter plots for various pairs of models
pairs = [12 9; 12 10; 10 8; 9 8; 9 6; 10 7; 8 6];
explVarsTuned = explainedVariance;
explVarsTuned(all(explVarsTuned<minEV,2),:) = [];
for p = 1:size(pairs,1)
    figure
    scatter(explVarsTuned(:,pairs(p,1)),explVarsTuned(:,pairs(p,2)))
    mini = min(reshape(explVarsTuned(:,pairs(p,:)),[],1));
    maxi = max(reshape(explVarsTuned(:,pairs(p,:)),[],1));
    rng = maxi-mini;
    mini = mini-.05*rng;
    maxi = maxi+.05*rng;
    hold on
    plot([mini maxi],[mini maxi],'r')
    axis([mini maxi mini maxi])
    axis square
    xlabel(mods{pairs(p,1)})
    ylabel(mods{pairs(p,2)})
end

% plot parameters of neurons where best model is much better than full
% additive model
minDiff = -.1;
for iExp = 1:length(results)
    for iPlane = 1:length(results(iExp).plane)
        neurons = find(any(diffsFullModel{iExp,iPlane} < minDiff,2));
        for iCell = neurons'
            pars = results(iExp).plane(iPlane).modelPars(iCell).parameters;
            if isempty(pars) || length(pars.stimuli) <= 1
                continue
            end
            f = [];
            if ~isempty(pars.nonvisEachStim)
                f(end+1) = figure('Position',[781 688 560 420]);
                scatter(pars.nonvisEachStim,pars.stimuli)
                xlabel('Nonvisual parameters')
                ylabel('Stimulus parameters')
                title(sprintf('Neuron %d, %s %s exp %d, plane %d', ...
                    iCell, results(iExp).subject, results(iExp).date, ...
                    results(iExp).exp, results(iExp).planes(iPlane)), ...
                    'interpreter','none')
            end
            if ~isempty(pars.baseEachStim)
                f(end+1) = figure('Position',[1351 688 560 420]);
                scatter(pars.baseEachStim,pars.stimuli)
                xlabel('Baseline parameters')
                ylabel('Stimulus parameters')
                title(sprintf('Neuron %d, %s %s exp %d, plane %d', ...
                    iCell, results(iExp).subject, results(iExp).date, ...
                    results(iExp).exp, results(iExp).planes(iPlane)), ...
                    'interpreter','none')
            end
            pause
            close(f)
        end
    end
end

%% Look at parameters of models
% load appropriate results structure first!
nonvisTotal = {};
nonvisEach = {};
stimEach = {};
baselineTotal = {};
baselineEach = {};
gadAll = [];
responseSignAll = [];
kernelSignAll = [];
for iExp = 1:length(results)
    folder = [folderROIData filesep results(iExp).subject filesep ...
        results(iExp).date filesep num2str(results(iExp).exp)];
    fileStart = [results(iExp).date '_' num2str(results(iExp).exp) '_' ...
        results(iExp).subject];
    file = [fileStart '_2P_plane%03d_ROI.mat'];
    for iPlane = 1:length(results(iExp).plane)
        % load meta
        info=load(fullfile(folder, sprintf(file,results(iExp).planes(iPlane))));
        meta = info.meta;
        stimEnd = round(results(iExp).plane(iPlane).stimDuration / ...
            median(diff(results(iExp).plane(iPlane).kernelTime)));
        for iCell = 1:length(results(iExp).plane(iPlane).cellIDs)
            if isempty(results(iExp).plane(iPlane).modelPars(iCell).kernel)
                continue
            end
%             maxEV = max(results(iExp).plane(iPlane).crossVal(iCell).explVars);
%             if maxEV < minEV
%                 continue
%             end
            parameters = results(iExp).plane(iPlane).modelPars(iCell).parameters;
            if isempty(parameters)
                continue
            end
            stimEach{end+1,1} = parameters.stimuli;
%             nonvisTotal{end+1,1} = parameters.nonvisTotal;
%             nonvisEach{end+1,1} = parameters.nonvisEachStim;
%             baselineTotal{end+1,1} = parameters.baseTotal;
%             baselineEach{end+1,1} = parameters.baseEachStim;
            nonvisTotal{end+1,1} = parameters.runningTotal;
            nonvisEach{end+1,1} = parameters.runningEachStim;
            baselineTotal{end+1,1} = parameters.pupilTotal;
            baselineEach{end+1,1} = parameters.pupilEachStim;
            
            gadAll(end+1) = meta.ROI.isGad(results(iExp).plane(iPlane).cellIDs(iCell));
            
            kernelSign = sign(sum(results(iExp).plane(iPlane)...
                .modelPars(iCell).kernel(1:stimEnd)));
            kernelSignAll(end+1) = kernelSign;
            medianResp = median(results(iExp).plane(iPlane) ...
                .modelPars(iCell).alphaEachTrial, 1);
            [~,maxInd] = max(abs(medianResp));
            responseSignAll(end+1) = sign(medianResp(maxInd)) * kernelSign;
        end
    end
end
indStim = cellfun(@length, stimEach) == 1;
% absStim = cellfun(@abs, stimEach, 'UniformOutput', false);
% [maxResp,maxStim] = cellfun(@max, absStim);
% stimSigns = NaN(size(maxStim));
% for k=1:length(maxStim)
%     stimSigns(k) = sign(stimEach{k}(maxStim(k)));
% end
% stimSigns(indStim) = NaN;
ind = cellfun(@isempty,nonvisTotal);
[nonvisTotal{ind}] = deal(NaN);
nonvisTotal = cell2mat(nonvisTotal);
ind = cellfun(@isempty,baselineTotal);
[baselineTotal{ind}] = deal(NaN);
baselineTotal = cell2mat(baselineTotal);

cols = lines(4);
% Magnitude of parameters
gadGroups = [-1 1];
respGroups = [-1 1];
% force sign of kernel during stimulus to be positive
ind = kernelSignAll==-1;
% Nonvisual parameter
nv = nonvisTotal;
nv(ind) = -nv(ind);
% x = -.3:.01:.3;
x = -1.5:.1:1.5;
figure
hist(nv,x)
% hist(nonvisTotal,-1:.01:1)
xlim(x([1 end]))
% title('Nonvisual influence')
title('Running (total) parameters')
ylabel('# Neurons')
figure
hold on
bars = zeros(length(x),2);
for k=1:2
    bars(:,k) = hist(nv(gadAll==gadGroups(k)),x);
end
bars = bsxfun(@rdivide,bars,sum(bars,1));
for k = 1:2
    plot(x,bars(:,k), 'LineWidth',2,'Color',cols(k,:))
end
legend('exc','inh')
xlim(x([1 end]))
% title('Nonvisual influence')
title('Running (total) parameters')
ylabel('Proportion of neurons (per type)')
figure
hold on
bars = zeros(length(x),2);
for k=1:2
    bars(:,k) = hist(nv(responseSignAll==respGroups(k)),x);
end
bars = bsxfun(@rdivide,bars,sum(bars,1));
for k = 1:2
    plot(x,bars(:,k), 'LineWidth',2,'Color',cols(k+2,:))
end
legend('suppressed','excited')
xlim(x([1 end]))
% title('Nonvisual influence')
title('Running (total) parameters')
ylabel('Proportion of neurons (per type)')
% figure
% hist(nonvisTotal.*stimSigns,-1:.05:1)
% xlim([-1 1])
% title('Nonvisual influence (times sign of stimulus response)')
% ylabel('# Neurons')

% Baseline parameter
bl = baselineTotal;
bl(ind) = -bl(ind);
figure
x = -1.5:.1:1.5;
hist(bl,x)
xlim(x([1 end]))
% title('Baseline influence')
title('Pupil (total) parameters')
ylabel('# Neurons')
figure
hold on
bars = zeros(length(x),2);
for k=1:2
    bars(:,k) = hist(bl(gadAll==gadGroups(k)),x);
end
bars = bsxfun(@rdivide,bars,sum(bars,1));
for k = 1:2
    plot(x,bars(:,k), 'LineWidth',2,'Color',cols(k,:))
end
legend('exc','inh')
xlim(x([1 end]))
% title('Baseline influence')
title('Pupil (total) parameters')
ylabel('Proportion of neurons (per type)')
figure
hold on
bars = zeros(length(x),2);
for k=1:2
    bars(:,k) = hist(bl(responseSignAll==respGroups(k)),x);
end
bars = bsxfun(@rdivide,bars,sum(bars,1));
for k = 1:2
    plot(x,bars(:,k), 'LineWidth',2,'Color',cols(k+2,:))
end
legend('suppressed','excited')
xlim(x([1 end]))
% title('Baseline influence')
title('Pupil (total) parameters')
ylabel('Proportion of neurons (per type)')
% figure
% hist(baselineTotal.*stimSigns,-3:.25:3)
% xlim([-3 3])
% title('Baseline influence (times sign of stimulus response)')
% ylabel('# Neurons')
% Stimulus parameter
x = 0:2.5:40;
se = [stimEach{indStim}];
if ~isempty(se)
    se(ind) = -se(ind);
    figure
    hist(se,x)
    xlim(x([1 end]))
    title('Stimulus impact in untuned neurons')
    ylabel('# Neurons')
    figure
    bars = zeros(length(x),2);
    for k=1:2
        bars(:,k) = hist(se(gadAll==gadGroups(k)),x);
    end
    bars = bsxfun(@rdivide,bars,sum(bars,1));
    plot(x,bars, 'LineWidth',2)
    legend('exc','inh')
    xlim(x([1 end]))
    title('Stimulus influence (untuned neurons)')
    ylabel('Proportion of neurons (per type)')
    figure
    bars = zeros(length(x),2);
    for k=1:2
        bars(:,k) = hist(se(responseSignAll==respGroups(k)),x);
    end
    bars = bsxfun(@rdivide,bars,sum(bars,1));
    plot(x,bars, 'LineWidth',2)
    legend('suppressed','excited')
    xlim(x([1 end]))
    title('Stimulus influence (untuned neurons)')
    ylabel('Proportion of neurons (per type)')
end
x = -10:20;
se = cat(1,stimEach{~indStim});
if ~isempty(se)
    se(ind,:) = -se(ind,:);
    figure
    hist(se(:),x)
    xlim(x([1 end]))
    title('Stimulus impact in tuned neurons')
    ylabel('# Neurons')
    figure
    hold on
    bars = zeros(length(x),2);
    for k=1:2
        bars(:,k) = hist(reshape(se(gadAll==gadGroups(k),:),[],1),x);
    end
    bars = bsxfun(@rdivide,bars,sum(bars,1));
    for k = 1:2
        plot(x,bars(:,k), 'LineWidth',2,'Color',cols(k,:))
    end
    legend('exc','inh')
    xlim(x([1 end]))
    title('Stimulus influence')
    ylabel('Proportion of neurons (per type)')
    figure
    hold on
    bars = zeros(length(x),2);
    for k=1:2
        bars(:,k) = hist(reshape(se(responseSignAll==respGroups(k),:),[],1),x);
    end
    bars = bsxfun(@rdivide,bars,sum(bars,1));
    for k = 1:2
        plot(x,bars(:,k), 'LineWidth',2,'Color',cols(k+2,:))
    end
    legend('suppressed','excited')
    xlim(x([1 end]))
    title('Stimulus influence')
    ylabel('Proportion of neurons (per type)')
end

% extremeStim = maxResp.*stimSigns;
% figure
% scatter(extremeStim, baselineTotal)
% xlabel('Max. abs. stimulus response (signed)')
% ylabel('Baseline parameter')
% figure
% scatter(extremeStim, nonvisTotal)
% xlabel('Max. abs. stimulus response (signed)')
% ylabel('Nonvisual parameter')
figure
% scatter(baselineTotal, nonvisTotal)
scatter(bl, nv)
xlabel('Pupil parameter')
ylabel('Running parameter')

%% OLD
% parameters per dataset across all models
% (1) stimulus (intercepts)
for iExp = 1:length(results)
    figure
    x = -5:.5:15;
    hist(stimPars{iExp}(:), x);
    xlim(x([1 end]))
    xlabel('Intercepts (stimulus influence)')
    ylabel('#Neurons*stimuli')
    title(sprintf('Dataset %d: intercepts (each stimulus)', iExp))
end
% (2) Nonvisual signal
for iExp = 1:length(results)
    figure
    x = -1:.05:1;
    hist(nonvisPars{iExp}(:), x);
    xlim(x([1 end]))
    xlabel('Nonvisual influence')
    ylabel('#Neurons*stimuli')
    title(sprintf('Dataset %d: nonvisual gain (each stimulus)', iExp))
end
% (3) Baseline
for iExp = 1:length(results)
    figure
    x = -4:.1:4;
    hist(basePars{iExp}(:), x);
    xlim(x([1 end]))
    xlabel('Baseline influence')
    ylabel('#Neurons*stimuli')
    title(sprintf('Dataset %d: Baseline gain (each stimulus)', iExp))
end
% (4) Interaction
for iExp = 1:length(results)
    figure
    x = -.4:.05:.4;
    hist(interactPars{iExp}(:), x);
    xlim(x([1 end]))
    xlabel('Interaction term')
    ylabel('#Neurons')
    title(sprintf('Dataset %d: interaction', iExp))
end

% now pool across datasets
bestModel = cat(1, bestModel{:});
isGad = cat(1, isGad{:});
stimPars = cat(1, stimPars{:});
nonvisPars = cat(1, nonvisPars{:});
basePars = cat(1, basePars{:});
interactPars = cat(1, interactPars{:});
nonvisual = cat(1, nonvisual{:});
baseline = cat(1, baseline{:});
gadGroups = [-1 0 1];

% relationship between parameters and variables
x = -1:.1:30;
figure
subplot(2,1,1)
plot(baseline(:), basePars(:), 'k.')
xlim(x([1 end]))
ylabel('Baseline gain')
title('Baseline vs gain (each neuron and stimulus)')
subplot(2,1,2)
hist(baseline(:), x)
xlim(x([1 end]))
xlabel('Baseline (\DeltaF/F)')
ylabel('# Neurons*stimuli')

x = 0:.1:17;
figure
subplot(2,1,1)
plot(nonvisual(:), nonvisPars(:), 'k.')
xlim(x([1 end]))
ylabel('Nonvisual gain')
title('Nonvisual vs gain (each neuron and stimulus)')
subplot(2,1,2)
hist(nonvisual(:), x)
xlim(x([1 end]))
xlabel('Nonvisual signal')
ylabel('# Neurons*stimuli')

x = -1:.5:60;
figure
subplot(2,1,1)
plot(mean(nonvisual.*baseline,2), interactPars, 'k.')
xlim(x([1 end]))
ylabel('Interaction term')
title('Nonvisual*baseline vs gain (per neuron)')
subplot(2,1,2)
hist(mean(nonvisual.*baseline,2), x)
xlim(x([1 end]))
xlabel('Nonvisual signal * baseline')
ylabel('# Neurons')

% relationship between parameters
figure
plot(nonvisPars(:), basePars(:), 'k.')
xlabel('Nonvisual gain')
ylabel('Baseline gain')
xlim([-5 5])
ylim([-10 10])

figure
plot(stimPars(:), nonvisPars(:), 'k.')
xlabel('Stimulus influence (intercepts)')
ylabel('Nonvisual gain')
xlim([-10 35])
ylim([-5 5])

figure
plot(stimPars(:), basePars(:), 'k.')
xlabel('Stimulus influence (intercepts)')
ylabel('Baseline gain')
xlim([-10 35])
ylim([-10 10])

% relationship between variables
figure
plot(nonvisual(:), baseline(:), 'k.')
xlabel('Nonvisual signal')
ylabel('Baseline')
xlim([-5 5])
ylim([0 30])

% parameters for each model
for iMod = 1%1:length(mods)
    ind = true(size(bestModel)); %bestModel == iMod;
    gad = isGad(ind);
    
    stim = stimPars(ind,:);
    if ~all(isnan(stim(:)))
        figure
        x = -2:.1:5;
%         x = -10:30;
        n = zeros(length(x), 3);
        for g = 1:3
            n(:,g) = hist(reshape(stim(gad==gadGroups(g),:),[],1), x);
        end
        bar(x, n, 'stacked')
        legend('exc','NA','inh')
        xlim(x([1 end]))
        xlabel('Intercepts (stimulus influence)')
        ylabel('# Neurons*stimuli')
        title(sprintf('%s: intercepts (each stimulus)', mods{iMod}))
    end
    nv = nonvisPars(ind,:);
    if ~all(isnan(nv(:)))
        figure
        x = -.8:.05:.8;
        n = zeros(length(x), 3);
        for g = 1:3
            n(:,g) = hist(reshape(nv(gad==gadGroups(g),:),[],1), x);
        end
        bar(x, n, 'stacked')
        legend('exc','NA','inh')
        xlim(x([1 end]))
        xlabel('Nonvisual influence')
        ylabel('# Neurons*stimuli')
        title(sprintf('%s: nonvisual gain (each stimulus)', mods{iMod}))
        figure
        n = zeros(length(x), 3);
        for g = 1:3
            n(:,g) = hist(mean(nv(gad==gadGroups(g),:),2), x);
        end
        bar(x, n, 'stacked')
        xlim(x([1 end]))
        xlabel('Nonvisual influence')
        ylabel('# Neurons')
        title(sprintf('%s: nonvisual gain (stimulus average)', mods{iMod}))
    end
    bl = basePars(ind,:);
    if ~all(isnan(bl(:)))
        figure
        x = -3:.1:3;
        n = zeros(length(x), 3);
        for g = 1:3
            n(:,g) = hist(reshape(bl(gad==gadGroups(g),:),[],1), x);
        end
        bar(x, n, 'stacked')
        legend('exc','NA','inh')
        xlim(x([1 end]))
        xlabel('Baseline influence')
        ylabel('# Neurons*stimuli')
        title(sprintf('%s: baseline gain (each stimulus)', mods{iMod}))
        figure
        n = zeros(length(x), 3);
        for g = 1:3
            n(:,g) = hist(mean(bl(gad==gadGroups(g),:),2), x);
        end
        bar(x, n, 'stacked')
        legend('exc','NA','inh')
        xlim(x([1 end]))
        xlabel('Baseline influence')
        ylabel('# Neurons')
        title(sprintf('%s: baseline gain (stimulus average)', mods{iMod}))
    end
    ia = interactPars(ind);
    if ~all(isnan(ia))
        figure
        x = -.4:.05:.4;
        n = zeros(length(x), size(ia,2), 3);
        for g = 1:3
            n(:,:,g) = hist(ia(gad==gadGroups(g),:), x);
        end
        bar(x, squeeze(sum(n,2)), 'stacked')
        legend('exc','NA','inh')
        xlim(x([1 end]))
        xlabel('Interaction term')
        ylabel('# Neurons')
        title(sprintf('%s: interaction nonvisual and baseline', mods{iMod}))
    end
end
