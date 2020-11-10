%% Load database
% db_driftingGratings
db_boutons_driftingGratings_blanks

%% Folders
% loading data
folderROIData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\InfoStructs';
% saving data
% folderResults = 'C:\RESULTS\nonVisualEffects\modelGratingResp\';
folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\boutons\nonVisualEffects\';

%% Parameters
% nonVisualSignal = 'running'; % 'running' or 'pupil'
doShuffle = 1; % test whether calcium trace is better fit using a kernel
               % decribing stimulus response or just a baseline
doSave = 1; % if 1, save model results and plots (if plotting)
doPlot = 1; % if 1, plot results for each cell (kernel fit and relationship
            % between nonvisual signal, stimulus response and baseline
doPause = 0; % if 1, pause after each cell

% determine which model to use
model = 'iterative'; % 'SVD' or 'iterative'
kernelLength = 15; %14; % length of kernel in sec (only for 'iterative' model)
basisLength = 0; %25; % length of cosine basis functions to model baseline in 
                  % sec (only for 'iterative' model)

% response times to consider before and after stimulus (in s)
timeBeforeOnset = 2;
timeAfterOffset = 3;

% model evaluation
minR2 = 0.3;
minEV = 0;

% find neurons that are well correlated with non-visual data
% minRho = 0.2; % use positive and negative correlations
% maxP = .05;

% fitting tuning curves
% minRSquared = .2;

results = db;

%% Loop across datasets to fit models to single cells
% adapted to incorporate new subset of cells (isDuplicate & isSwitchOn)
% data = load(fullfile(folderResults, 'kernelFit', 'results.mat'));
% results = data.results;
%
for k=1:length(db)
    fprintf('Dataset %d: %s %s exp.: %d\n', k, db(k).subject, ...
        db(k).date, db(k).expGratings);
    folder = [folderROIData filesep db(k).subject filesep ...
        db(k).date filesep num2str(db(k).expGratings)];
    fileStart = [db(k).date '_' num2str(db(k).expGratings) '_' ...
        db(k).subject];
    file = [fileStart '_2P_plane%03d_ROI.mat'];
    for iPlane = 1:length(db(k).planes)
        % load meta
        data=load(fullfile(folder, sprintf(file,db(k).planes(iPlane))));
        meta = data.meta;

        [stimTimes, stimSeq, stimMatrix, frameTimes, samplingRate] = ...
            ssLocal.getStimulusResponseInfo(meta);
        stimDuration = mean(stimTimes.offset - stimTimes.onset);
        [~,blankStim] = gratings.getOrientations(stimSeq);
        
        % only consider ROIs that are unique
        if isfield(meta.ROI, 'isDuplicate')
            cellIDs = find(meta.ROI.isDuplicate == 0 & meta.ROI.isSwitchOn == 0 & ...
                ~all(isnan(meta.F_final),1)');
        else
            cellIDs = find(~all(isnan(meta.F_final),1)');
        end
        calciumTraces = meta.F_final(:,cellIDs);
%         
%         old = results(k).plane(iPlane).cellIDs;
%         new = find(meta.ROI.isDuplicate == 0 & meta.ROI.isSwitchOn == 0);
%         toAdd = new(~ismember(new, old));
%         toRemove = find(~ismember(old, new));
%         cellIDs = new;
%         calciumTraces = meta.F_final(:,toAdd);
%         cellIDs = find(meta.ROI.isDuplicate == 0);

        offsets = round([timeBeforeOnset, timeAfterOffset] .* samplingRate);
        respBothEffects = ssLocal.getTracesPerStimulus(calciumTraces, ...
            stimMatrix, offsets);
%         % NOTE: NaN values occur in responses if offsets for first or last
%         % stimuli are outside of the recording time -> set them to zero
%         respBothEffects(isnan(respBothEffects)) = 0;
        responseTime = ((0:size(respBothEffects,4)-1)-offsets(1)) ./ samplingRate;
        % respBothEffects: [neurons x stimuli x trials x time]
        
        kernelSamples = round(kernelLength * samplingRate);
        kernelTime = (1:kernelSamples) ./ samplingRate;
        basisSamples = round(basisLength * samplingRate);
        stimSamples = round(stimDuration * samplingRate);
        
        results(k).plane(iPlane).cellIDs = cellIDs;
        results(k).plane(iPlane).stimDuration = stimDuration;
        results(k).plane(iPlane).kernelTime = kernelTime;
%         
%         results(k).plane(iPlane).responses = cat(1, ...
%             results(k).plane(iPlane).responses, ...
%             permute(respBothEffects, [1 4 3 2]));
        results(k).plane(iPlane).responses = permute(respBothEffects, [1 4 3 2]);
%
        results(k).plane(iPlane).responseTime = responseTime;
%         results(k).plane(iPlane).isGad = meta.ROI.isGad(cellIDs);
        
        if doSave==1 && doPlot==1
            fPlots = fullfile(folderResults, 'plots_kernelFit', ...
                [db(k).subject '_' db(k).date '_' num2str(db(k).expGratings)], ...
                num2str(db(k).planes(iPlane)));
            if ~exist(fPlots, 'dir')
                mkdir(fPlots);
            end
        end
        
        display(['Plane ' num2str(db(k).planes(iPlane))])
        
        for iCell = 1:size(calciumTraces, 2)
            nChars = fprintf('  Cell %d of %d', iCell, size(calciumTraces, 2));
            
            cellResponse = permute(respBothEffects(iCell,:,:,:), [4 3 2 1]);
            % [time x trials x stimuli]
            switch model
                case 'SVD'
                    fitResults = models.fitKernelWithSVD(cellResponse, ...
                        offsets(1));
                    % calculate adjusted R^2 for kernel fit
                    adjR2 = models.getAdjR2(cellResponse, fitResults);
                    % calculate prediction for each trial
                    predKernel = bsxfun(@plus, bsxfun(@times, fitResults.kernel, ...
                        permute(fitResults.alphaEachTrial, [3 1 2])), ...
                        permute(fitResults.baselineEachTrial, [3 1 2]));
                    baseL = [];
                case 'iterative'
                    [fitResults, r2] = ...
                        models.fitKernelIteratively(calciumTraces(:,iCell), ...
                        stimMatrix, blankStim, kernelSamples, basisSamples, ...
                        offsets(1), doPlot, doShuffle, stimSamples);
                    if doPlot == 1
                        if doSave == 1
                            fig = gcf;
                            fig.PaperPositionMode = 'auto';
                            print(fullfile(fPlots, sprintf('kernelModel%03d.jpg', iCell)), ...
                                '-djpeg','-r0')
                        end
                        if doPause == 1
                            pause
                        end
                        close gcf
                    end
                    if isempty(fitResults.kernel)
                        predKernel = [];
                    else
                        predKernel = ssLocal.getTracesPerStimulus( ...
                            fitResults.prediction(:), stimMatrix, offsets);
                        predKernel = permute(predKernel, [4 3 2 1]);
                    end
                    if isempty(fitResults.baseline) || basisLength == 0
                        baseL = [];
                    else
                        baseL = ssLocal.getTracesPerStimulus( ...
                            fitResults.baseline(:), stimMatrix, offsets);
                        baseL = permute(baseL, [4 3 2 1]);
                    end
                    adjR2.kernel = r2;
            end
            
            
            % store all results
%             
%             results(k).plane(iPlane).kernelFit(iCell+length(old)) = fitResults;
%             results(k).plane(iPlane).adjR2(iCell+length(old)) = adjR2;
            results(k).plane(iPlane).kernelFit(iCell) = fitResults;
            results(k).plane(iPlane).adjR2(iCell) = adjR2;
%
            
            if doPlot == 1
                % check validity of kernel approach
                nonVis.plotRawModelledResponses(cellResponse, predKernel, baseL, ...
                    responseTime, stimDuration, sprintf('Adjusted R2: %.2f',adjR2.kernel));
                
                if doSave==1
                    fig = gcf;
                    fig.PaperPositionMode = 'auto';
                    print(fullfile(fPlots, sprintf('kernelFit%03d.jpg', iCell)), ...
                        '-djpeg','-r0')
                end
                close all
            end
            fprintf(repmat('\b', 1, nChars));
        end
        
%         
%         results(k).plane(iPlane).kernelFit(toRemove) = [];
%         results(k).plane(iPlane).adjR2(toRemove) = [];
%         results(k).plane(iPlane).responses(toRemove,:,:,:) = [];
%         cellIDs = old;
%         cellIDs(toRemove) = [];
%         [~,order] = sort([cellIDs; toAdd], 'ascend');
%         results(k).plane(iPlane).kernelFit = results(k).plane(iPlane).kernelFit(order);
%         results(k).plane(iPlane).adjR2 = results(k).plane(iPlane).adjR2(order);
%         results(k).plane(iPlane).responses = results(k).plane(iPlane).responses(order,:,:,:);
%         
        
        if doSave == 1
            if ~exist(fullfile(folderResults, 'kernelFit'), 'dir')
                mkdir(fullfile(folderResults, 'kernelFit'));
            end
            save(fullfile(folderResults, 'kernelFit', 'results.mat'), 'results');
        end
    end
end

%% Add cell type information to results structure
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
        results(k).plane(iPlane).isGad = meta.ROI.isGad( ...
            results(k).plane(iPlane).cellIDs);
    end
end