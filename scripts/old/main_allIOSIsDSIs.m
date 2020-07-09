%% Load database
db_driftingGratings

%% Folders
folderROIData = 'C:\Temp2p';
folderTuningModelRunning = ['C:\Temp2p\Results\nonVisualEffects\tuningModelledRunning_' stimType];
folderTuningModelPupil = ['C:\Temp2p\Results\nonVisualEffects\tuningModelledPupil_' stimType];

%% Loop across datasets
osis = [];
dsis = [];
for k=1:length(db)
    
    %% Load data
%     folder = [folderROIData filesep db(k).subject filesep ...
%         db(k).date filesep num2str(db(k).exp)];
%     fileStart = [db(k).date '_' num2str(db(k).exp) '_' ...
%         db(k).subject];
%     file = [fileStart '_2P_plane%03d_ROI.mat'];
%     clear info
%     neurons = cell(1, length(db(k).planes));
%     neuronIDs = [];
%     for iPlane=1:length(db(k).planes)
%         % load meta
%         data=load(fullfile(folder, sprintf(file,db(k).planes(iPlane))));
%         info(iPlane) = orderfields(data.meta);
%         ind = find(strcmp(data.meta.ROI.CellClasses, 's'));
%         neurons{iPlane} = ind;
%     end
    
    % load results, parameters, respondingNeurons
    load(fullfile(folderTuningModelRunning, [db(k).subject '_' ...
        db(k).date '_' num2str(db(k).exp)], 'results.mat'));
%     window = parameters.window;
%     percentile = parameters.percentile;
    
    %% Preprocess calcium traces
    % interpolate neural responses to get a single series of time points
%     [calciumTraces, calciumTime] = ssLocal.matchSampleTimesAcrossPlanes(info);
%     sampleInt = median(diff(calciumTime));
%     for iPlane = 1:length(db(k).planes)
%         calciumTraces{iPlane} = calciumTraces{iPlane}(:,neurons{iPlane});
%     end
%     calciumTraces = cat(2, calciumTraces{:});
%     % caclulate F_0
%     [~,F_0] = ssLocal.removeSlowDrift(calciumTraces, calciumTime, window, percentile);
%     calciumTraces = (calciumTraces - F_0) ./ F_0;
%     
%     %% Get stimulus information
%     stimTimes = ppbox.getStimTimes(info(1));
%     stimSequence = ppbox.getStimSequence(info(1));
%     stimMatrix = ppbox.buildStimMatrix(stimSequence, stimTimes, calciumTime);
%     
%     %% Get response traces for each trial
%     timeLimits = [2 3]; % in sec; (1) before stim onset, (2) after stim offset
%     frameLimits = round(timeLimits ./ sampleInt);
%     stimTraces = ssLocal.getTracesPerStimulus(calciumTraces, stimMatrix, ...
%         frameLimits); % [neurons x stimuli x repetitions x timePoints]
%     sampleTimes = ((0:size(stimTraces,4)-1)-frameLimits(1)) .* sampleInt;
%     trialDur = (size(stimTraces,4)-1)*sampleInt;
%     stimDur = stimTimes.offset(1) - stimTimes.onset(1);
%     
%     trialResponses = mean(stimTraces,4);
%     respondingNeurons = zeros(size(trialResponses,1),1);
%     for iCell = 1:size(trialResponses,1)
% %         p = anova1(squeeze(trialResponses(iCell,:,:))');
%         p = anova1(squeeze(trialResponses(iCell,:,:))',[],'off');
%         if p < 0.05
%             respondingNeurons(iCell) = 1;
%         end
% %         pause
% %         close all
%     end
%     respondingNeurons = respondingNeurons(results.corr.order);

    ind = respondingNeurons==1 & results.errors(:,1)>parameters.minRSquared;
    o = NaN(size(results.osi,1),1);
    o(ind) = results.osi(ind,1);
    d = NaN(size(results.dsi,1),1);
    d(ind) = results.dsi(ind,1);
    osis = [osis;o];
    dsis = [dsis;d];
    
%     save(fullfile(folderTuningModelRunning, [db(k).subject '_' ...
%         db(k).date '_' num2str(db(k).exp)], 'results.mat'), ...
%         'respondingNeurons', '-append');
end

adjOsis = osis;
% adjOsis(osis==-1) = 0;
adjDsis = dsis;
% adjDsis(dsis==-1) = 0;
figure
counts = hist(adjOsis,0:0.1:1);
bar(0:0.1:1,counts,'k')
xlim([-.1 1.1])
set(gca,'box','off')
xlabel('OSI')
ylabel('#Neuron')
title(sprintf('n = %d (of %d in total)',sum(~isnan(adjOsis)),length(adjOsis)))
figure
counts = hist(adjDsis,0:0.1:1);
bar(0:0.1:1,counts,'k')
xlim([-.1 1.1])
set(gca,'box','off')
xlabel('DSI')
ylabel('#Neuron')
title(sprintf('n = %d (of %d in total)',sum(~isnan(adjDsis)),length(adjDsis)))