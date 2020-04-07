% label = 'neurons';
label = 'boutons';

%% Load database
if strcmp(label, 'neurons')
    db_sparseNoise
else
    db_boutons_sparseNoise
end

%% Folders
% file location depending on PCs
folderPC = 'C:\STORAGE\OneDrive - University College London';
% folderPC = 'C:\Users\skgtchr\OneDrive - University College London'; % ZEN notebook
% loading data
folderROIData = fullfile(folderPC, 'Lab\DATA\InfoStructs');
% saving data
if strcmp(label, 'neurons')
    folderResults = fullfile(folderPC, 'Lab\RESULTS\receptiveFields\SC neurons');
else
    folderResults = fullfile(folderPC, 'Lab\RESULTS\receptiveFields\boutons');
end

%% Parameters
% for correcting baseline drifts of calcium traces at start of experiments
driftWin = 20; % in s, window to test whether baseline is higher than normal
driftThresh = 1.5; % in std, threshold for drift
correctWin = 150; % in s, window to fit exponential

% for receptive field estimates
% used for fitting 2 RFs (ON and OFF simultaneously), and fitting running
% kernels and RFs simultaneously
lambdasStim = logspace(-4, 1, 6);
lambdasRun = logspace(0, 6, 7);
% RFlimits = [0.1 0.6];
RFlimits = [0.2 0.4];
crossFolds = 10;
rfTypes = {'ON','OFF'};

% parameters for running speed as predictor
runKrnlLimits = [-5 5];

% colormap
grad = linspace(0,1,100)';
reds = [ones(100,1),grad,grad];
blues = [grad,grad,ones(100,1)];
cm = [blues; flip(reds(1:end-1,:),1)];

% for evaluation of receptive fields (significance/goodness)
minExplainedVarianceStim = 0.01;
minPVal = 0.05;
maxLambda = 1;

%% Fit RFs and get cross-validated explained variance

RFs = db;
RF_predictions = db;
% gDev = gpuDevice;
for k = 1:length(db)
%     reset(gDev);
    fprintf('Dataset %d: %s %s exp.: %d\n', k, db(k).subject, ...
        db(k).date, db(k).exp);
    folder = fullfile(folderROIData, db(k).subject, ...
        db(k).date, num2str(db(k).exp));
    file = [sprintf('%s_%d_%s',db(k).date, db(k).exp, ...
        db(k).subject) '_2P_plane%03d_ROI.mat'];
    
    traces = [];
    numNeurons = NaN(1, length(db(k).planes));
    for iPlane = 1:length(db(k).planes)
        % load meta
        data=load(fullfile(folder, sprintf(file,db(k).planes(iPlane))));
        meta = data.meta;
        meta.folderTL = strrep(meta.folderTL, 'ioo.ucl.ac.uk', ...
            'cortexlab.net');
        
        % get stimulus information and running data
        if iPlane == 1
            stimTimes = ppbox.getStimTimes(meta);
            [stimFrames, stimPosition] = whiteNoise.getStimulusFrames(meta);
            stimFrameDur = mean(stimTimes.offset - stimTimes.onset) / size(stimFrames,3);
            stimFrameTimes = ((1:size(stimFrames,3))-1) .* stimFrameDur;
            allStimFrameTimes = reshape((stimTimes.onset + stimFrameTimes)', [], 1);
            RFtimesInFrames = floor(RFlimits(1) / stimFrameDur) : ...
                ceil(RFlimits(2) / stimFrameDur);
            
            % get trace time
            traceTimes = ppbox.getFrameTimes(meta);
            traceBin = median(diff(traceTimes));
            
            % get running speed
            ballData = nonVis.getRunningSpeed(meta);
            if isempty(ballData)
                disp('No ball data!')
                return
            end
            runSpeed = ballData.total / median(diff(ballData.t)) / 53;
            runTime = ballData.t;
        end
        
        %get calcium traces
        tr = meta.F_final;
        % set data of duplicate neurons or unhealthy neurons to NaN
        if isfield(meta.ROI, 'isDuplicate')
            tr(:, meta.ROI.isDuplicate==1) = NaN;
        end
        if isfield(meta.ROI, 'isSwitchOn')
            tr(:, meta.ROI.isSwitchOn==1) = NaN;
        end
        if iPlane == 1
            traces = [traces, tr];
            numNeurons(1) = size(tr,2);
        else
            tr_int = NaN(length(traceTimes), size(tr,2));
            t = ppbox.getFrameTimes(meta);
            for n = 1:size(tr,2)
                if all(isnan(tr(:,n)))
                    continue
                end
                nanInd1 = isnan(tr(:,n));
                tr_int(:,n) = interp1(t(~nanInd1), tr(~nanInd1,n), traceTimes, 'pchip');
                nanInd2 = hist(t(nanInd1), traceTimes) > 0;
                tr_int(nanInd2,n) = NaN;
            end
            traces = [traces, tr_int];
            numNeurons(iPlane) = size(tr,2);
        end
    end
    % remove strong baseline decay at start of experiment in cells that
    % show it
    indCells = find(nanmean(traces(1:round(driftWin/traceBin),:),1) > ...
        nanmean(traces,1) + driftThresh.*nanstd(traces,0,1));
    ind = round(correctWin / traceBin);
    for iCell = 1:length(indCells)
        y = traces(:,indCells(iCell));
        y = fillmissing(y, 'linear');
        f = fit((1:length(y))', y, ...
            @(a,b,c,d,e,x) a+b.*exp(-x./c)+d.*exp(-x./e), ...
            'Lower', [0 0 0 0 0], ...
            'Upper', [max(y) max(y) 500 max(y) 500], ...
            'StartPoint', [min(y) mean(y) 50 mean(y) 5]);
        traces(:,indCells(iCell)) = y - f(1:size(traces,1)) + f.a;
    end
    
    % only show left side of stimulus
    stimFrames = stimFrames(:,1:min(34, size(stimFrames,2)),:);
    
    [rFields, runKernels, runWin, ev, ev_run, ev_stim, ...
        preds, preds_run, time] = ...
        whiteNoise.getReceptiveField_OnOff_smooth_running( ...
        traces, traceTimes, stimFrames, stimFrameTimes, stimTimes, ...
        RFtimesInFrames, runSpeed, runTime, runKrnlLimits, ...
        {lambdasStim, lambdasRun}, crossFolds, false);
    RFs(k).RFTimes = RFtimesInFrames * stimFrameDur;
    RFs(k).stimPosition = stimPosition;
    RFs(k).runWindow = runWin;
    RF_predictions(k).lambdasStim = lambdasStim;
    RF_predictions(k).lambdasRun = lambdasRun;
    RF_predictions(k).time = time;
    
    n = 0;
    for iPlane = 1:length(db(k).planes)
        indNeurons = (1:numNeurons(iPlane))+n;
        
        v = squeeze(mean(ev(indNeurons,:,:,:),2)); % [neuron x lamStim x lamRun], average across cross-folds
        [maxEV,maxStimLam] = max(v,[],2);
        maxEV = squeeze(maxEV); % [neuron x lamRun];
        maxStimLam = squeeze(maxStimLam); % [neuron x lamRun];
        [maxEV, maxRunLam] = max(maxEV, [], 2); % [neuron x 1]
        indLam = sub2ind(size(maxStimLam), (1:size(maxStimLam,1))', maxRunLam);
        maxStimLam = maxStimLam(indLam); % [neuron x 1]
        
        vRun = squeeze(mean(ev_run(indNeurons,:,:,:),2)); % [neuron x lamStim x lamRun], average across cross-folds
        vStim = squeeze(mean(ev_stim(indNeurons,:,:,:),2)); % [neuron x lamStim x lamRun]
        inds = sub2ind(size(vRun), (1:size(vRun,1))', maxStimLam, maxRunLam);
        maxEVRun = vRun(inds); % [neuron x 1]
        maxEVStim = vStim(inds); % [neuron x 1]
        
        RF_predictions(k).plane(iPlane).explainedVariances = ev(indNeurons,:,:,:);
        RF_predictions(k).plane(iPlane).predictions = preds(:,indNeurons);
        RF_predictions(k).plane(iPlane).explainedVariances_runOnly = ev_run(indNeurons,:,:,:);
        RF_predictions(k).plane(iPlane).predictions_runOnly = preds_run(:,indNeurons);
        RF_predictions(k).plane(iPlane).explainedVariances_stimOnly = ev_stim(indNeurons,:,:,:);
        
        RFs(k).plane(iPlane).receptiveFields = rFields(:,:,:,:,indNeurons);
        RFs(k).plane(iPlane).runningKernels = runKernels(:,indNeurons);
        RFs(k).plane(iPlane).explainedVariances = maxEV;
        RFs(k).plane(iPlane).lambdasStim = lambdasStim(maxStimLam);
        RFs(k).plane(iPlane).lambdasRun = lambdasRun(maxRunLam);
        RFs(k).plane(iPlane).explainedVariances_runOnly = maxEVRun;
        RFs(k).plane(iPlane).explainedVariances_stimOnly = maxEVStim;
        
        n = n + numNeurons(iPlane);
    end

    save(fullfile(folderResults, 'receptiveFields.mat'), 'RFs')
    save(fullfile(folderResults, 'receptiveFields_predictions.mat'), 'RF_predictions')
    fprintf('\n')
end

%% Test significance of receptive field
data = load(fullfile(folderResults, 'receptiveFields.mat'));
RFs = data.RFs;

gDev = gpuDevice;
reset(gDev);

for k = 1:length(RFs)
    fprintf('Dataset %d: %s %s exp.: %d\n', k, RFs(k).subject, ...
        RFs(k).date, RFs(k).exp);
    folder = fullfile(folderROIData, RFs(k).subject, ...
        RFs(k).date, num2str(RFs(k).exp));
    file = [sprintf('%s_%d_%s',RFs(k).date, RFs(k).exp, ...
        RFs(k).subject) '_2P_plane%03d_ROI.mat'];
    
    fprintf('  Plane (of %d):\n', length(RFs(k).planes))
    for iPlane = 1:length(RFs(k).planes)
        % load meta
        data=load(fullfile(folder, sprintf(file,RFs(k).planes(iPlane))));
        meta = data.meta;
        meta.folderTL = strrep(meta.folderTL, 'ioo.ucl.ac.uk', ...
            'cortexlab.net');
        
        % get stimulus information
        if iPlane == 1
            stimTimes = ppbox.getStimTimes(meta);
            [stimFrames, stimPosition] = whiteNoise.getStimulusFrames(meta);
            stimFrameDur = mean(stimTimes.offset - stimTimes.onset) / size(stimFrames,3);
            stimFrameTimes = ((1:size(stimFrames,3))-1) .* stimFrameDur;
            allStimFrameTimes = reshape((stimTimes.onset + stimFrameTimes)', [], 1);
            RFtimesInFrames = floor(RFlimits(1) / stimFrameDur) : ...
                ceil(RFlimits(2) / stimFrameDur);
        
            % only show left side of stimulus
            stimFrames = stimFrames(:,1:min(34, size(stimFrames,2)),:);
            
            % get trace time
            traceTimes = ppbox.getFrameTimes(meta);
            traceBin = median(diff(traceTimes));
            
            % get running speed
            ballData = nonVis.getRunningSpeed(meta);
            if isempty(ballData)
                disp('No ball data!')
                return
            end
            runSpeed = ballData.total / median(diff(ballData.t)) / 53;
            runTime = ballData.t;
            fprintf('    ')
        end
        
        fprintf(' %d', iPlane)
        
        %get calcium traces
        tr = meta.F_final;
        % set data of duplicate neurons or unhealthy neurons to NaN
        if isfield(meta.ROI, 'isDuplicate')
            tr(:, meta.ROI.isDuplicate==1) = NaN;
        end
        if isfield(meta.ROI, 'isSwitchOn')
            tr(:, meta.ROI.isSwitchOn==1) = NaN;
        end
        if all(isnan(tr(:)))
            RFs(k).plane(iPlane).pVal_RFonly = NaN(size(tr,2),1);
            continue
        end
        if iPlane == 1
            traces = tr;
        else
            tr_int = NaN(length(traceTimes), size(tr,2));
            t = ppbox.getFrameTimes(meta);
            for n = 1:size(tr,2)
                if all(isnan(tr(:,n)))
                    continue
                end
                nanInd1 = isnan(tr(:,n));
                tr_int(:,n) = interp1(t(~nanInd1), tr(~nanInd1,n), traceTimes, 'pchip');
                nanInd2 = hist(t(nanInd1), traceTimes) > 0;
                tr_int(nanInd2,n) = NaN;
            end
            traces = tr_int;
        end
        
        % remove strong baseline decay at start of experiment in cells that
        % show it
        indCells = find(nanmean(traces(1:round(driftWin/traceBin),:),1) > ...
            nanmean(traces,1) + driftThresh.*nanstd(traces,0,1));
        ind = round(correctWin / traceBin);
        for iCell = 1:length(indCells)
            y = traces(:,indCells(iCell));
            y = fillmissing(y, 'linear');
            f = fit((1:length(y))', y, ...
                @(a,b,c,d,e,x) a+b.*exp(-x./c)+d.*exp(-x./e), ...
                'Lower', [0 0 0 0 0], ...
                'Upper', [max(y) max(y) 500 max(y) 500], ...
                'StartPoint', [min(y) mean(y) 50 mean(y) 5]);
            traces(:,indCells(iCell)) = y - f(1:size(traces,1)) + f.a;
        end
        
        [ev, ev_shift] = ...
            whiteNoise.receptiveFieldShiftTest_OnOff( ...
            traces, traceTimes, stimFrames, stimFrameTimes, stimTimes, ...
            RFtimesInFrames, runSpeed, runTime, ...
            RFs(k).plane(iPlane).runningKernels, RFs(k).runWindow, ...
            RFs(k).plane(iPlane).receptiveFields, RFs(k).plane(iPlane).lambdasStim, 500);
        pvals = sum(ev_shift > ev, 2) ./ size(ev_shift,2);
        pvals(isnan(RFs(k).plane(iPlane).explainedVariances)) = NaN;
        RFs(k).plane(iPlane).pVal_RFonly = pvals;
    end

    save(fullfile(folderResults, 'receptiveFields_wSignificance.mat'), 'RFs')
    fprintf('\n')
end

%% Plot RF of each cell
data = load(fullfile(folderResults, 'receptiveFields.mat'));
RFs = data.RFs;
data = load(fullfile(folderResults, 'receptiveFields_predictions.mat'));
RF_preds = data.RF_predictions;

colors = lines(2);

for k = 1:length(RFs)
    fprintf('Dataset %d: %s %s exp.: %d\n', k, RFs(k).subject, ...
        RFs(k).date, RFs(k).exp);
    folder = fullfile(folderROIData, RFs(k).subject, ...
        RFs(k).date, num2str(RFs(k).exp));
    file = [sprintf('%s_%d_%s',RFs(k).date, RFs(k).exp, ...
        RFs(k).subject) '_2P_plane%03d_ROI.mat'];
    for iPlane = 1:length(RFs(k).plane)
        % load meta
        data=load(fullfile(folder, sprintf(file,RFs(k).planes(iPlane))));
        meta = data.meta;
        meta.folderTL = strrep(meta.folderTL, 'ioo.ucl.ac.uk', ...
            'cortexlab.net');
        traces = meta.F_final;
        traceTimes = ppbox.getFrameTimes(meta);
        traceBin = median(diff(traceTimes));
        traces = (traces - nanmean(traces)) ./ nanstd(traces);
        % remove strong baseline decay at start of experiment in cells that
        % show it
        indCells = find(nanmean(traces(1:round(driftWin/traceBin),:),1) > ...
            nanmean(traces,1) + driftThresh.*nanstd(traces,0,1));
        for iCell = 1:length(indCells)
            y = traces(:,indCells(iCell));
            y = fillmissing(y, 'linear');
            f = fit((1:length(y))', y, ...
                @(a,b,c,d,e,x) a+b.*exp(-x./c)+d.*exp(-x./e), ...
                'Lower', [0 0 0 0 0], ...
                'Upper', [max(y) max(y) 500 max(y) 500], ...
                'StartPoint', [min(y) mean(y) 50 mean(y) 5]);
            traces(:,indCells(iCell)) = y - f(1:size(traces,1)) + f.a;
        end
        
        rFields = RFs(k).plane(iPlane).receptiveFields;
        runKernels = RFs(k).plane(iPlane).runningKernels;
        ev = RFs(k).plane(iPlane).explainedVariances;
        evRun = RFs(k).plane(iPlane).explainedVariances_runOnly;
        evStim = RFs(k).plane(iPlane).explainedVariances_stimOnly;
        pValues = RFs(k).plane(iPlane).pVal_RFonly;
        runWin = RFs(k).runWindow;
        time = RF_preds(k).time;
        preds = RF_preds(k).plane(iPlane).predictions;
        predsRun = RF_preds(k).plane(iPlane).predictions_runOnly;
        
        fPlots = fullfile(folderResults, 'plots_kernels', ...
            sprintf('%s_%s', RFs(k).subject, RFs(k).date), ...
            sprintf('Plane%02d',RFs(k).planes(iPlane)));
        if ~isfolder(fPlots)
            mkdir(fPlots)
        end
        for iCell = 1:length(RFs(k).plane(iPlane).explainedVariances)
            if isnan(RFs(k).plane(iPlane).explainedVariances(iCell))
                continue
            end
            figure('Position',[1921 1 1920 1123])
            cols = size(rFields,2);
            subplot(3,5,10+(1:5))
            plot(traceTimes, traces(:,iCell), 'k')
            hold on
            rlims = [min(reshape(runKernels(:,iCell),[],1)), ...
                max(reshape(runKernels(:,iCell),[],1))];
            rf = rFields(:,:,:,:,iCell);
            clim = max(abs(rf(:)));
            for m = 1:2
                r = [rf(:,:,:,m), zeros(size(rf,1),1,size(rf,3))];
                r = flip(r,3);
                r = reshape(r,size(r,1),[]);
                r(:,end) = [];
                subplot(3,5,(m-1)*5+(1:4))
                imagesc(r,[-clim clim])
                colorbar
                axis image
                set(gca,'box','off','XTick',[],'YTick',[])
                if m == 1
                    title(sprintf('explained variance: %.2f%% total, %.2f%% stim only (p = %.4f)', ...
                        ev(iCell)*100, evStim(iCell)*100, pValues(iCell)))
                else
                    set(gca,'XTick',ceil(cols/2):(cols+1):(cols+1)*size(rFields,3), ...
                        'XTickLabel',sprintf('%.2f\n',-flip(RFs(k).RFTimes)))
                    xlabel('Time from response (in s)')
                end
                ylabel(rfTypes{m})
            end
            colormap(cm)
            
            subplot(3,5,5)
            plot(runWin, runKernels(:,iCell), 'Color', colors(2,:), 'LineWidth', 2)
            axis tight
            ylim(rlims)
            set(gca, 'box', 'off')
            title(sprintf('explained variance: %.2f%% running only', ...
                evRun(iCell)*100))
            xlabel('Time from response (in s)')
            
            subplot(3,5,10+(1:5))
            h = [0 0];
            h(1) = plot(time, preds(1:length(time), iCell), 'Color', colors(1,:), 'LineWidth', 1);
            h(2) = plot(time, predsRun(1:length(time), iCell), 'Color', colors(2,:), 'LineWidth', 1);
            axis tight
            set(gca,'box','off')
            legend(h, {'stim + run','run only'})
            xlabel('Time (in s)')
            ylabel('\DeltaF/F')
            annotation('textbox', [0 .95 1 .03], 'String', sprintf('Neuron %d', iCell), ...
                'FontSize', 14, 'FontWeight', 'bold', 'LineStyle', 'none', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(fullfile(fPlots, sprintf('Neuron%03d.jpg', iCell)), ...
                '-djpeg','-r0')
            close gcf
        end
    end
end

%% Plot distribution of ON-OFF-ratios
data = load(fullfile(folderResults, 'receptiveFields.mat'));
RFs = data.RFs;

datasets = [];
planes = [];
neurons = [];
evStim = [];
lambdasStim = [];
pValues = [];
OnOffValues = [];

for k = 1:length(RFs)
    for iPlane = 1:length(RFs(k).plane)
        n = length(RFs(k).plane(iPlane).explainedVariances);
        datasets = [datasets; ones(n,1).*k];
        planes = [planes; ones(n,1).*iPlane];
        neurons = [neurons; (1:n)'];
        evStim = [evStim; RFs(k).plane(iPlane).explainedVariances_stimOnly];
        lambdasStim = [lambdasStim; RFs(k).plane(iPlane).lambdasStim'];
        pValues = [pValues; RFs(k).plane(iPlane).pVal_RFonly];
        
        oov = NaN(n,2);
        for iCell = 1:n
            rf = RFs(k).plane(iPlane).receptiveFields(:,:,:,:,iCell);
            [~,t] = max(max(reshape(permute(abs(rf),[1 2 4 3]), [], ...
                size(rf,3)), [], 1));
            rf = squeeze(rf(:,:,t,:));
            [mx,row] = max(max(abs(rf),[],3),[],1);
            [~,col] = max(mx);
            row = row(col);
            rf = squeeze(rf(row,col,:));
            rf(2) = -rf(2);
            oov(iCell,:) = rf;
        end
        OnOffValues = [OnOffValues; oov];
    end
end

[~,type] = max(abs(OnOffValues),[],2);
ind = sub2ind(size(OnOffValues), (1:size(OnOffValues,1))', type);
signs = sign(OnOffValues(ind));
ratios = OnOffValues;
ratios(signs<0,:) = -ratios(signs<0,:);
ratios(ratios<0) = 0;
ratios = (ratios(:,1)-ratios(:,2))./sum(ratios,2);

% Plot histogram of On-Off-ratios
binSize = 0.05;
ind = pValues < minPVal & evStim > minExplainedVarianceStim & lambdasStim < maxLambda;
edges = -1:binSize:1;
bins = edges(1:end-1)+binSize/2;
figure
histogram(ratios(ind), edges, 'FaceColor', 'k')
set(gca, 'box', 'off')
xlim([-1 1])
xlabel('ON/OFF-ratio ((ON-OFF)/(ON+OFF))')
ylabel(sprintf('#%s',label))
title(sprintf('ON/OFF ratios (n=%d)', sum(ind)))

% Plot On peak versus Off peak values
ind = [pValues > -1, pValues < 0.05, pValues < 0.05 & evStim > 0.01, ...
    pValues < 0.05 & evStim > 0.015 & lambdasStim < 1];
groups = {'all', 'p<0.05', 'p<0.05 & ev(stim)>0.015', ...
    'p<0.05 & ev(stim)>0.015 & lambda(stim)<1'};
for g = 1:length(groups)
    figure
    scatter(OnOffValues(ind(:,g),1), OnOffValues(ind(:,g),2), 'filled', ...
        'MarkerFaceAlpha', 0.2)
    hold on
    plot([0 0],[-.05 .2],'k')
    plot([-.05 .2],[0 0],'k')
    axis([-0.05 0.2 -0.05 0.2])
    xlabel('ON peak value')
    ylabel('OFF peak value')
    title(sprintf('%s, %s (n=%d)', label, groups{g}, sum(ind)))
end