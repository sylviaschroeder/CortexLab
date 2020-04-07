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
% lambda = [0.05 0.1 0.5 1 5 10];
% lambda = [0.5 1 5 10 50 100];
% lambda = [0.5 1 5 10 50 100 500 1000];
% lambda = [0.1 0.5 1 5 10 50 100];
% lambdas = logspace(-2, 8, 11);
% used for fitting 4 RF models separately, and fitting running kernel and
% RF separately:
lambdasStim = logspace(-4, 1, 6);
lambdasRun = logspace(-2, 3, 6);

% RFlimits = [0.1 0.6];
RFlimits = [0.2 0.4];
crossFolds = 10;
models = {'linear','absolute','white','black'};
rfTypes = {'ON','OFF'};

% parameters for running speed as predictor
runKrnlLimits = [-5 5];

% colormap
grad = linspace(0,1,100)';
reds = [ones(100,1),grad,grad];
blues = [grad,grad,ones(100,1)];
cm = [blues; flip(reds(1:end-1,:),1)];

%% Fit RFs

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
    for iPlane = 1:min(2,length(db(k).planes))
        % load meta
        data=load(fullfile(folder, sprintf(file,db(k).planes(iPlane))));
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
            stimTimeStep = median(diff(stimFrameTimes));
            RFtimesInFrames = floor(RFlimits(1) / stimTimeStep) : ...
                ceil(RFlimits(2) / stimTimeStep);
            
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
    
    [rFields, runKernels, runWin, ev, ev_run, ev_stim, ev_shift, preds, preds_run, time] = ...
        whiteNoise.getReceptiveField_smooth_running( ...
        traces, traceTimes, stimFrames, stimFrameTimes, stimTimes, ...
        RFtimesInFrames, runSpeed, runTime, runKrnlLimits, ...
        {lambdasStim, lambdasRun}, crossFolds, models, true, true);
    RFs(k).RFTimes = RFtimesInFrames * stimTimeStep;
    RFs(k).stimPosition = stimPosition;
    RFs(k).models = models;
    RFs(k).runWindow = runWin;
    RF_predictions(k).lambdasStim = lambdasStim;
    RF_predictions(k).lambdasRun = lambdasRun;

    n = 0;
    for iPlane = 1:length(db(k).planes)
        ind = (1:numNeurons(iPlane))+n;
        
        v = squeeze(mean(ev(ind,:,:,:),2)); % [neuron x lambda x model], average across cross-folds
        [maxEV,maxLam] = max(v,[],2);
        maxEV = squeeze(maxEV); % [neuron x model];
        maxLam = squeeze(maxLam); % [neuron x model];
        
        vRun = squeeze(mean(ev_run(ind,:,:),2)); % [neuron x lambda], average across cross-folds
        if size(runKernels,3) > 1
            maxEVRun = NaN(numNeurons(iPlane), length(models)); % [neuron x model];
            for iCell = 1:numNeurons(iPlane)
                for m = 1:length(models)
                    maxEVRun(iCell,m) = vRun(iCell, maxLam(iCell,m), m);
                end
            end
        else
            [maxEVRun, maxLamRun] = max(vRun, [], 2);
        end
        p = NaN(size(preds,1),length(ind),length(models));
        pr = NaN(size(preds,1),length(ind));
        for iCell = 1:length(ind)
            for m = 1:length(models)
                p(:,iCell,m) = preds(:,ind(iCell),maxLam(iCell,m));
            end
            pr(:,iCell) = preds_run(:,ind(iCell),maxLamRun(iCell));
        end
        
        vStim = squeeze(mean(ev_stim(ind,:,:,:),2)); % [neuron x lambda x model]
        maxEVStim = NaN(numNeurons(iPlane),length(models)); % [neuron x model]
        for iCell = 1:numNeurons(iPlane)
            for m = 1:length(models)
                maxEVStim(iCell,m) = vStim(iCell,maxLam(iCell,m),m);
            end
        end
        vShift = squeeze(mean(ev_shift(ind,:,:,:),2)); % [neuron x model x shift]
        ciEVStim = prctile(vShift,[2.5 97.5],3); % [neuron x model x [2.5% 97.5%]]
        
        RF_predictions(k).plane(iPlane).explainedVariances = ev(ind,:,:,:);
        RF_predictions(k).plane(iPlane).predictions = p;
        RF_predictions(k).plane(iPlane).explainedVariances_runOnly = ev_run(ind,:,:);
        RF_predictions(k).plane(iPlane).predictions_runOnly = pr;
        RF_predictions(k).plane(iPlane).explainedVariances_stimOnly = ev_stim(ind,:,:,:);
        RF_predictions(k).plane(iPlane).explainedVariances_stimOnly_shifted = ev_shift(ind,:,:,:);
        
        RFs(k).plane(iPlane).receptiveFields = rFields(:,:,:,ind,:);
        RFs(k).plane(iPlane).explainedVariances = maxEV;
        RFs(k).plane(iPlane).lambdas = lambdasStim(maxLam);
        RFs(k).plane(iPlane).runningKernels = runKernels(:,ind,:);
        RFs(k).plane(iPlane).explainedVariances_runOnly = maxEVRun;
        RFs(k).plane(iPlane).lambdas_runOnly = lambdasRun(maxLamRun);
        RFs(k).plane(iPlane).explainedVariances_stimOnly = maxEVStim;
        RFs(k).plane(iPlane).explainedVariances_stimOnly_confInt = ciEVStim;
        
        % plot RFs
        fPlots = fullfile(folderResults, 'plots_kernels', ...
            sprintf('%s_%s', db(k).subject, db(k).date), ...
            sprintf('Plane%02d',db(k).planes(iPlane)));
        if ~isfolder(fPlots)
            mkdir(fPlots)
        end
        tr = traces(:,ind);
        tr = (tr - nanmean(tr)) ./ nanstd(tr);
        for iCell = 1:length(ind)
            if all(isnan(traces(:,ind(iCell))))
                continue
            end
            figure('Position',[1921 1 1920 1123])
%             clim = max(abs(reshape(rFields(:,:,:,ind(iCell),:),[],1)));
            cols = size(rFields,2);
            subplot(length(models)+1,5,length(models)*5+(1:5))
            plot(traceTimes, tr(:,iCell), 'k')
            hold on
            rlims = [min(reshape(runKernels(:,ind(iCell),:),[],1)), ...
                max(reshape(runKernels(:,ind(iCell),:),[],1))];
            for m = 1:length(models)
                r = [rFields(:,:,:,ind(iCell),m),zeros(size(rFields,1),1,size(rFields,3))];
                r = flip(r,3);
                r = reshape(r,size(r,1),[]);
                r(:,end) = [];
                clim = max(abs(r(:)));
                subplot(length(models)+1,5,(m-1)*5+(1:4))
                imagesc(r,[-clim clim])
                colorbar
                axis image
                set(gca,'box','off','XTick',[],'YTick',[])
                if m == length(models)
                    set(gca,'XTick',ceil(cols/2):(cols+1):(cols+1)*size(rFields,3), ...
                        'XTickLabel',sprintf('%.2f\n',-flip(RFs(k).RFTimes)))
                    xlabel('Time from response (in s)')
                end
                ylabel(models{m})
                title(sprintf('explained variance: %.2f%% total, %.2f%% stim only [%.2f%% %.2f%%]', ...
                    maxEV(iCell,m)*100, maxEVStim(iCell,m)*100, ...
                    ciEVStim(iCell,m,1)*100, ciEVStim(iCell,m,2)*100))
                
                if size(runKernels,3) > 1
                    subplot(length(models)+1,5,5*m)
                    plot(runWin, runKernels(:,ind(iCell),m), 'r', 'LineWidth', 2)
                    axis tight
                    ylim(rlims)
                    set(gca, 'box', 'off')
                    if m == 1
                        title('Running kernels')
                    elseif m == length(models)
                        xlabel('Time from response (in s)')
                    end
                end
            end
            colormap(cm)
            if size(runKernels,3) == 1
                subplot(length(models)+1,5,5)
                plot(runWin, runKernels(:,ind(iCell)), 'r', 'LineWidth', 2)
                axis tight
                ylim(rlims)
                set(gca, 'box', 'off')
                title(sprintf('explained variance: %.2f%% running only', ...
                    maxEVRun(iCell)*100))
                xlabel('Time from response (in s)')
            end
            [~,bestM] = max(maxEV(iCell,:),[],2);
            subplot(length(models)+1,5,length(models)*5+(1:5))
            h = [0 0];
            h(1) = plot(time, preds(:, ind(iCell), maxLam(iCell,bestM), ...
                bestM), 'LineWidth', 2);
            if size(runKernels,3) > 1
                h(2) = plot(time, preds_run(:, ind(iCell), maxLam(iCell,bestM), bestM), ...
                    'r', 'LineWidth', 2);
            else
                h(2) = plot(time, preds_run(:, ind(iCell), maxLamRun(iCell)), ...
                    'r', 'LineWidth', 2);
            end
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
        n = n + numNeurons(iPlane);
    end
    
    save(fullfile(folderResults, 'receptiveFields.mat'), 'RFs')
    save(fullfile(folderResults, 'receptiveFields_predictions.mat'), 'RF_predictions')
    fprintf('\n')
end

%% Check distribution of best lambda values
data = load(fullfile(folderResults, 'receptiveFields1-4.mat'));
RFs = data.RFs;
lambdasStim = [];
lambdasRun = [];
for k = 1:length(RFs)
    for iPlane = 1:length(RFs(k).plane)
        valid = ~all(isnan(RFs(k).plane(iPlane).explainedVariances),2);
        lambdasStim = [lambdasStim; RFs(k).plane(iPlane).lambdas(valid,:)];
        lambdasRun = [lambdasRun; RFs(k).plane(iPlane).lambdas_runOnly(valid)'];
    end
end

bins = unique(lambdasStim(:));
edges = [bins; 2*bins(end)] - diff(bins(1:2))/2;
figure
hold on
for m = 1:size(lambdasStim,2)
    n = histcounts(lambdasStim(:,m), edges);
    plot(1:length(bins), n, 'LineWidth', 2)
end
set(gca,'XTick',1:length(bins),'XTickLabel',num2str(bins,'%1.0e\n'))
xlabel('Lambda for RF')
ylabel('#Units')
legend(models)

bins = unique(lambdasRun);
edges = [bins; 2*bins(end)] - diff(bins(1:2))/2;
figure
hold on
n = histcounts(lambdasRun, edges);
plot(1:length(bins), n, 'k', 'LineWidth', 2)
set(gca,'XTick',1:length(bins),'XTickLabel',num2str(bins,'%1.0e\n'))
xlabel('Lambda for Running')
ylabel('#Units')
