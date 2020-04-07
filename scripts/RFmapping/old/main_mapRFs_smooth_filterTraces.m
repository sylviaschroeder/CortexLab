%% Load database
db_sparseNoise
% db_boutons_sparseNoise

%% Folders
% loading data
folderROIData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\InfoStructs';
% saving data
folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\receptiveFields\SC neurons';
% folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\boutons\receptiveFields';

%% Parameters
% lambda = [0.05 0.1 0.5 1 5 10];
% lambda = [0.5 1 5 10 50 100];
% lambda = [0.5 1 5 10 50 100 500 1000];
% lambda = [0.1 0.5 1 5 10 50 100];
lambdas = logspace(-2, 8, 11);
RFlimits = [0.1 0.6];
crossFolds = 10;
models = {'linear','absolute','white','black'};

% parameters for high-pass filtering neural trace
sigma = 3; % in s; sigma of Gaussian used to filter traces

% colormap
grad = linspace(0,1,40)';
reds = [ones(40,1),grad,grad];
blues = [grad,grad,ones(40,1)];
cm = [[1 1 1]; blues; flip(reds(1:end-1,:),1)];

%% Fit RFs

RFs = db;
RF_fitting = db;
gDev = gpuDevice;
reset(gDev);
for k = 1:length(db)
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
        end
        
        %get calcium traces
        tr = meta.F_final;
        % set data of duplicate neurons or unhealthy neurons to NaN
        tr(:, meta.ROI.isDuplicate==1 | meta.ROI.isSwitchOn == 1) = NaN;
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
                tr_int(:,n) = interp1(t, tr, traceTimes, 'pchip');
            end
            traces = [traces, tr_int];
            numNeurons(iPlane) = size(tr,2);
        end
    end
    % high-pass filter traces by convolving them with smooth Gaussian and
    % subtracting the result
    sigInFr = ceil(sigma / traceBin);
    win = normpdf(-3*sigInFr : 3*sigInFr, 0, sigInFr);
    len = size(traces,1);
    convMat = convmtx(win, len);
    filtered = convMat' * traces;
    filtered = filtered(3*sigInFr + (1:len),:);
    tracesFil = traces - filtered;
    
    [rFields, ev, tracesTest, preds] = ...
        whiteNoise.getReceptiveField_smooth( ...
        traces, traceTimes, stimFrames, stimFrameTimes, stimTimes, ...
        RFtimesInFrames, lambdas, crossFolds, models);
%     [rFields, ev, tracesTest, preds] = ...
%         whiteNoise.getReceptiveField_smooth( ...
%         traces, traceTimes, stimFrames, stimFrameTimes, stimTimes, ...
%         RFtimesInFrames, lambda, crossFolds, models);
%     [rFields, intercepts, ev, tracesTest, preds] = ...
%         whiteNoise.getReceptiveField_smooth_intercept( ...
%         traces, traceTimes, stimFrames, stimFrameTimes, stimTimes, ...
%         RFtimesInFrames, lambda, crossFolds, models);
    RFs(k).RFTimes = RFtimesInFrames * stimTimeStep;
    RFs(k).stimPosition = stimPosition;
    RFs(k).models = models;
    RF_fitting(k).RFTimes = RFtimesInFrames * stimTimeStep;
    RF_fitting(k).stimPosition = stimPosition;
    RF_fitting(k).models = models;
    RF_fitting(k).lambdas = lambdas;

    n = 0;
    for iPlane = 1:length(db(k).planes)
        ind = (1:numNeurons(iPlane))+n;
        RF_fitting(k).plane(iPlane).receptiveFields = rFields(:,:,:,ind,:,:,:);
        RF_fitting(k).plane(iPlane).explainedVariances = ev(ind,:,:,:);
        RF_fitting(k).plane(iPlane).testTraces = tracesTest(:,ind,:);
        RF_fitting(k).plane(iPlane).predictions = preds(:,ind,:,:,:);
        % in RFs only store neurons' best RFs for each model
        v = squeeze(mean(ev(ind,:,:,:),2)); % [neuron x lambda x model], average across cross folds
        [maxEV,maxLam] = max(v,[],2);
        rfs = NaN(size(rFields,1), size(rFields,2), size(rFields,3), ...
            numNeurons(iPlane), length(models));
        for m = 1:length(models)
            for iCell = 1:length(ind)
                rfs(:,:,:,iCell,m) = squeeze(mean(...
                    rFields(:,:,:,ind(iCell),:,maxLam(iCell,1,m),m), 5));
            end
        end
        for m = 1:length(models)
            RFs(k).plane(iPlane).receptiveFields = rfs;
            RFs(k).plane(iPlane).explainedVariances = squeeze(maxEV);
            RFs(k).plane(iPlane).lambdas = lambdas(squeeze(maxLam));
        end
        save(fullfile(folderResults, 'receptiveFields.mat'), 'RFs')
        save(fullfile(folderResults, 'receptiveFields_fitting.mat'), 'RF_fitting')
        
        % plot RFs
        fPlots = fullfile(folderResults, 'plots_kernels', ...
            sprintf('%s_%s', db(k).subject, db(k).date), ...
            sprintf('Plane%02d',db(k).planes(iPlane)));
        if ~isfolder(fPlots)
            mkdir(fPlots)
        end
        for iCell = 1:length(ind)
            if all(isnan(traces(:,ind(iCell))))
                continue
            end
            figure('Position',[2230 2 1600 1110])
            clim = max(abs(reshape(rfs(:,:,:,iCell,:),[],1)));
            cols = size(rfs,2);
            subplot(length(models)+1,1,length(models)+1)
            tr_test = reshape(tracesTest(:,ind(iCell),:),[],1);
            plot((1:length(tr_test)) .* stimTimeStep, tr_test, 'k')
            hold on
            for m = 1:length(models)
                r = [rfs(:,:,:,iCell,m),NaN(size(rfs,1),1,size(rfs,3))];
                r = reshape(r,size(r,1),[]);
                r(:,end) = [];
                subplot(length(models)+1,1,m)
                imagesc(r,[-clim clim])
                axis image
                set(gca,'box','off','XTick',[],'YTick',[])
                if m == length(models)
                    set(gca,'XTick',ceil(cols/2):(cols+1):(cols+1)*size(rfs,3), ...
                        'XTickLabel',sprintf('%.2f\n',RF_fitting(k).RFTimes))
                    xlabel('Time from response (in s)')
                end
                ylabel(models{m})
                title(sprintf('explained variance: %.2f%%', maxEV(iCell,1,m)*100))
                
                subplot(length(models)+1,1,length(models)+1)
                pred = reshape(preds(:,ind(iCell),:,maxLam(iCell,1,m),m),[],1);
                plot((1:length(tr_test)) .* stimTimeStep, pred)
            end
            colormap(cm)
            axis tight
            set(gca,'box','off')
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
    fprintf('\n')
end