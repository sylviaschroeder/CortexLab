%% Dataset
db_ephys_opticTract

%% Define folders
% protocolFolder = '\\ZSERVER.cortexlab.net\Data\trodes';

% protocolFolder = '\\ZSERVER.cortexlab.net\Data\expInfo';
% hardwareFolder = '\\ZSERVER.cortexlab.net\Data\expInfo';
% timelineFolder = '\\ZSERVER.cortexlab.net\Data\expInfo';
% subjectsFolder = '\\ZSERVER.cortexlab.net\Data\Subjects';

protocolFolder = '\\ZUBJECTS.cortexlab.net\Subjects';
hardwareFolder = '\\ZUBJECTS.cortexlab.net\Subjects';
subjectsFolder = '\\ZUBJECTS.cortexlab.net\Subjects';
% subjectsFolder = 'J:\Ephys';

% file location depending on PCs
folderPC = 'C:\STORAGE\OneDrive - University College London';
folderResults = fullfile(folderPC, 'Lab\RESULTS\receptiveFields\axons');

%% Parameters
lambdasStim = [0 0.0001];
lambdasRun = [0 logspace(0, 6, 7)];
RFlimits = [0 0.1];
crossFolds = 10;
rfTypes = {'ON','OFF'};
runKrnlLimits = [-5 5];
% colormap
grad = linspace(0,1,100)';
reds = [ones(100,1),grad,grad];
blues = [grad,grad,ones(100,1)];
cm = [blues; flip(reds(1:end-1,:),1)];

%% Fit RFs and get cross-validated explained variance

if ~isfolder(folderResults)
    mkdir(folderResults)
end

RFs = db;
RF_predictions = db;
for k = 1:length(db)
    fprintf('\nProcessing: %s %s\n', db(k).subject, db(k).date);
    alignDir = fullfile(subjectsFolder, db(k).subject, db(k).date, 'alignments');
    % Load spike data
    sp = loadAllKsDir(db(k).subject, db(k).date);
    
    [expNums, ~, ~, ~, ~, tl, hasTimeline] = ...
        dat.whichExpNums(db(k).subject, db(k).date);
    TLexp = expNums(hasTimeline);
    TLexp = TLexp(end);
    tl = tl{end};
    bTLtoMaster = readNPY(fullfile(alignDir, ...
        sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, db(k).probeNames{1})));
    tlTime = applyCorrection(tl.rawDAQTimestamps, bTLtoMaster);
    
    % Load stimulus information
    % data = load(fullfile(protocolFolder, subject, date, num2str(exp), ...
    %     'Protocol.mat'));
    data = load(fullfile(protocolFolder, db(k).subject, db(k).date, ...
        num2str(db(k).expNoise), sprintf('%s_%d_%s_parameters.mat', ...
        db(k).date, db(k).expNoise, db(k).subject)));
    parsNoise = data.parameters.Protocol;
    stimFile = str2func(strtok(parsNoise.xfile, '.'));
    % load myScreenInfo
    load(fullfile(hardwareFolder, db(k).subject, db(k).date, ...
        num2str(db(k).expNoise), sprintf('%s_%d_%s_hardwareInfo.mat', ...
        db(k).date, db(k).expNoise, db(k).subject)));
    myScreenInfo.windowPtr = NaN;
    stimOnTL = readNPY(fullfile(alignDir, ...
        sprintf('mpep_%d_onsets_in_timeline_%d.npy', db(k).expNoise, TLexp)));
    stimOffTL = readNPY(fullfile(alignDir, ...
        sprintf('mpep_%d_offsets_in_timeline_%d.npy', db(k).expNoise, TLexp)));
    stimTimes.onset = applyCorrection(stimOnTL, bTLtoMaster);
    stimTimes.offset = applyCorrection(stimOffTL, bTLtoMaster);
    
    % call x-file to create stimuli
    SS = stimFile(myScreenInfo, parsNoise.pars);
    stimFrames = cat(3, SS.ImageTextures{:});
    
    stimFrameDur = mean(stimTimes.offset - stimTimes.onset) / size(stimFrames,3);
    stimFrameTimes = ((1:size(stimFrames,3))-1) .* stimFrameDur;
    noiseFrameTimes = (stimFrameTimes + stimTimes.onset)';
    RFtimesInFrames = floor(RFlimits(1) / stimFrameDur) : ...
        ceil(RFlimits(2) / stimFrameDur);
    
    stimPosition = parsNoise.pars(2:5)./10;
    
    % if number of presented frames is smaller than number of parameters to
    % estimate, use only contralateral portion of stimulus
    if size(stimFrames,3)*length(stimTimes) < ...
            size(stimFrames,1)*size(stimFrames,2)*2*length(RFtimesInFrames)
        pxSize = diff(stimPosition(1:2)) / size(stimFrames,2);
        half = size(stimFrames,2) / 2;
        if strcmp(db(k).hemispheres{db(k).OTprobe}, 'right')
            half = ceil(half);
            stimFrames = stimFrames(:,1:half,:);
            stimPosition(2) = stimPosition(1) + half*pxSize;
        else
            half = floor(half)+1;
            stimFrames = stimFrames(:,half:end,:);
            stimPosition(1) = stimPosition(2) - size(stimFrames,2)*pxSize;
        end
    end        
    
    % Load and prepare data for running correlation
    rotary = double(tl.rawDAQData(:, strcmp({tl.hw.inputs.name}, ...
        'rotaryEncoder')));
    runSpeed = nonVis.getRunningSpeed_wheel(rotary, tlTime, 0);
    runTime = runSpeed.t;
    runSpeed = runSpeed.total;
    cmPerUnit = 2*pi * 8.75 / (4 * 1024);
    runSpeed = runSpeed * cmPerUnit;
    
    traces = [];
    for probe = 1:length(sp)
        if length(sp) > 1 && db(k).OTprobe~=probe
            continue
        end
        units_OT = db(k).OTunits{1};
        
        for iCell = 1:length(units_OT)
            if db(k).OTgood{1}(iCell) == 0
                continue
            end
            unitID = units_OT(iCell);
            
            % get spiketimes for this cell
            st = sp(probe).st(sp(probe).clu == unitID);
            
            % Visual noise: Find pixels that trigger response, then plot raster
            % for these pixels
            respPerFrame = histcounts(st, [noiseFrameTimes(:); ...
                noiseFrameTimes(end) + stimFrameDur]);
            traces = [traces, respPerFrame'];
        end
    end
    
    [rFields, runKernels, runWin, ev, ev_run, ev_stim, ...
        preds, preds_run, time] = ...
        whiteNoise.getReceptiveField_OnOff_smooth_running( ...
        traces, noiseFrameTimes, stimFrames, stimFrameTimes, stimTimes, ...
        RFtimesInFrames, runSpeed, runTime, runKrnlLimits, ...
        {lambdasStim, lambdasRun}, crossFolds, false);
    RFs(k).RFTimes = RFtimesInFrames * stimFrameDur;
    RFs(k).stimPosition = stimPosition;
    RFs(k).runWindow = runWin;
    RF_predictions(k).lambdasStim = lambdasStim;
    RF_predictions(k).lambdasRun = lambdasRun;
    RF_predictions(k).time = time;
    
    v = permute(mean(ev,2),[1 3 4 2]); % [neuron x lamStim x lamRun], average across cross-folds
    [maxEV,maxStimLam] = max(v,[],2);
    maxEV = permute(maxEV, [1 3 2]); % [neuron x lamRun];
    maxStimLam = permute(maxStimLam, [1 3 2]); % [neuron x lamRun];
    [maxEV, maxRunLam] = max(maxEV, [], 2); % [neuron x 1]
    indLam = sub2ind(size(maxStimLam), (1:size(maxStimLam,1))', maxRunLam);
    maxStimLam = maxStimLam(indLam); % [neuron x 1]
    
    vRun = permute(mean(ev_run,2), [1 3 4 2]); % [neuron x lamStim x lamRun], average across cross-folds
    vStim = permute(mean(ev_stim,2), [1 3 4 2]); % [neuron x lamStim x lamRun]
    inds = sub2ind(size(vRun), (1:size(vRun,1))', maxStimLam, maxRunLam);
    maxEVRun = vRun(inds); % [neuron x 1]
    maxEVStim = vStim(inds); % [neuron x 1]
    
    RF_predictions(k).explainedVariances = ev;
    RF_predictions(k).predictions = preds;
    RF_predictions(k).explainedVariances_runOnly = ev_run;
    RF_predictions(k).predictions_runOnly = preds_run;
    RF_predictions(k).explainedVariances_stimOnly = ev_stim;
    
    RFs(k).receptiveFields = rFields;
    RFs(k).runningKernels = runKernels;
    RFs(k).explainedVariances = maxEV;
    RFs(k).lambdasStim = lambdasStim(maxStimLam);
    RFs(k).lambdasRun = lambdasRun(maxRunLam);
    RFs(k).explainedVariances_runOnly = maxEVRun;
    RFs(k).explainedVariances_stimOnly = maxEVStim;
        
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
        RFs(k).date, RFs(k).expNoise);
    
    alignDir = fullfile(subjectsFolder, db(k).subject, db(k).date, 'alignments');
    % Load spike data
    sp = loadAllKsDir(db(k).subject, db(k).date);
    
    [expNums, ~, ~, ~, ~, tl, hasTimeline] = ...
        dat.whichExpNums(db(k).subject, db(k).date);
    TLexp = expNums(hasTimeline);
    TLexp = TLexp(end);
    tl = tl{end};
    bTLtoMaster = readNPY(fullfile(alignDir, ...
        sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, db(k).probeNames{1})));
    tlTime = applyCorrection(tl.rawDAQTimestamps, bTLtoMaster);
    
    % Load stimulus information
    % data = load(fullfile(protocolFolder, subject, date, num2str(exp), ...
    %     'Protocol.mat'));
    data = load(fullfile(protocolFolder, db(k).subject, db(k).date, ...
        num2str(db(k).expNoise), sprintf('%s_%d_%s_parameters.mat', ...
        db(k).date, db(k).expNoise, db(k).subject)));
    parsNoise = data.parameters.Protocol;
    stimFile = str2func(strtok(parsNoise.xfile, '.'));
    % load myScreenInfo
    load(fullfile(hardwareFolder, db(k).subject, db(k).date, ...
        num2str(db(k).expNoise), sprintf('%s_%d_%s_hardwareInfo.mat', ...
        db(k).date, db(k).expNoise, db(k).subject)));
    myScreenInfo.windowPtr = NaN;
    stimOnTL = readNPY(fullfile(alignDir, ...
        sprintf('mpep_%d_onsets_in_timeline_%d.npy', db(k).expNoise, TLexp)));
    stimOffTL = readNPY(fullfile(alignDir, ...
        sprintf('mpep_%d_offsets_in_timeline_%d.npy', db(k).expNoise, TLexp)));
    stimTimes.onset = applyCorrection(stimOnTL, bTLtoMaster);
    stimTimes.offset = applyCorrection(stimOffTL, bTLtoMaster);
    
    % call x-file to create stimuli
    SS = stimFile(myScreenInfo, parsNoise.pars);
    stimFrames = cat(3, SS.ImageTextures{:});
    
    stimFrameDur = mean(stimTimes.offset - stimTimes.onset) / size(stimFrames,3);
    stimFrameTimes = ((1:size(stimFrames,3))-1) .* stimFrameDur;
    noiseFrameTimes = (stimFrameTimes + stimTimes.onset)';
    RFtimesInFrames = floor(RFlimits(1) / stimFrameDur) : ...
        ceil(RFlimits(2) / stimFrameDur);
    
    stimPosition = parsNoise.pars(2:5)./10;
    
    % if number of presented frames is smaller than number of parameters to
    % estimate, use only contralateral portion of stimulus
    if size(stimFrames,3)*length(stimTimes) < ...
            size(stimFrames,1)*size(stimFrames,2)*2*length(RFtimesInFrames)
        pxSize = diff(stimPosition(1:2)) / size(stimFrames,2);
        half = size(stimFrames,2) / 2;
        if strcmp(db(k).hemispheres{db(k).OTprobe}, 'right')
            half = ceil(half);
            stimFrames = stimFrames(:,1:half,:);
            stimPosition(2) = stimPosition(1) + half*pxSize;
        else
            half = floor(half)+1;
            stimFrames = stimFrames(:,half:end,:);
            stimPosition(1) = stimPosition(2) - size(stimFrames,2)*pxSize;
        end
    end
    
    % Load and prepare data for running correlation
    rotary = double(tl.rawDAQData(:, strcmp({tl.hw.inputs.name}, ...
        'rotaryEncoder')));
    runSpeed = nonVis.getRunningSpeed_wheel(rotary, tlTime, 0);
    runTime = runSpeed.t;
    runSpeed = runSpeed.total;
    cmPerUnit = 2*pi * 8.75 / (4 * 1024);
    runSpeed = runSpeed * cmPerUnit;
    
    traces = [];
    for probe = 1:length(sp)
        if length(sp) > 1 && db(k).OTprobe~=probe
            continue
        end
        units_OT = db(k).OTunits{1};
        
        for iCell = 1:length(units_OT)
            if db(k).OTgood{1}(iCell) == 0
                continue
            end
            unitID = units_OT(iCell);
            
            % get spiketimes for this cell
            st = sp(probe).st(sp(probe).clu == unitID);
            
            % Visual noise: Find pixels that trigger response, then plot raster
            % for these pixels
            respPerFrame = histcounts(st, [noiseFrameTimes(:); ...
                noiseFrameTimes(end) + stimFrameDur]);
            traces = [traces, respPerFrame'];
        end
    end
    
    [ev, ev_shift] = ...
        whiteNoise.receptiveFieldShiftTest_OnOff( ...
        traces, noiseFrameTimes, stimFrames, stimFrameTimes, stimTimes, ...
        RFtimesInFrames, runSpeed, runTime, ...
        RFs(k).runningKernels, RFs(k).runWindow, ...
        RFs(k).receptiveFields, RFs(k).lambdasStim, 500);
    pvals = sum(ev_shift > ev, 2) ./ size(ev_shift,2);
    pvals(isnan(RFs(k).explainedVariances)) = NaN;
    RFs(k).pVal_RFonly = pvals;

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
        RFs(k).date, RFs(k).expNoise);
    
    alignDir = fullfile(subjectsFolder, db(k).subject, db(k).date, 'alignments');
    % Load spike data
    sp = loadAllKsDir(db(k).subject, db(k).date);
    
    [expNums, ~, ~, ~, ~, ~, hasTimeline] = ...
        dat.whichExpNums(db(k).subject, db(k).date);
    TLexp = expNums(hasTimeline);
    TLexp = TLexp(end);
    bTLtoMaster = readNPY(fullfile(alignDir, ...
        sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, db(k).probeNames{1})));
    
    % Load stimulus information
    % data = load(fullfile(protocolFolder, subject, date, num2str(exp), ...
    %     'Protocol.mat'));
    data = load(fullfile(protocolFolder, db(k).subject, db(k).date, ...
        num2str(db(k).expNoise), sprintf('%s_%d_%s_parameters.mat', ...
        db(k).date, db(k).expNoise, db(k).subject)));
    parsNoise = data.parameters.Protocol;
    stimFile = str2func(strtok(parsNoise.xfile, '.'));
    % load myScreenInfo
    load(fullfile(hardwareFolder, db(k).subject, db(k).date, ...
        num2str(db(k).expNoise), sprintf('%s_%d_%s_hardwareInfo.mat', ...
        db(k).date, db(k).expNoise, db(k).subject)));
    myScreenInfo.windowPtr = NaN;
    stimOnTL = readNPY(fullfile(alignDir, ...
        sprintf('mpep_%d_onsets_in_timeline_%d.npy', db(k).expNoise, TLexp)));
    stimOffTL = readNPY(fullfile(alignDir, ...
        sprintf('mpep_%d_offsets_in_timeline_%d.npy', db(k).expNoise, TLexp)));
    stimTimes.onset = applyCorrection(stimOnTL, bTLtoMaster);
    stimTimes.offset = applyCorrection(stimOffTL, bTLtoMaster);
    
    % call x-file to create stimuli
    SS = stimFile(myScreenInfo, parsNoise.pars);
    stimFrames = cat(3, SS.ImageTextures{:});
    
    stimFrameDur = mean(stimTimes.offset - stimTimes.onset) / size(stimFrames,3);
    stimFrameTimes = ((1:size(stimFrames,3))-1) .* stimFrameDur;
    traceTimes = (stimFrameTimes + stimTimes.onset)';
    
    traces = [];
    IDs = [];
    for probe = 1:length(sp)
        if length(sp) > 1 && RFs(k).OTprobe~=probe
            continue
        end
        units_OT = db(k).OTunits{1};
        
        for iCell = 1:length(units_OT)
            if RFs(k).OTgood{1}(iCell) == 0
                continue
            end
            IDs(end+1,1) = units_OT(iCell);
            
            % get spiketimes for this cell
            st = sp(probe).st(sp(probe).clu == IDs(end));
            
            % Visual noise: Find pixels that trigger response, then plot raster
            % for these pixels
            respPerFrame = histcounts(st, [traceTimes(:); ...
                traceTimes(end) + stimFrameDur]);
            traces = [traces, respPerFrame'];
        end
    end
    
    traces = (traces - nanmean(traces)) ./ nanstd(traces);
    
    rFields = RFs(k).receptiveFields;
    runKernels = RFs(k).runningKernels;
    ev = RFs(k).explainedVariances;
    evRun = RFs(k).explainedVariances_runOnly;
    evStim = RFs(k).explainedVariances_stimOnly;
    pValues = RFs(k).pVal_RFonly;
    runWin = RFs(k).runWindow;
    time = RF_preds(k).time;
    preds = RF_preds(k).predictions;
    predsRun = RF_preds(k).predictions_runOnly;
    
    fPlots = fullfile(folderResults, 'plots_kernels', ...
        sprintf('%s_%s', RFs(k).subject, RFs(k).date));
    if ~isfolder(fPlots)
        mkdir(fPlots)
    end
    for iCell = 1:length(RFs(k).explainedVariances)
        if isnan(RFs(k).explainedVariances(iCell))
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
        annotation('textbox', [0 .95 1 .03], 'String', sprintf('Neuron %d', IDs(iCell)), ...
            'FontSize', 14, 'FontWeight', 'bold', 'LineStyle', 'none', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(fullfile(fPlots, sprintf('Neuron%03d.jpg', IDs(iCell))), ...
            '-djpeg','-r0')
        close gcf
    end
end