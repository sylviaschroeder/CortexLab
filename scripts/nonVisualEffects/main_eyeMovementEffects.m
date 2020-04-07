folderROIData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\InfoStructs';
folderKernels = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\nonVisualEffects\modelGratingResp\';
folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\nonVisualEffects\eyeMovements';

%% Parameters
timeWin = [-1 2];
numShifts = 500;

%% Load data
db_driftingGratings_blank

% data = load(fullfile(folderKernels, 'kernelFit', 'results.mat'));
% results = data.results;

%% Saccade-triggered average of calcium during gray screens

if ~isfolder(fullfile(folderResults, 'grayScreen'))
    mkdir(fullfile(folderResults, 'grayScreen'))
end

for iExp = 1:length(db)
    fprintf('Dataset %d (of %d)\n', iExp, length(db))
    if isempty(db(iExp).expGrayScreen)
        continue
    end
    folder = fullfile(folderROIData, db(iExp).subject, ...
        db(iExp).date, num2str(db(iExp).expGrayScreen));
    file = [sprintf('%s_%d_%s', db(iExp).date, ...
        db(iExp).expGrayScreen, db(iExp).subject) '_2P_plane%03d_ROI.mat'];
    % load meta
    data = load(fullfile(folder, sprintf(file,db(iExp).planes(1))));
    meta = data.meta;
    meta.folderTL = strrep(meta.folderTL, 'ioo.ucl.ac.uk', ...
        'cortexlab.net');
    frameTimes = ppbox.getFrameTimes(meta);
    dt = median(diff(frameTimes));
    window = floor(timeWin(1)/dt) : ceil(timeWin(2)/dt);
    win_t = window .* dt;
    trigTimes = win_t;
    % load pupil data
    [pupilData, pupilTime] = nonVis.loadPupilData(meta);
    if isempty(pupilData)
        continue
    end
    pupilTime(length(pupilData.x)+1:end) = [];
    % get saccade times
    x = pupilData.x;
    y = pupilData.y;
    ind = isnan(x) | isnan(y) | pupilData.blink | ~pupilData.goodFit;
    x(ind) = NaN;
    y(ind) = NaN;
    dt = median(diff(pupilTime));
    N = min(9, ceil(floor(0.4/dt)/2)*2-1);
    x = medfilt1(x,N);
    y = medfilt1(y,N);
    saccadeTimes = eye.findSaccades(x,y,ceil(0.2/dt),0.75,1);
    annotation('textbox', [0 .95 1 .05], 'String', ...
        sprintf('%s %s, gray screen (exp %d)', db(iExp).subject, ...
        db(iExp).date, db(iExp).expGrayScreen), 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'interpreter', 'none', 'FontWeight', 'bold')
    figFolder = fullfile(folderResults, 'grayScreen', sprintf('%s %s', ...
        db(iExp).subject, db(iExp).date));
    if ~isfolder(figFolder)
        mkdir(figFolder)
    end 
    savefig(gcf, fullfile(figFolder, 'saccades.fig'))
    close gcf
    
    saccadeTimes = pupilTime(saccadeTimes);
    % get triggered traces for each neuron
    fprintf('  Plane (of %d):', length(db(iExp).planes))
    for iPlane = 1:length(db(iExp).planes)
        fprintf(' %d', iPlane)
        if iPlane > 1
            % load meta
            data = load(fullfile(folder, sprintf(file,db(iExp).planes(iPlane))));
            meta = data.meta;
            meta.folderTL = strrep(meta.folderTL, 'ioo.ucl.ac.uk', ...
                'cortexlab.net');
            frameTimes = ppbox.getFrameTimes(meta);
        end
        [~,saccFr] = find(cumsum(frameTimes - saccadeTimes( ...
            saccadeTimes>frameTimes(1) & saccadeTimes<frameTimes(end)) ...
            > 0, 2) == 1);
        pseudoSaccs = saccFr + randi(length(frameTimes), 1, numShifts);
        pseudoSaccs = mod(pseudoSaccs-1, length(frameTimes)) + 1;
        
        % only consider ROIs that are unique
        if isfield(meta.ROI, 'isDuplicate')
            ids = find(meta.ROI.isDuplicate == 0 & meta.ROI.isSwitchOn == 0 & ...
                ~all(isnan(meta.F_final),1)');
        else
            ids = find(~all(isnan(meta.F_final),1)');
        end
        calciumTraces = meta.F_final(:,ids);
        
        fr = permute(saccFr, [2 3 1]) + window';
        ind = fr < 1 | fr > size(calciumTraces,1);
        fr(ind) = 1;
        fr = repmat(fr, 1, length(ids), 1, 1) + ...
            (0:length(ids)-1).* size(calciumTraces,1);
        tt = calciumTraces(fr);
        tt(repmat(ind,1,length(ids),1,1)) = NaN;
        trigTraces = tt;
        
        ps = permute(pseudoSaccs, [3 4 1 2]);
        ps = ps + window';
        ind = ps < 1 | ps > size(calciumTraces,1);
        ps(ind) = 1;
        ps = repmat(ps, 1, length(ids), 1, 1) + ...
            (0:length(ids)-1).*size(calciumTraces,1);
        tt = calciumTraces(ps);
        tt(repmat(ind,1,length(ids),1,1)) = NaN;
        pseudoTraces = tt;
        
        ind0 = find(trigTimes < 0, 3, 'last');
        tr = trigTraces - nanmean(trigTraces(ind0,:,:),1);
        tr_avg = nanmean(tr,3); % [time x neuron]
        nTr = pseudoTraces - nanmean(pseudoTraces(ind0,:,:,:),1);
        nTr_avg = nanmean(nTr,3);
        % calculate significance level for multiple comparisons: test all time
        % points in the interval [0 1.1)
        testInds = trigTimes>=0 & trigTimes < 1.1;
        numTests = sum(testInds);
        alpha = 0.05 / numTests;
        confInt = squeeze(prctile(nTr_avg, alpha/2*100.*[1 -1] + [0 100], 4));
        sgnf = tr_avg < confInt(:,:,1) | tr_avg > confInt(:,:,2);
        
        sgnfNeurons = find(any(sgnf(testInds,:),1));
        for iCell = 1:length(sgnfNeurons)
            figure
            hold on
            fill(trigTimes([1:end,end:-1:1]), [confInt(:,sgnfNeurons(iCell),1)' ...
                flip(confInt(:,sgnfNeurons(iCell),2))'], 'k', 'EdgeColor', 'none', 'FaceColor', [1 1 1].*0.9)
            %         plot(trigTimes{iExp},squeeze(tr(:,sgnfNeurons(iCell),:)), 'Color', [1 1 1].*0.8)
            plot(trigTimes,tr_avg(:,sgnfNeurons(iCell)), 'k', 'LineWidth', 2)
            plot(trigTimes(sgnf(:,sgnfNeurons(iCell))&testInds'), ...
                tr_avg(sgnf(:,sgnfNeurons(iCell))&testInds',sgnfNeurons(iCell)), ...
                'ro', 'MarkerFaceColor', 'r')
            axis tight
            xlabel('Time from saccade (in s)')
            ylabel('\DeltaF/F')
            title(sprintf('%s %s plane %d neuron %d', db(iExp).subject, ...
                db(iExp).date, db(iExp).planes(iPlane), ...
                ids(sgnfNeurons(iCell))), 'Interpreter', 'none')
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(fullfile(figFolder, sprintf('trace_plane%02d_neuron%03d.jpg', ...
                db(iExp).planes(iPlane), ids(sgnfNeurons(iCell)))), ...
                '-djpeg','-r0')
            close(fig)
        end
    end
    fprintf('\n')
end

%% Plot saccade triggered traces for each neuron
for iExp = 1:length(traces)
    if isempty(traces{iExp})
        continue
    end
    ind0 = find(trigTimes{iExp} < 0, 3, 'last');
    tr = traces{iExp} - nanmean(traces{iExp}(ind0,:,:),1);
    tr_avg = nanmean(tr,3); % [time x neuron]
    nTr = nullTraces{iExp} - nanmean(nullTraces{iExp}(ind0,:,:,:),1);
    nTr_avg = nanmean(nTr,3);
    % calculate significance level for multiple comparisons: test all time
    % points in the interval [0 1.1)
    testInds = trigTimes{iExp}>=0 & trigTimes{iExp} < 1.1;
    numTests = sum(testInds);
    alpha = 0.05 / numTests;
    confInt = squeeze(prctile(nTr_avg, alpha/2*100.*[1 -1] + [0 100], 4));
    sgnf = tr_avg < confInt(:,:,1) | tr_avg > confInt(:,:,2);
    
    folder = fullfile(folderResults, 'grayScreen', sprintf('%s %s', ...
        db(iExp).subject, db(iExp).date));
    
    sgnfNeurons = find(any(sgnf(testInds,:),1));
    figure
    plot(trigTimes{iExp}, tr_avg(:,sgnfNeurons), 'k')
    axis tight
    set(gca, 'box', 'off')
    xlabel('Time (in s)')
    ylabel('\DeltaF/F')
    title(sprintf('Significant saccade triggered traces (n = %d)',length(sgnfNeurons)))
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(fullfile(folder, 'allSaccadeTriggeredTraces.jpg'), ...
        '-djpeg','-r0')
    close fig
    for iCell = 1:length(sgnfNeurons)
        figure
        hold on
        fill(trigTimes{iExp}([1:end,end:-1:1]), [confInt(:,sgnfNeurons(iCell),1)' ...
            flip(confInt(:,sgnfNeurons(iCell),2))'], 'k', 'EdgeColor', 'none', 'FaceColor', [1 1 1].*0.9)
%         plot(trigTimes{iExp},squeeze(tr(:,sgnfNeurons(iCell),:)), 'Color', [1 1 1].*0.8)
        plot(trigTimes{iExp},tr_avg(:,sgnfNeurons(iCell)), 'k', 'LineWidth', 2)
        plot(trigTimes{iExp}(sgnf(:,sgnfNeurons(iCell))&testInds'), ...
            tr_avg(sgnf(:,sgnfNeurons(iCell))&testInds',sgnfNeurons(iCell)), ...
            'ro', 'MarkerFaceColor', 'r')
        axis tight
        xlabel('Time from saccade (in s)')
        ylabel('\DeltaF/F')
        title(sprintf('%s %s plane %d neuron %d', db(iExp).subject, ...
            db(iExp).date, cellIDs{iExp}(sgnfNeurons(iCell),1), ...
            cellIDs{iExp}(sgnfNeurons(iCell),2)), ...
            'Interpreter', 'none')
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(fullfile(folder, sprintf('trace_plane%02d_neuron%03d.jpg', ...
            cellIDs{iExp}(sgnfNeurons(iCell),1), cellIDs{iExp}(sgnfNeurons(iCell),2))), ...
            '-djpeg','-r0')
        close fig
    end
end
