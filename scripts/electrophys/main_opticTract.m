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
timelineFolder = '\\ZUBJECTS.cortexlab.net\Subjects';
subjectsFolder = '\\ZUBJECTS.cortexlab.net\Subjects';
% subjectsFolder = 'J:\Ephys';

plotFolder = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\OpticTract';
waveformFolder = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\OpticTract\Waveforms';
rawSpikeFolder = 'J:\Ephys';
% rawSpikeFolder = '\\zubjects.cortexlab.net\Subjects';

%% Parameters
% for RFs
paramsNoisePsth.makePlots = false;
paramsNoisePsth.useSVD = true;
paramsNoisePsth.countWindow = [-0.05 0.18];
paramsNoisePsth.binSize = 0.001;
paramsNoiseRF.makePlots = false;
paramsNoiseRF.useSVD = true;
paramsNoiseRF.countWindow = [-0.05 0.18];
paramsNoiseRF.binSize = 0.01;

RFtypes = {'Absolute'}; %, 'White', 'Black', 'Linear'};
grad = linspace(0,1,40)';
reds = [ones(40,1),grad,grad];
blues = [grad,grad,ones(40,1)];
cm = [blues; flip(reds(1:end-1,:),1)];
tmp = lines(3);
colors = [tmp(1,:); tmp(3,:); 0 0 0];

nPseudo = 20;

% for gratings
binSizeGrating = 0.001;
winLength = .04; % in s
winHamm = hamming(winLength / binSizeGrating);
winHamm = winHamm ./ sum(winHamm);
thresholdFactorGratings = 3; % to detect response delay (for gratings); x std of derivative
thresholdFactorNoise = 3;
minThresh = 3;

% for flicker stim
binSizeFlicker = 0.005;

%% (2) Make plot of features for each potential unit in optic tract

for k = 1:length(db)
    fprintf('\nProcessing: %s %s\n', db(k).subject, db(k).date);
    alignDir = fullfile(subjectsFolder, db(k).subject, db(k).date, 'alignments');
    %% Load spike data
    sp = loadAllKsDir(db(k).subject, db(k).date);
    
    [expNums, ~, ~, ~, ~, ~, hasTimeline] = ...
        dat.whichExpNums(db(k).subject, db(k).date);
    TLexp = expNums(hasTimeline);
    TLexp = TLexp(end);
    bTLtoMaster = readNPY(fullfile(alignDir, ...
        sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, db(k).probeNames{1})));
    
    %% Load and prepare data for flickering screen analysis
    if isfield(db, 'expFlicker') && ~isempty(db(k).expFlicker)
        data = load(fullfile(protocolFolder, db(k).subject, db(k).date, ...
            num2str(db(k).expFlicker), sprintf('%s_%d_%s_parameters.mat', ...
            db(k).date, db(k).expFlicker, db(k).subject)));
        parsFlicker = data.parameters.Protocol;
        data = load(fullfile(protocolFolder, db(k).subject, db(k).date, ...
            num2str(db(k).expFlicker), sprintf('%s_%d_%s_hardwareInfo.mat', ...
            db(k).date, db(k).expFlicker, db(k).subject)));
        monitorRefreshRate = data.myScreenInfo.FrameRate;
        flickerDurs = parsFlicker.pars(strcmp(parsFlicker.parnames, ...
            'nfr'),:) / monitorRefreshRate;
        
        stimOnTL = readNPY(fullfile(alignDir, ...
            sprintf('mpep_%d_onsets_in_timeline_%d.npy', db(k).expFlicker, TLexp)));
        stimOffTL = readNPY(fullfile(alignDir, ...
            sprintf('mpep_%d_offsets_in_timeline_%d.npy', db(k).expFlicker, TLexp)));
        flickerWhite = applyCorrection(stimOnTL, bTLtoMaster);
        flickerBlack = applyCorrection(stimOffTL, bTLtoMaster);
        flickerTimes = cell(size(parsFlicker.seqnums));
        if length(flickerWhite) == size(parsFlicker.pars,2)*parsFlicker.nrepeats % only stim onset was recorded, not every flick
            for stim = 1:size(flickerTimes,1)
                for trial = 1:size(flickerTimes,2)
                    ind = parsFlicker.seqnums(stim,trial);
                    flickerTimes{stim,trial} = (flickerWhite(ind) : ...
                        flickerDurs(stim) : flickerBlack(ind))' + ...
                        [0 flickerDurs(stim)/2];
                    if flickerBlack(ind) < flickerTimes{stim,trial}(end,2)
                        flickerTimes{stim,trial}(end,:) = [];
                    end
                end
            end
            flickerAll = [flickerWhite(1) flickerBlack(end)];
        else % each flicker was recorded
            flickerAll = reshape([flickerWhite, flickerBlack]', [], 1);
            longestFlicker = max(parsFlicker.pars( ...
                strcmp(parsFlicker.parnames, 'nfr'),:)) / monitorRefreshRate;
            durs = diff(flickerAll);
            flickerStimStartInds = [1; find(durs > ...
                longestFlicker * 1.5) + 1; length(flickerAll)+1];
            for stim = 1:size(flickerTimes,1)
                for trial = 1:size(flickerTimes,2)
                    ind = parsFlicker.seqnums(stim,trial);
                    flickerTimes{stim,trial} = reshape(flickerAll( ...
                        flickerStimStartInds(ind) : ...
                        flickerStimStartInds(ind+1)-1), 2, [])';
                end
            end
        end
    else
        flickerTimes = [];
    end
    
    %% Load and prepare data for receptive fields
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
    
    noiseOn = applyCorrection(stimOnTL, bTLtoMaster);
    
    % call x-file to create stimuli
    SS = stimFile(myScreenInfo, parsNoise.pars);
    stimFrames = cat(3, SS.ImageTextures{:});
    
    framesPerImage = parsNoise.pars(6,1);
    frameTimes = (0 : size(stimFrames, 3)-1) * framesPerImage / myScreenInfo.FrameRate;
    noiseFrameTimes = (frameTimes + noiseOn)';
    
    noisePosition = parsNoise.pars(2:5)./10;
    
    %% Load and prepare data for tuning curves
    data = load(fullfile(protocolFolder, db(k).subject, db(k).date, ...
        num2str(db(k).expOri), sprintf('%s_%d_%s_parameters.mat', ...
        db(k).date, db(k).expOri, db(k).subject)));
    parsGratings = data.parameters.Protocol;
    
    stimOnTL = readNPY(fullfile(alignDir, ...
        sprintf('mpep_%d_onsets_in_timeline_%d.npy', db(k).expOri, TLexp)));
    stimOffTL = readNPY(fullfile(alignDir, ...
        sprintf('mpep_%d_offsets_in_timeline_%d.npy', db(k).expOri, TLexp)));
    stimOn = applyCorrection(stimOnTL, bTLtoMaster);
    stimOff = applyCorrection(stimOffTL, bTLtoMaster);
    stimDur = mean(stimOff - stimOn);
    window = [-0.2 stimDur+0.2];
    binBorders = window(1) : binSizeGrating : window(2);
    binsGrating = binBorders(1:end-1) + binSizeGrating/2;
    stimIDs = repmat((1:parsGratings.npfilestimuli)',parsGratings.nrepeats,1);
    [~,order] = sort(parsGratings.seqnums(:));
    stimSeq = stimIDs(order);
    stimBins = binsGrating>0 & binsGrating<stimDur;
    blank = parsGratings.pars(15,:) == 0;
    directions = parsGratings.pars(6,:);
    directions(blank) = [];
    
    if ~isempty(db(k).expOriTropic)
        stimOnTL = readNPY(fullfile(alignDir, ...
            sprintf('mpep_%d_onsets_in_timeline_%d.npy', db(k).expOriTropic, TLexp)));
        stimOffTL = readNPY(fullfile(alignDir, ...
            sprintf('mpep_%d_offsets_in_timeline_%d.npy', db(k).expOriTropic, TLexp)));
        stimOnTropic = applyCorrection(stimOnTL, bTLtoMaster);
        stimOffTropic = applyCorrection(stimOffTL, bTLtoMaster);
    end
    
    %% Make plot for each neuron potentially in optic tract
    plotRows = 9;
    plotCols = 6;
    plotsTotal = plotRows * plotCols;
    
    binW = .02;
    winHamLarge = hamming(10 / binW);
    winHamLarge = winHamLarge ./ sum(winHamLarge);
    
    numShifts = 500;
    
    params_fft.tapers = [3 5];
    params_fft.pad = 0;
    params_fft.Fs = 1 / binSizeGrating;
    params_fft.fpass = [0 floor(params_fft.Fs/2)];
    tf_gratings = parsGratings.pars(strcmp('tf',parsGratings.parnames),1) / 10;
%     params_fft.fpass = [0 tf_gratings * 2];
    params_fft.trialave = 1;
    t = 0:1/params_fft.Fs:stimDur;
    N=length(t); % number of points in grid for dpss
    nfft=max(2^(nextpow2(N)+params_fft.pad),N); % number of points in fft of prolates
    [f,findx]=getfgrid(params_fft.Fs,nfft,params_fft.fpass); % get frequency grid for evaluation
    tap=dpsschk(params_fft.tapers,N,params_fft.Fs); % check tapers
    
    cols = lines(4);
    if ~isempty(flickerTimes)
        colsFlicker = cool(length(flickerDurs));
    end
    for probe = 1:length(sp)
        OTdepth = db(k).OTdepth{probe};
        if isempty(OTdepth) % probe did not pass OT
            continue
        end
        rawF = fullfile(rawSpikeFolder, db(k).subject, ...
            db(k).date, ['ephys_' sp(probe).name], 'sorting');

        if isfolder(rawF)
            [~, ~, spikeDepths] = ksDriftmap(rawF);
        end
    
        folder = fullfile(plotFolder, 'plots_basicFeatures_kernelRF', ...
            [db(k).subject '_' db(k).date '_' sp(probe).name '_KS2_Th_8_4']);
%         folder = fullfile(plotFolder, 'plots_basicFeatures_kernelRF', ...
%             [db(k).subject '_' db(k).date '_' sp(probe).name '_new']);
        if ~isfolder(folder)
            mkdir(folder)
        end
        
        units = unique(sp(probe).clu);
        units_good = sp(probe).cids(sp(probe).cgs == 2);
        units_mua = sp(probe).cids(sp(probe).cgs == 1);
        units_OT = db(k).OTunits;
        if iscell(units_OT)
            units_OT = units_OT{probe};
        end
        templates = findTempForEachClu(sp(probe).clu, sp(probe).spikeTemplates);
        
        t_total = (floor(sp(probe).st(1)/binW):ceil(sp(probe).st(end)/binW)).*binW;
        
        ind_dark = t_total > db(k).darkTime(1) & t_total < db(k).darkTime(2);
        
        for iCell = 1:length(units)
            unitID = units(iCell);
            OT_ID = find(units_OT == unitID);
            
            % get spiketimes for this cell
            st = sp(probe).st(sp(probe).clu == unitID);
            sa = sp(probe).spikeAmps(sp(probe).clu == unitID);
            
            sd = spikeDepths(sp(probe).clu == unitID);
            depth = nanmean(sd);
            if depth < OTdepth(1) || depth > OTdepth(2) % unit not in OT
                continue
            end
            
            % complete time course of spike rate
            spikeRate = hist(st, t_total) ./ binW;
            spikeRate = conv([ones(1,ceil((length(winHamLarge)-1)/2)) .* ...
                mean(spikeRate(1:100)), spikeRate, ...
                ones(1, floor((length(winHamLarge)-1)/2)) .* ...
                mean(spikeRate(end-99:end))], winHamLarge, 'valid');
            % complete time course of spike amplitude
            duplicates = [];
                % remove duplicates
            uni = unique(st);
            n = hist(st,uni);
            if any(n>1)
                dup = find(n>1);
                for d = 1:length(dup)
                    ind = find(st == uni(dup(d)));
                    duplicates = [duplicates; ind(2:end)];
                end
            end
            valid = true(size(st));
            valid(duplicates) = false;
            spikeAmp = interp1(st(valid), sa(valid), t_total, 'pchip');
            spikeAmp = conv([ones(1,ceil((length(winHamLarge)-1)/2)) .* ...
                mean(spikeAmp(1:100)), spikeAmp, ...
                ones(1, floor((length(winHamLarge)-1)/2)) .* ...
                mean(spikeAmp(end-99:end))], winHamLarge, 'valid');
            % corr of spikeRate and spikeAmp
            if isempty(OT_ID) || isempty(db(k).OTtimes{probe}{OT_ID})
                valid = true(1, length(t_total));
            else
                valid = false(1, length(t_total));
                for j = 1:size(db(k).OTtimes{probe}{OT_ID},1)
                    valid = valid | (t_total > db(k).OTtimes{probe}{OT_ID}(j,1) & ...
                        t_total < db(k).OTtimes{probe}{OT_ID}(j,2));
                end
            end
            indD = ind_dark & valid;
            spikeRate_z = zscore(spikeRate(indD));
            spikeAmp_z = zscore(spikeAmp(indD));
            [xc,lags] = xcorr(spikeRate_z, spikeAmp_z, 500/binW, 'unbiased');
            xc_control = zeros(numShifts, length(lags));
            shifts = randi(sum(indD), 1, numShifts);
            for j = 1:numShifts
                xc_control(j,:) = xcorr(circshift(spikeRate_z, shifts(j)), ...
                    spikeAmp_z, 500/binW, 'unbiased');
            end
            percentiles = prctile(xc_control,[2.5 50 97.5],1);
            
            % Visual noise: Find pixels that trigger response, then plot raster
            % for these pixels
            respPerFrame = histcounts(st, [noiseFrameTimes(:); ...
                noiseFrameTimes(end)+ framesPerImage / myScreenInfo.FrameRate]) / ...
                (framesPerImage / myScreenInfo.FrameRate);
            respPerFrame = respPerFrame - mean(respPerFrame);
            stims = reshape(stimFrames, [], size(stimFrames,3))';
            stims = repmat(stims, size(noiseFrameTimes,2), 1);
            kernels = abs(stims) \ respPerFrame';
            pix = reshape(kernels, size(stimFrames,1), size(stimFrames,2));
            [RFmaxi, ind] = max(abs(kernels));
            pixTimes = noiseFrameTimes(stims(:,ind) ~= 0);
            [psth, bins, rasterX, rasterY] = psthAndBA(st, pixTimes, [-.1 .3], .001);
        
            figure('Position', [1 41 1920 1083])
            
            % plot corr of spike rate and amplitude
            subplot(plotRows, plotCols, plotCols .* [1 2 3])
            hold on
            fill([lags, flip(lags)].*binW, [percentiles(1,:), flip(percentiles(3,:))], ...
                'k', 'EdgeColor', 'none', 'FaceColor', [.5 .5 .5], 'FaceAlpha', .3)
            plot(lags.*binW, percentiles(2,:), ':', 'Color', [.5 .5 .5], 'LineWidth', 1)
            plot(lags.*binW, xc, 'k', 'LineWidth', 1)
            plot(lags([1 end]).*binW, [0 0], 'k:')
            ylimits = get(gca, 'YLim');
            plot([0 0], ylimits, 'k:')
            ylabel('Corr(firing rate, spike amplitude)')
            xlabel('Lag of spike rate (t + lag, in sec)')
            set(gca, 'box', 'off', 'Position', [.82 .68 .1 .29])
        
            % plot spike rate
            h = [0 0];
            legLabels = {'vis. noise', 'gratings'};
            subplot(plotRows, plotCols, 1:plotCols-1)
            hold on
            plot(t_total, spikeRate, 'k')
            maxi = max(spikeRate);
            mini = min(spikeRate);
            h(1) = plot([1 1].*noiseFrameTimes(1), [mini maxi], 'Color', cols(1,:), ...
                'LineWidth', 2);
            plot([1 1].*noiseFrameTimes(end), [mini maxi], ':', 'Color', cols(1,:), ...
                'LineWidth', 2)
            h(2) = plot([1 1].*stimOn(1), [mini maxi], 'Color', cols(2,:), ...
                'LineWidth', 2);
            plot([1 1].*stimOff(end), [mini maxi], ':', 'Color', cols(2,:), ...
                'LineWidth', 2)
            if ~isempty(flickerTimes)
                h(end+1) = plot([1 1].*flickerAll(1), [mini maxi], 'Color', cols(3,:), ...
                    'LineWidth', 2);
                plot([1 1].*flickerAll(end), [mini maxi], ':', 'Color', cols(3,:), ...
                    'LineWidth', 2)
                legLabels{end+1} = 'flicker';
            end
            if ~isempty(db(k).darkTime)
                h(end+1) = plot([1 1].*db(k).darkTime(1), [mini maxi], 'Color', cols(4,:), ...
                    'LineWidth', 2);
                plot([1 1].*db(k).darkTime(2), [mini maxi], ':', 'Color', cols(4,:), ...
                    'LineWidth', 2)
                legLabels{end+1} = 'darkness';
            end
            if ~isempty(db(k).darkTime_IR_off)
                plot([1 1].*db(k).darkTime_IR_off(1), [mini maxi], 'Color', cols(4,:), ...
                    'LineWidth', 2);
                plot([1 1].*db(k).darkTime_IR_off(2), [mini maxi], ':', 'Color', cols(4,:), ...
                    'LineWidth', 2)
            end
            if ~isempty(db(k).expOriTropic)
                plot([1 1].*stimOnTropic(1), [mini maxi], 'Color', cols(2,:), ...
                    'LineWidth', 2);
                plot([1 1].*stimOffTropic(end), [mini maxi], ':', 'Color', cols(2,:), ...
                    'LineWidth', 2)
            end
            leg = legend(h, legLabels, 'Location', 'NorthWest');
            leg.Position = [.03 .9 .07 .07];
            axis tight
            ylabel('Firing rate')
            set(gca, 'XTickLabel', [], 'Position', [.13 .9 .65 .07])
            
            % plot spike amplitude
            subplot(plotRows, plotCols, plotCols+(1:plotCols-1))
            h = plot(st, sa, 'k.', 'MarkerSize', .5);
            axis([t_total([1 end]) 0 max(sa)])
            ylabel('amp.')
            set(gca, 'box', 'off', 'XTickLabel', [], 'Position', [.13 .82 .65 .0673])
            
            % density plot of spike amplitudes
            subplot(plotRows, plotCols, reshape([2;3].*plotCols+(1:plotCols-1),1,[]))
            edges = {t_total(1):10:t_total(end), linspace(min(sa),prctile(sa,99.9),40)};
            n = hist3([st, sa], 'Edges', edges);
            imagesc(t_total([1 end]), [min(sa) max(sa)], n')
            ylabel('amp.')
            colormap hot
            set(gca,'YDir','normal', 'XTickLabel', [], 'Position', [.13 .68 .65 .12])
            
            % plot spike depth
            subplot(plotRows, plotCols, 4*plotCols+(1:plotCols-1))
            plot(st, sd, 'k.', 'MarkerSize', 1)
            axis tight
            ylabel('depth')
            xlabel('Time (s)')
            set(gca, 'box', 'off', 'Position', [.13 .6 .65 .0673])
        
            % plot waveform of this cell
%             subplot(plotRows, plotCols, 3*plotCols+1:plotCols:plotsTotal)
            subplot(plotRows, plotCols, reshape((4:7)'.*plotCols+1,1,[]))
            hold on
            thisTemp = squeeze(sp(probe).tempsUnW(templates(unitID+1)+1,:,:));
            [~,peakChan] = max(max(abs(thisTemp),[],1),[],2);
            plotChans = max(1,peakChan-10) : min(length(sp(probe).xcoords),peakChan+10);
            arrayfun(@(x)plot(sp(probe).xcoords(x)+0.3*(1:size(thisTemp,1))', ...
                sp(probe).ycoords(x)+5*thisTemp(:,x), 'k'), plotChans)
            ycoords = unique(sp(probe).ycoords(plotChans));
            mini = min(ycoords);
            maxi = max(ycoords);
            chanDist = median(diff(ycoords));
            ylim([mini - 0.5*chanDist, maxi + 0.5 * chanDist])
            ylabel('Depth (um)')
            set(gca, 'XTick', [])
            
            % plot autocorrelogram
            subplot(plotRows, plotCols, reshape(8.*plotCols+1,1,[]))
            [xLin, nLin] = myACG(st,[],gca);
            ylabel('Autocorr.')
            xlabel('Lag (ms)')
            set(gca,'Position',[.13 .08 .1001 .1])
            
            if ~isempty(flickerTimes)
                % plot responses to flickering monitor
                psth_f = cell(size(flickerTimes,1), size(flickerTimes,2), 2);
                rasterX_f = cell(size(flickerTimes,1), size(flickerTimes,2), 2);
                rasterY_f = cell(size(flickerTimes,1), size(flickerTimes,2), 2);
                if ~isempty(flickerTimes)
                    for stim = 1:size(flickerTimes,1)
                        for rep = 1:size(flickerTimes,2)
                            flicks = flickerTimes{stim,rep};
                            for lum = 1:2 % white and black)
                                [psth_f{stim,rep,lum}, bins_f, rasterX_f{stim,rep,lum}, ...
                                    rasterY_f{stim,rep,lum}] = ...
                                    psthAndBA(st, flicks(:,lum), [-.05 .3], .001);
                            end
                        end
                    end
                end
                
                subplot(plotRows, plotCols, reshape(8*plotCols+(2:3),1,[]))
                rates = NaN(length(bins_f),size(flickerTimes,1));
                for stim = 1:size(flickerTimes,1)
                    rates(:,stim) = conv(psth_f{stim,2}, winHamm, 'same');
                end
                fill(bins_f([1 end end 1]), [0 0 [1 1].* ...
                    max(rates(:))*1.05], 'k', ...
                    'EdgeColor', 'none', 'FaceAlpha', .1)
                hold on
                h = zeros(1,size(flickerTimes,1));
                for stim = 1:size(flickerTimes,1)
                    h(stim) = plot(bins_f, rates(:,stim), 'Color', colsFlicker(stim,:));
                end
                plot([0 0],[0 max(rates(:))*1.05], 'k:')
                axis tight
                xlabel('Time from onset')
                ylabel('Firing rate')
                set(gca, 'Position', [0.2645 0.08 0.2368 0.0673])
                
                subplot(plotRows, plotCols, reshape(7*plotCols+(2:3),1,[]))
                hold on
                rates = NaN(length(bins_f),size(flickerTimes,1));
                for stim = 1:size(flickerTimes,1)
                    rates(:,stim) = conv(mean(cat(1,psth_f{stim,:,1})), winHamm, 'same');
                    plot(bins_f, rates(:,stim), 'Color', colsFlicker(stim,:));
                end
                plot([0 0],[0 max(rates(:))*1.05], 'k:')
                c = legend(num2str(flickerDurs'));
                c.Position = [0.4462 0.09 0.0479 0.13];
                axis tight
                ylabel('Firing rate')
                set(gca, 'XTickLabel', [], 'Position', [0.2645 0.17 0.2368 0.0673])
                
                subplot(plotRows, plotCols, reshape((4:6)'.*plotCols+(2:3),1,[]))
                hold on
                m = 0;
                yLabPos = zeros(size(flickerTimes,1),1);
                for stim = 1:size(flickerTimes,1)
                    for lum = 1:2
                        if lum == 2
                            if ~isempty(cat(2,rasterY_f{stim,:,lum}))
                                y = sum(cellfun(@max,rasterY_f(stim,:,lum)));
                                fill(bins_f([1 end end 1]), ...
                                    m + [0 0 [1 1].* y], 'k', ...
                                    'EdgeColor', 'none', 'FaceAlpha', .1)
                            end
                            yLabPos(stim) = m;
                        end
                        for rep = 1:size(flickerTimes,2)
                            plot(rasterX_f{stim,rep,lum}, rasterY_f{stim,rep,lum}+m, 'k')
                            m = m + max([rasterY_f{stim,lum}(:);0]);
                        end
                    end
                end
                plot([0 0], [0 m], 'k')
                axis([bins_f([1 end]) 0 max(m,1)])
                yLab = strsplit(sprintf('%.2f ', flickerDurs));
                yLab = yLab(1:end-1);
                noResp = all(cellfun(@isempty,rasterY_f),2);
                yLabPos(noResp) = [];
                yLab(noResp) = [];
                set(gca,'YDir','reverse', 'YTick', yLabPos, 'YTickLabel', ...
                    yLab, 'XTickLabel', [], 'Position', [0.2645 0.2569 0.2368 0.2942])
                ylabel('Stim duration')
                title('Flickering monitor')
            end
            
%             p = cell(1,size(flickerTimes,2));
% %             ind = 7:5:length(bins_f)-2;
% %             ind = 4:3:length(bins_f)-1;
%             ind = 1:length(bins_f);
%             for rep = 1:size(flickerTimes,2)
% %                 p{rep} = mean(reshape(psth_f{1,rep,1}(5:ind(end)+2),5,[]),1);
% %                 p{rep} = mean(reshape(psth_f{1,rep,1}(3:ind(end)+1),3,[]),1);
%                 p{rep} = psth_f{1,rep,1};
%             end
%             figure
%             hold on
%             m = max(cellfun(@max,p));
%             b = bins_f(ind);
%             ft = median(diff(flickerTimes{1,1}(:,1)))/2;
%             ft = (-4:21).*ft;
%             stairs(ft, reshape(repmat([1 0]',1,length(ft)/2),[],1) + m, 'k')
%             for rep = 1:size(flickerTimes,2)
%                 stairs(b, p{rep}-(rep-1)*m, 'k')
%             end
%             xlim([-0.05 0.3])
%             xlabel('Time (s)')
%             ylabel('Firing rate (sp/s)')
            
            % plot responses to visual noise
            if RFmaxi > 0
                % spatial RF
                if isempty(flickerTimes)
                    subplot(plotRows, plotCols, reshape((4:5)'*plotCols+(2:3),1,[]))
                else
                    subplot(plotRows, plotCols, reshape((4:5)'.*plotCols+(4:6),1,[]))
                end
                imagesc(noisePosition([1 2]), noisePosition([3 4]), ...
                    pix, [-RFmaxi RFmaxi])
                c = colorbar;
                c.Label.String = 'spikes/s';
                title('Visual noise (absolute RF)')
                xlabel('Azimuth (vis. deg.)')
                ylabel('Altitude (vis. deg.)')
                colormap(gca,cm)
                axis image
                
                if isempty(flickerTimes)
                    % raster plot for most driving pixel
                    subplot(plotRows, plotCols, reshape((6:7)'.*plotCols+(2:3),1,[]))
                    plot(rasterX, rasterY, 'k')
                    xlim(bins([1 end]))
                    if ~isempty(rasterY)
                        ylim([1 max(rasterY)])
                    end
                    set(gca, 'YDir', 'reverse', 'box', 'off', 'XTickLabel', [])
                    ylabel('Spike raster for max pixel')
                end
            end
            % PSTH
            if isempty(flickerTimes)
                subplot(plotRows, plotCols, reshape([8].*plotCols+(2:3),1,[]))
                plot(bins, conv(psth, winHamm, 'same'), 'k')
                yLimits = get(gca, 'YLim');
                hold on
                plot([0 0], yLimits, 'k:')
                plot([1 1] .* framesPerImage / myScreenInfo.FrameRate, yLimits, 'k:')
                xlim(bins([1 end]))
                set(gca, 'box', 'off', 'Position', [.2645 .08 .2368 .1])
                ylabel('PSTH')
                xlabel('Time from frame onset (s)')
            end
        
            % plot responses to gratings
            gratingPsths = NaN(parsGratings.npfilestimuli, length(binsGrating));
            gratingSpikes = struct('times', cell(parsGratings.npfilestimuli, ...
                parsGratings.nrepeats));
            for stim = 1:parsGratings.npfilestimuli
                ind = find(stimSeq == stim);
                gratingPsths(stim,:) = psthAndBA(st, stimOn(ind), window, ...
                    binSizeGrating);
                gratingPsths(stim,:) = conv(gratingPsths(stim,:), winHamm, 'same');
                for tr = 1:parsGratings.nrepeats
                    s = st - stimOn(ind(tr));
                    gratingSpikes(stim,tr).times = s(s>=0 & s <= stimDur);
                end
            end
            m_spont_gratings = mean(gratingPsths(blank,stimBins));
            % mean PSTH
            if isempty(flickerTimes)
                subplot(plotRows, plotCols, reshape([4;5].*plotCols+(4:6),1,[]))
                hold on
                plot(binsGrating, gratingPsths, 'Color', [1 1 1].*.5)
                plot(binsGrating, mean(gratingPsths), 'k', 'LineWidth', 2)
                plot([0 0], [min(gratingPsths(:)), max(gratingPsths(:))], ...
                    'k:', 'LineWidth', 2)
                plot([1 1].*stimDur, [min(gratingPsths(:)), max(gratingPsths(:))], ...
                    'k:', 'LineWidth', 2)
                axis tight
                set(gca, 'box', 'off', 'XTickLabel', [])
                title('Gratings')
                ylabel('PSTH')
            end
            % timecourse of response to each stim
            if isempty(flickerTimes)
                subplot(plotRows, plotCols, reshape([6;7].*plotCols+(4:6),1,[]))
            else
                subplot(plotRows, plotCols, reshape([6;7].*plotCols+(4:6),1,[]))
                title('Gratings')
            end
            imagesc(binsGrating([1 end]),[1 parsGratings.npfilestimuli], gratingPsths)
            hold on
            plot([0 0], [0.5 parsGratings.npfilestimuli+.5], 'w:', 'LineWidth', 2)
            plot([1 1].*stimDur, [0.5 parsGratings.npfilestimuli+.5], 'w:', 'LineWidth', 2)
            axis tight
            set(gca, 'box', 'off', 'YTick', [1:2:length(directions), ...
                parsGratings.npfilestimuli], ...
                'YTickLabel', [cellstr(num2str(directions(1:2:end)'));'\0'])
            ylabel('Stim. responses')
            colormap(gca, 'hot')
            c = colorbar;
            c.Label.String = 'spikes/s';
            c.Position = [0.93 0.2031 0.0139 0.1607];
            % tuning curve
            F1 = NaN(size(gratingPsths,1),1);
            for stim = 1:size(gratingPsths,1)
                [J,Msp,Nsp]=mtfftpt(gratingSpikes(stim,:),tap,nfft,t,f,findx); % mt fft for point process times
                J = mean(mean(abs(J/N),2),3);
                F1(stim) = 2 * interp1(f,J,tf_gratings);
            end
            F1 = F1(~blank);
            F0 = mean(gratingPsths(~blank,stimBins),2);
            if max(F0) > max(F1)
                gratingResp = F0;
                yLab = 'F0 (spikes/s)';
            else
                gratingResp = F1;
                yLab = 'F1 (spikes/s)';
            end
            if isempty(flickerTimes)
                subplot(plotRows, plotCols, reshape(8.*plotCols+(4:6),1,[]))
            else
                subplot(plotRows, plotCols, reshape(8.*plotCols+(4:6),1,[]))
            end
            hold on
            plot(directions, gratingResp, 'o-k', 'MarkerFaceColor', 'k')
            plot(directions([1 end]), [1 1] .* ...
                mean(gratingPsths(blank,stimBins)), 'k:')
            xlim([directions(1)-10 directions(end)+10])
            xlabel('Direction')
            ylabel(yLab)
            set(gca,'box','off', 'Position', [.5336 .08 .3714 .08])
        
            str = 'not sorted';
            if ismember(unitID, units_good)
                str = 'single unit';
            elseif ismember(unitID, units_mua)
                str = 'multi-unit';
            end
            annotation('textbox', [0 .958 1 .04], 'String', sprintf(...
                'Unit %d (%s)', unitID, str), ...
                'FontSize', 15, 'FontWeight', 'bold', ...
                'LineStyle', 'none', 'HorizontalAlignment', 'center')
%             depth = mean(sp(probe).spikeDepths(sp(probe).clu == unitID));
            annotation('textbox', [.03 .8 .09 .05], 'String', sprintf(...
                'Depth: %d um', round(depth)), ...
                'FontSize', 11, 'LineStyle', 'none', 'HorizontalAlignment', 'left')
        
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(fullfile(folder, sprintf('%04d_neuron%04d.jpg', ...
                round(depth), unitID)), '-djpeg','-r0')
            
            if sum(units_OT == unitID) == 0
                close(fig)
            end
        end
    end
    
%     goodSpikes = ismember(sp.clu, goodUnits);
%     psthViewer(sp.st(goodSpikes), sp.clu(goodSpikes), stimOn, window, ...
%         stimSeq)
end

%% (3) Correlation with running

% Calculate correlations

runningSigma = 0.25; % sec, for smoothing of running speed
% smoothingSigma = [.1 .2 .5 1 2 5 10 20 50]; % 0.5; % sec, for smoothing spike rate
smoothingSigma = 1; % 0.5; % sec, for smoothing spike rate
highPassWindow = 180; % in sec; to high-pass filter firing rates
prctileFilter = 8;

maxLag = 100; % in sec
numShifts = 500; % for significance test of correlations

runningThreshold = 1; %(cm/s)
% binSizeRun = 0.005;
binSizeRun = 1/7.5;
binSizeRun_kernel = 0.05;
reduce = round(binSizeRun_kernel / binSizeRun);

runningDeriv = false;

mLag = round(maxLag / binSizeRun);
sigma = round(smoothingSigma / binSizeRun);
wins = cell(1, length(sigma));
for j = 1:length(sigma)
    wins{j} = normpdf(-5*sigma(j) : 5*sigma(j), 0, sigma(j));
end
            
rows = 5;
cols = 3;
timeScales = 1; %[3 5 7];
% sigma_highPass = round(100 / binSizeRun);
% win_highPass = normpdf(-5*sigma_highPass : 5*sigma_highPass, 0, sigma_highPass);
a = gcp('nocreate');
if isempty(a)
    parpool('local', 4);
end

corrs = struct([]);
set = 1;

for k = 1:length(db)
    fprintf('\nProcessing: %s %s\n', db(k).subject, db(k).date);
    alignDir = fullfile(subjectsFolder, db(k).subject, db(k).date, 'alignments');
    %% Load spike data
    sp = loadAllKsDir(db(k).subject, db(k).date);
    
    [expNums, ~, ~, ~, ~, tl, hasTimeline] = ...
        dat.whichExpNums(db(k).subject, db(k).date);
    TLexp = expNums(hasTimeline);
    TLexp = TLexp(end);
    tl = tl{end}; % ?? Correct?
    bTLtoMaster = readNPY(fullfile(alignDir, ...
        sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, db(k).probeNames{1})));
    tlTime = applyCorrection(tl.rawDAQTimestamps, bTLtoMaster);
    
    % Load and prepare data for running correlation
    rotary = double(tl.rawDAQData(:, strcmp({tl.hw.inputs.name}, ...
        'rotaryEncoder')));
    runSpeed = nonVis.getRunningSpeed_wheel(rotary, tlTime, runningSigma);
    runningTime = runSpeed.t;
    runSpeed = runSpeed.total;
    cmPerUnit = 2*pi * 8.75 / (4 * 1024);
    runSpeed = runSpeed * cmPerUnit;
    
    binEdges = db(k).darkTime(1) : binSizeRun : db(k).darkTime(2);
    timeBins = binEdges(1:end-1) + binSizeRun/2;
    
    run = interp1(runningTime, runSpeed, timeBins, 'pchip');
    run = [ones(1,length(wins{end})) .* ...
        mean(run(1:min(length(run),length(wins{end})))), ...
        run, ones(1,length(wins{end})) .* ...
        mean(run(end-min(length(run),length(wins{end}))+1:end))];
    
    running = cell(length(sigma), 1);
    for j = 1:length(sigma)
        running{j} = conv(run, wins{j}, 'same');
        running{j} = running{j}(length(wins{end})+1 : end-length(wins{end}));
        if runningDeriv
            running{j} = interp1(binEdges(2:end-1), diff(running{j}), timeBins, 'pchip');
            running{j} = conv(running{j}, wins{j}, 'same');
        else
            [running{j}, smoothed] = preproc.removeSlowDrift(running{j}', 1/binSizeRun, ...
                highPassWindow, prctileFilter);
            running{j} = running{j}' + mean(smoothed);
        end
    end
    
    for probe = 1:length(sp)
        corrs(set).subject = db(k).subject;
        corrs(set).date = db(k).date;
        corrs(set).probeName = sp(probe).name;
        corrs(set).time = timeBins;
        corrs(set).running = running{1};
        
        plotF = fullfile(plotFolder, 'plots_runningCorr', ...
            [db(k).subject '_' db(k).date '_' sp(probe).name]);
        if ~isfolder(plotF)
            mkdir(plotF)
        end
        
        units_OT = db(k).OTunits{probe};
        unit = 1;
        for iCell = 1:length(units_OT)
            if db(k).OTgood{probe}(iCell) == 0
                continue
            end
            fprintf('  Unit %d of %d\n', iCell, length(units_OT))
            corrs(set).units(unit).ID = units_OT(iCell);
            % get spiketimes for this cell
            st = sp(probe).st(sp(probe).clu == units_OT(iCell));
            
            % determine valid times of darkness
            darkTimes = db(k).darkTime;
            validTime = db(k).OTtimes{probe}{iCell};
            if ~isempty(validTime)
                darkTimes = [];
            end
            for j = 1:size(validTime,1)
                if (validTime(j,1)>db(k).darkTime(1) && ...
                        validTime(j,1)<db(k).darkTime(2)) || ...
                        (validTime(j,2)>db(k).darkTime(1) && ...
                        validTime(j,2)<db(k).darkTime(2))
                    darkTimes(end+1,:) = [max(validTime(j,1),db(k).darkTime(1)), ...
                        min(validTime(j,2),db(k).darkTime(2))];
                end
            end
            
            % get periods of running and spiking data and filter at
            % different time scales
            numD = size(darkTimes,1);
            validInds = cell(1, numD);
            binsRun = cell(1, numD);
            time_parts = cell(1, numD);
            run = cell(length(sigma), numD);
%             running_kernel = cell(1, numD);
            spikerate = cell(length(sigma), numD);
            spiketimes = cell(1, numD);
%             toeplitz = cell(1, numD);
%             spikerate_kernel = cell(1, numD);
%             spikerate_kernel_good = cell(1, numD);
            corrs(set).units(unit).spikerate = NaN(1, length(timeBins));
            for d = 1:numD
%                 binsRun{d} = darkTimes(d,1) : binSizeRun : darkTimes(d,2);
                validInds{d} = find(binEdges>=darkTimes(d,1),1) : ...
                    find(binEdges<=darkTimes(d,2),1,'last');
                binsRun{d} = binEdges(validInds{d});
                time_parts{d} = timeBins(validInds{d}(1:end-1));
                spiketimes{d} = st(st >= binsRun{d}(1) & st <= binsRun{d}(end))';
                [sr,~,bin] = histcounts(spiketimes{d}, binsRun{d});
                sr = [ones(1,length(wins{end})) .* ...
                    mean(sr(1:min(length(sr),length(wins{end})))), ...
                    sr, ones(1,length(wins{end})) .* ...
                    mean(sr(end-min(length(sr),length(wins{end}))+1:end))] ./ binSizeRun;
                
                for j = 1:length(sigma)
                    run{j,d} = running{j}(validInds{d}(1:end-1));
                    spikerate{j,d} = conv(sr, wins{j}, 'same');
                    spikerate{j,d} = spikerate{j,d}(length(wins{end})+1 : ...
                        end-length(wins{end}));
                    [spikerate{j,d}, smoothed] = preproc.removeSlowDrift( ...
                        spikerate{j,d}', 1/binSizeRun, ...
                        highPassWindow, prctileFilter);
                    spikerate{j,d} = spikerate{j,d}' + mean(smoothed);
                end
                %         time_parts{d}(end+1) = NaN;
                %         if d>1
                %             time_parts{d} = time_parts{d} - (time_parts{d}(1)-time_parts{d-1}(end-1)) + 100;
                %         end
                
                %             running_kernel{d} = decimate(running{d}, reduce, 'fir');
%                 sr = [ones(1,length(win_highPass)).*mean(spikerate{d}(1:length(win_highPass))), ...
%                     spikerate{d}, ...
%                     ones(1,length(win_highPass)).*mean(spikerate{d}(end-length(win_highPass)+1:end))];
%                 spikerate_slow = conv(sr, win_highPass, 'same');
%                 spikerate_slow = spikerate_slow(length(win_highPass)+1 : ...
%                     end-length(win_highPass));
%                 spikerate{d} = spikerate{d} - spikerate_slow + min(spikerate_slow);
%                 spikerate_kernel{d} = decimate(spikerate{d}, reduce, 'fir');
                %             [rho(d),pVal(d)] = corr(run{d}', spikerateSmooth{d}');
                %             run{d}(end+1) = NaN;
                %             spikerateSmooth{d}(end+1) = NaN;
                
                corrs(set).units(unit).spikerate(validInds{d}(1:end-1)) = ...
                    spikerate{1,d};
            end
            
            % get means and stds across data parts, use for z-scoring
            time = cat(2, time_parts{:});
            running_total = cell(length(sigma),1);
            run_mean = zeros(length(sigma),1);
            run_std = zeros(length(sigma),1);
            spikerate_total = cell(length(sigma),1);
            sr_mean = zeros(length(sigma),1);
            sr_std = zeros(length(sigma),1);
            for j = 1:length(sigma)
                running_total{j} = cat(2, run{j,:});
                run_mean(j) = mean(running_total{j});
                run_std(j) = std(running_total{j});
                spikerate_total{j} = cat(2, spikerate{j,:});
                sr_mean(j) = mean(spikerate_total{j});
                sr_std(j) = std(spikerate_total{j});
            end
            corrs(set).units(unit).runningTime = sum(running_total{j} > runningThreshold) / ...
                length(running_total{j});
            %         time_kernel = cellfun(@downsample, time_parts, ...
            %             num2cell(repmat(reduce,1,numD)), 'UniformOutput', false);
            
            % calculate cross-correlation at different time scales
            crossCorr = cell(length(sigma), numD);
            crossCorr_cv = cell(length(sigma), 1);
            shifts = randi(length(time), 1, numShifts);
            for j = 1:length(sigma)
                for d = 1:numD
                    [crossCorr{j,d},lags] = ...
                        xcorr((spikerate{j,d} - sr_mean(j)) ./ sr_std(j), ...
                        (run{j,d} - run_mean(j)) ./ run_std(j), ...
                        mLag, 'unbiased');
                end
                cc = zeros(numShifts, length(lags));
                sr_total = spikerate_total{j};
                parfor sh = 1:numShifts
                    sr = circshift(sr_total, shifts(sh));
                    cc(sh,:) = xcorr((sr - sr_mean(j)) ./ sr_std(j), ...
                        (running_total{j} - run_mean(j)) ./ run_std(j), ...
                        mLag, 'unbiased');
                end
                crossCorr_cv{j} = cc;
                
%                 [toeplitz{d}, numSamples, windowTimes]  = krnl.getToeplitz( ...
%                     time_kernel{d}(t1:t2), [], [], {running_kernel{d}(t1:t2)}, {[-3 3]});
%                 spikerate_kernel{d} = spikerate_kernel{d}(t1:t2);
            end
%             toeplitz = cat(1, toeplitz{:});
%             spikerate_kernel = cat(2, spikerate_kernel{:})';
%             kernel = toeplitz \ spikerate_kernel;
%             prediction = toeplitz * kernel;
            
            spiketimes = cat(2, spiketimes{:});
            spikerate_total = cell(length(sigma),1);
            running_total = cell(length(sigma),1);
            time = [];
            for d = 1:numD
                for j = 1:length(sigma)
                    spikerate_total{j} = [spikerate_total{j}, spikerate{j,d}, NaN];
                    running_total{j} = [running_total{j}, run{j,d}, NaN];
                end
                time = [time, time_parts{d}, NaN];
            end
            crossCorr_total = cell(length(sigma),1);
            for j = 1:length(sigma)
                crossCorr_total{j} = sum(cat(1, crossCorr{j,:}) .* ...
                    cellfun(@length,time_parts)' ./ length(time),1);
            end
            lags = lags .* binSizeRun;
            
            corrs(set).units(unit).crosscorr = crossCorr_total{1};
            corrs(set).units(unit).nullCrossCorr = crossCorr_cv{1};
            corrs(set).lags = lags;
            
            figure('Position', [3 41 1920 1080])
            
            colors = lines(length(timeScales));
            lw = 0.5 .* (1:length(timeScales));
            ax = zeros(1,2);
            subplot(rows,cols,1:cols)
            hold on
            for s = 1:length(timeScales)
                plot(time, running_total{timeScales(s)}, 'Color', colors(s,:), ...
                    'LineWidth', lw(s))
            end
            a = gca;
            a.Box = 'off';
            ax(1) = a;
            ylabel('Running speed (cm/s)')
            
            subplot(rows,cols,cols+(1:cols))
            hold on
            for s = 1:length(timeScales)
                plot(time, spikerate_total{timeScales(s)}, 'Color', colors(s,:), ...
                    'LineWidth', lw(s))
            end
            a = gca;
            a.Box = 'off';
            ax(2) = a;
            ylabel('Spikes/s')
            axis tight
            xlabel('Time (s)')
            linkaxes(ax, 'x')
            xlim(time([1 end-1]))
            l = legend(cellstr([num2str(smoothingSigma(timeScales)'), ...
                repmat(' s',length(timeScales),1)]));
            l.Position = [.92 .7 .04 .04];
            
            maxis = zeros(1, length(sigma));
            for j = 1:length(sigma)
%                 subplot(rows,cols,2*cols+j)
                subplot(rows,cols,3*cols+(1:(rows-3)*cols))
                hold on
                pc = prctile(crossCorr_cv{j},[2.5 50 97.5],1);
                maxis(j) = max([reshape(abs(pc([1 3],:)),1,[]), ...
                    abs(crossCorr_total{j})]);
                fill([lags,flip(lags)], [pc(1,:),flip(pc(3,:))], 'k', ...
                    'EdgeColor', 'none', 'FaceColor', ones(1,3).*.5, ...
                    'FaceAlpha', .3)
                plot(lags, pc(2,:), ':', 'Color', ones(1,3).*.5, 'LineWidth', 1)
                plot(lags, crossCorr_total{j}, 'k', 'LineWidth', 1)
                plot(lags([1 end]), [0 0], 'k:')
                if j == 1 %8
                    xlabel('Time lag (s) (of spike rate)')
                else
                    set(gca, 'XTickLabel',[])
                end
                if j == 1
                    ylabel('Cross-corr.')
                end
                if mod(j,3) ~= 1
                    set(gca, 'YTickLabel', [])
                end
                title(sprintf('%.1f sec', smoothingSigma(j)))
            end
            for j = 1:length(sigma)
%                 m = max(maxis(ceil(j/3)*3 + (-2:0)));
                m = maxis;
%                 subplot(rows,cols,2*cols+j)
                subplot(rows,cols,3*cols+(1:(rows-3)*cols))
                plot([0 0],[-m m], 'k:')
                ylim([-m m])
            end
            
            annotation('textbox', [.01 .8 .1 .1], 'String', sprintf(...
                '%s %s %s\nUnit %d', db(k).subject, db(k).date, ...
                sp(probe).name, units_OT(iCell)), ...
                'FontSize', 12, 'FontWeight', 'bold', ...
                'LineStyle', 'none', 'HorizontalAlignment', 'left')
            
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(fullfile(plotF, sprintf('unit%04d.jpg', ...
                units_OT(iCell))), '-djpeg','-r0')
            savefig(fig, fullfile(plotF, sprintf('unit%04d.fig', units_OT(iCell))))
            close(fig)
            
            unit = unit + 1;
        end
        set = set + 1;
    end
end
save(fullfile(plotFolder, 'runningCorrelation_dark', 'correlations_runningDerivative.mat'), ...
    'corrs')

% Plot correlations
data = load(fullfile(plotFolder, 'runningCorrelation_dark', 'correlations_runningFiltered.mat'));
corrs = data.corrs;

rhos = [];
pVals = [];
signif = [];
runTime = [];
for k = 1:length(corrs)
    [~,zerolag] = min(abs(corrs(k).lags));
    for n = 1:length(corrs(k).units)
        rhos(end+1) = corrs(k).units(n).crosscorr(zerolag);
        nRhos = sort(corrs(k).units(n).nullCrossCorr(:,zerolag));
        score = (find(nRhos >= rhos(end), 1)-1) / ...
            size(corrs(k).units(n).nullCrossCorr,1);
        if isempty(score), score = 0; end
        if score > 0.5, score = 1-score; end
        pVals(end+1) = 2 * score;
        prct = prctile(corrs(k).units(n).nullCrossCorr(:,zerolag), ...
            [2.5 97.5]);
        signif(end+1) = rhos(end)<prct(1) || rhos(end)>prct(2);
        runTime(end+1) = corrs(k).units(n).runningTime;
    end
end
rhos(runTime < 0.05) = [];
pVals(runTime < 0.05) = [];
signif(runTime < 0.05) = [];
signif = signif == 1;

% Fisher's method to combine p-values
pVals(pVals==0) = 1/size(corrs(1).units(1).nullCrossCorr,1);
chi_vals = -2.*log(pVals);
group_pval = 1 - chi2cdf(sum(chi_vals),2*length(pVals));

bins = -0.4:.1:.3;
n1 = hist(rhos(signif),bins)';
n2 = hist(rhos(~signif),bins)';
figure
b = bar(bins,[n1,n2], 'stacked');
b(1).FaceColor = 'k';
b(2).FaceColor = 'w';
xlabel('Correlation coeff.')
ylabel('#Synapses')
title(sprintf('Darkness - running (n = %d, %d%% signif., 1 s filter)', ...
    length(rhos), round(sum(signif) / length(rhos)*100)))
xlim([bins(1)-.05 bins(end)+.05])
ax = gca;
ax.Box = 'off';
legend('p < 0.05', 'p >= 0.05')

%% (4) Running-triggered average (not finished because didn't look promising)

% Parameters
runningThreshold = 1; % cm/s
noRunBefore = 1; % in sec
minRunPeriod = 1; % in sec

runningSigma = 0.25; % sec, for smoothing of running speed
binSizeRun = 0.01;
% binSizeRun = 1/7.5;

smoothingSigma = 1; % 0.5; % sec, for smoothing spike rate
highPassWindow = 180; % in sec; to high-pass filter firing rates
prctileFilter = 8;


rta = struct([]);
count = 1;

for k = 1:length(db)
    fprintf('\nProcessing: %s %s\n', db(k).subject, db(k).date);
    alignDir = fullfile(subjectsFolder, db(k).subject, db(k).date, 'alignments');
    %% Load spike data
    sp = loadAllKsDir(db(k).subject, db(k).date);
    
    [expNums, ~, ~, ~, ~, tl, hasTimeline] = ...
        dat.whichExpNums(db(k).subject, db(k).date);
    TLexp = expNums(hasTimeline);
    TLexp = TLexp(end);
    tl = tl{end}; % ?? Correct?
    bTLtoMaster = readNPY(fullfile(alignDir, ...
        sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, db(k).probeNames{1})));
    tlTime = applyCorrection(tl.rawDAQTimestamps, bTLtoMaster);
    
    % Load and prepare data for running correlation
    rotary = double(tl.rawDAQData(:, strcmp({tl.hw.inputs.name}, ...
        'rotaryEncoder')));
    runSpeed = nonVis.getRunningSpeed_wheel(rotary, tlTime, runningSigma);
    runningTime = runSpeed.t;
    runSpeed = runSpeed.total;
    cmPerUnit = 2*pi * 8.75 / (4 * 1024);
    runSpeed = runSpeed * cmPerUnit;
    
    binEdges = db(k).darkTime(1) : binSizeRun : db(k).darkTime(2);
    timeBins = binEdges(1:end-1) + binSizeRun/2;
    
    run = interp1(runningTime, runSpeed, timeBins, 'pchip');
    running_binary = run > runningThreshold;
    sta = find(diff(running_binary)==1); % start
    sto = find(diff(running_binary)==-1); % stop
    % if animal is running at start of recording, discard end of 1st
    % running period
    if sta(1)>sto(1)
        sta = [0 sta];
    end
    if length(sta)>length(sto)
        sto(end+1) = length(running_binary)+1;
    end
    % discard periods of too short running
    discard = (sto - sta) < minRunPeriod/binSizeRun;
    % discard running period for which preceding stationary period is too
    % short
    discard = discard | ...
        sta - [1 sto(1:length(sta)-1)] < noRunBefore/binSizeRun;
    sta(discard) = [];
    sto(discard) = [];
    
    runStarts = timeBins(sta);
    
    for probe = 1:length(db(k).OTprobe)
        rta(count).subject = db(k).subject;
        rta(count).date = db(k).date;
        rta(count).probeName = sp(db(k).OTprobe(probe)).name;
        rta(count).time = timeBins;
        rta(count).running = run;
        rta(count).runStarts = runStarts;
        
        units_OT = db(k).OTunits{probe};
        unit = 1;
        for iCell = 1:length(units_OT)
            if db(k).OTgood{probe}(iCell) == 0
                continue
            end
            fprintf('  Unit %d of %d\n', iCell, length(units_OT))
            rta(count).units(unit).ID = units_OT(iCell);
            % get spiketimes for this cell
            st = sp(db(k).OTprobe(probe)).st(sp(db(k).OTprobe(probe)).clu ...
                == units_OT(iCell));
            st(st<db(k).darkTime(1) | st>db(k).darkTime(2)) = [];
            
            % discard running onsets, which are outside valid times of the
            % unit
            sta = runStarts;
            validTime = db(k).OTtimes{probe}{iCell};
            if ~isempty(validTime)
                invalid = [];
                for s = 1:length(sta)
                    if ~any(sta(s) - noRunBefore > validTime(:,1) & ...
                            sta(s) + minRunPeriod < validTime(:,2))
                        invalid = [invalid s];
                    end
                end
                sta(invalid) = [];
            end
            
            rtSpikes = cell(length(sta),1);
            for s = 1:length(sta)
                t = st - sta(s);
                rtSpikes{s} = t(t>=-noRunBefore & t<=minRunPeriod);
            end
            
            % stopped here
            
            unit = unit + 1;
        end
        set = set + 1;
    end
end


%% OLD

        % RF
%         [~, stats, psth] = sparseNoiseRF(st, stimTimes, stimPositions, ...
%             paramsNoiseRF);
%         pseudoPsths = cell(1,nPseudo);
%         for j = 1:nPseudo
%             [~,~, pseudoPsths{j}] = sparseNoiseRF(st, stimPseudoTimes{j}, ...
%                 stimPositions, paramsNoiseRF);
%         end
%         psth_corr = psth - mean(cat(3, pseudoPsths{:}), 3);
%         [U,S,V] = svd(psth_corr - nanmean(psth_corr(:)), 'econ');
%         spatialRF = reshape(U(:,1), size(stimFrames,1), ...
%             size(stimFrames,2)) .* S(1,1);
%         temporalRF = V(:,1);
%         [~,peakInd] = max(abs(U(:,1)));
%         if U(peakInd,1) < 0
%             spatialRF = -spatialRF;
%             temporalRF = -temporalRF;
%         end

%         % mean PSTH
%         hold on
%         plot(binsNoise, noisePsths{unit}, 'k')
%         plot(binsNoise, noisePsths_smooth{unit}, 'k', 'LineWidth', 2)
%         plot([1 1].*delayNoise(unit), [min(noisePsths{unit}), ...
%             max(noisePsths{unit})], 'r', 'LineWidth', 2)
%         axis tight
%         set(gca, 'box', 'off')
%         ylabel('PSTH')
%         % derivative
%         subplot(plotRows, plotCols, reshape([3;4].*plotCols+(2:3),1,[]))
%         hold on
%         plot(binsNoise(1:end-1), noiseDerivs{unit}, 'k')
%         plot([1 1].*delayNoise(unit), [min(noiseDerivs{unit}), ...
%             max(noiseDerivs{unit})], 'r', 'LineWidth', 2)
%         axis tight
%         set(gca, 'box', 'off')
%         ylabel('Derivative')
%         % time course of RF
%         subplot(plotRows, plotCols, reshape([5;6].*plotCols+(2:3),1,[]))
%         plot(stats.timeBins, temporalRF, 'k')
%         axis tight
%         set(gca, 'box', 'off')
%         xlabel('Time (s)')
%         ylabel('Temporal RF')
%         % spatial RF
%         maxi = max(abs(spatialRF(:)));
%         subplot(plotRows, plotCols, reshape([7;8].*plotCols+(2:3),1,[]))
%         imagesc(noisePosition([1 2]), noisePosition([3 4]), ...
%             spatialRF, [-maxi maxi])
%         colorbar
%         ylabel('Spatial RF')
%         colormap(gca,cm)

%% (1) Depth plots for each recording

% for k = 1:length(db)
%     fprintf('\nProcessing: %s %s\n', db(k).subject, db(k).date);
%     alignDir = fullfile(subjectsFolder, db(k).subject, db(k).date, 'alignments');
%     %% Load spike data
%     sp = loadAllKsDir(db(k).subject, db(k).date);
%     % NOTE! if several electrodes were used for recording, times need to be
%     % aligned! (see scripts\electrophys\extractData
%     
%     [expNums, ~, ~, ~, ~, tl, hasTimeline] = ...
%         dat.whichExpNums(db(k).subject, db(k).date);
%     TLexp = expNums(hasTimeline);
%     tl = tl{1};
%     
%     
%     %% Load and prepare data for flickering screen analysis
%     if isfield(db, 'expFlicker') && ~isempty(db(k).expFlicker)
%         data = load(fullfile(protocolFolder, db(k).subject, db(k).date, ...
%             num2str(db(k).expFlicker), sprintf('%s_%d_%s_parameters.mat', ...
%             db(k).date, db(k).expFlicker, db(k).subject)));
%         parsFlicker = data.parameters.Protocol;
%         
%         bTLtoMaster = readNPY(fullfile(alignDir, ...
%             sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, sp(1).name)));
%         stimOnTL = readNPY(fullfile(alignDir, ...
%             sprintf('mpep_%d_onsets_in_timeline_%d.npy', db(k).expFlicker, TLexp)));
%         stimOffTL = readNPY(fullfile(alignDir, ...
%             sprintf('mpep_%d_offsets_in_timeline_%d.npy', db(k).expFlicker, TLexp)));
%         flickerWhite = applyCorrection(stimOnTL, bTLtoMaster);
%         flickerBlack = applyCorrection(stimOffTL, bTLtoMaster);
%         flickerAll = reshape([flickerWhite, flickerBlack]', [], 1);
%         
%         data = load(fullfile(protocolFolder, db(k).subject, db(k).date, ...
%             num2str(db(k).expFlicker), sprintf('%s_%d_%s_hardwareInfo.mat', ...
%             db(k).date, db(k).expFlicker, db(k).subject)));
%         monitorRefreshRate = data.myScreenInfo.FrameRate;
%         longestFlicker = max(parsFlicker.pars( ...
%             strcmp(parsFlicker.parnames, 'nfr'),:)) / monitorRefreshRate;
%         flickerDurs = diff(flickerAll);
%         flickerStimStartInds = [1; find(flickerDurs > ...
%             longestFlicker * 1.5) + 1; length(flickerAll)+1];
%         flickerTimes = cell(size(parsFlicker.seqnums));
%         for stim = 1:size(flickerTimes,1)
%             for trial = 1:size(flickerTimes,2)
%                 ind = parsFlicker.seqnums(stim,trial);
%                 flickerTimes{stim,trial} = reshape(flickerAll( ...
%                     flickerStimStartInds(ind) : ...
%                     flickerStimStartInds(ind+1)-1), 2, [])';
%             end
%         end
%         flickerDurs = parsFlicker.pars(strcmp(parsFlicker.parnames, ...
%             'nfr'),:) / monitorRefreshRate;
%     end
%     
%     %% Load and prepare data for receptive fields
%     % data = load(fullfile(protocolFolder, subject, date, num2str(exp), ...
%     %     'Protocol.mat'));
%     data = load(fullfile(protocolFolder, db(k).subject, db(k).date, ...
%         num2str(db(k).expNoise), sprintf('%s_%d_%s_parameters.mat', ...
%         db(k).date, db(k).expNoise, db(k).subject)));
%     parsNoise = data.parameters.Protocol;
%     stimFile = str2func(strtok(parsNoise.xfile, '.'));
%     % load myScreenInfo
%     load(fullfile(hardwareFolder, db(k).subject, db(k).date, ...
%         num2str(db(k).expNoise), sprintf('%s_%d_%s_hardwareInfo.mat', ...
%         db(k).date, db(k).expNoise, db(k).subject)));
%     myScreenInfo.windowPtr = NaN;
%     
%     bTLtoMaster = readNPY(fullfile(alignDir, ...
%         sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, sp.name)));
%     stimOnTL = readNPY(fullfile(alignDir, ...
%         sprintf('mpep_%d_onsets_in_timeline_%d.npy', db(k).expNoise, TLexp)));
%     
%     noiseOn = applyCorrection(stimOnTL, bTLtoMaster);
%     
%     % call x-file to create stimuli
%     SS = stimFile(myScreenInfo, parsNoise.pars);
%     stimFrames = cat(3, SS.ImageTextures{:});
%     
%     framesPerImage = parsNoise.pars(6,1);
%     frameTimes = (0 : size(stimFrames, 3)-1) * framesPerImage / myScreenInfo.FrameRate;
%     noiseFrameTimes = (frameTimes + noiseOn)';
%     
%     noisePosition = parsNoise.pars(2:5);
%     
%     %% Load and prepare data for tuning curves
%     data = load(fullfile(protocolFolder, db(k).subject, db(k).date, ...
%         num2str(db(k).expOri), sprintf('%s_%d_%s_parameters.mat', ...
%         db(k).date, db(k).expOri, db(k).subject)));
%     parsGratings = data.parameters.Protocol;
%     
%     bTLtoMaster = readNPY(fullfile(alignDir, ...
%         sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, sp.name)));
%     stimOnTL = readNPY(fullfile(alignDir, ...
%         sprintf('mpep_%d_onsets_in_timeline_%d.npy', db(k).expOri, TLexp)));
%     stimOffTL = readNPY(fullfile(alignDir, ...
%         sprintf('mpep_%d_offsets_in_timeline_%d.npy', db(k).expOri, TLexp)));
%     stimOn = applyCorrection(stimOnTL, bTLtoMaster);
%     stimOff = applyCorrection(stimOffTL, bTLtoMaster);
%     stimDur = mean(stimOff - stimOn);
%     window = [-0.2 stimDur+0.2];
%     binBorders = window(1) : binSizeGrating : window(2);
%     binsGrating = binBorders(1:end-1) + binSizeGrating/2;
%     stimIDs = repmat((1:parsGratings.npfilestimuli)',parsGratings.nrepeats,1);
%     [~,order] = sort(parsGratings.seqnums(:));
%     stimSeq = stimIDs(order);
%     stimBins = binsGrating>0 & binsGrating<stimDur;
%     blank = parsGratings.pars(15,:) == 0;
%     directions = parsGratings.pars(6,:);
%     directions(blank) = [];
%     
%      %% Make plot for whole dataset (recording depth vs. firing rate/visual drive)
%     
%     folder = fullfile(plotFolder, 'plots_basicFeatures_kernelRF', ...
%         [db(k).subject '_' db(k).date]);
%     if ~isdir(folder)
%         mkdir(folder)
%     end
%     
%     minAmplitude = 40;
%     minSpikes = 100;
%     % (1+2) Plot power of spiking in stimulus refresh frequency + mean
%     % firing rate
%     units = unique(sp.clu);
%     units_good = sp.cids(sp.cgs == 2);
%     units_mua = sp.cids(sp.cgs == 1);
%     templates = findTempForEachClu(sp.clu, sp.spikeTemplates);
%     
% %     stimPowers = NaN(length(units), 1);
%     meanRate = NaN(length(units), 1);
%     amplitude = NaN(length(units), 1);
% %     delayNoise = NaN(length(units), 1);
%     delayGratings = NaN(length(units), 1);
%     depths = NaN(length(units), 1);
%     
% %     noisePsths_eachPixel = cell(length(units), 1);
% %     noisePsths = cell(length(units), 1);
% %     noisePsths_smooth = cell(length(units), 1);
% %     noiseDerivs = cell(length(units), 1);
%     gratingPsths = cell(length(units), 1);
%     gratingDerivs = cell(length(units), 1);
%     
%     % determine frequencies for power analysis
%     totalT = sp.st(end)-sp.st(1);
% %     fs = 1000;
% %     timeNoise = noiseOn + (frameTimes(1) : 1/fs : frameTimes(end));
% %     if mod(size(timeNoise,2), 2) > 0
% %         timeNoise(:,end+1) = timeNoise(:,end)+1/fs;
% %     end
% %     fr = (0:(size(timeNoise,2)/2))./size(timeNoise,2).*fs;
% %     frequencies = myScreenInfo.FrameRate / framesPerImage + [0 -.5 .5];
% %     freqInds = zeros(size(frequencies));
% %     for f = 1:length(freqInds)
% %         freqInds(f) = find(fr >= frequencies(f), 1);
% %         if diff(abs(fr(freqInds(f)-[1 0]) - frequencies(f))) > 0
% %             freqInds(f) = freqInds(f)-1;
% %         end
% %     end
%     
%     % get time bins for visual noise
% %     [~, stats] = sparseNoiseRF(sp.st(sp.clu == units(1)), ...
% %         stimTimes, stimPositions, paramsNoisePsth);
% %     binsNoise = stats.timeBins;
% %     bin0 = find(binsNoise < 0, 1, 'last');
% %     bin20 = find(binsNoise < 0.02, 1, 'last');
% %     timeCourses_noise = NaN(length(units), length(binsNoise));
%     
%     % get time bins for gratings
%     maxBin = find(binsGrating > .4, 1);
%     minBin = find(stimBins, 1);
%     timeCourses_gratings = NaN(length(units), length(binsGrating));
%     
%     % loop through all units to get measures
%     fprintf('Getting measures of %d units: ', length(units))
%     for n = 1:length(units)
%         if mod(n,5) == 0
%             fprintf('%d ', n)
%         end
%         % amplitude 
%         amplitude(n) = mean(sp.spikeAmps(sp.clu == units(n)));
%         % firing rate
%         if amplitude(n) >= minAmplitude || ismember(units(n), units_good)
%             meanRate(n) = sum(sp.clu == units(n)) / totalT;
%         end
%         
%         % for the rest, only include units with minimal amplitude and no.
%         % of spikes
%         if (sum(sp.clu == units(n)) < minSpikes || amplitude(n) < minAmplitude) ...
%                 && ~ismember(units(n), units_good)
%             continue
%         end
%         
% %         % power at visual noise frequency (divided by power at neighbouring
% %         % frequencies)
% %         pow = NaN(length(noiseOn), size(timeNoise,2)/2+1);
% %         for trial = 1:length(noiseOn)
% %             ind = sp.clu == units(n) & sp.st > timeNoise(trial,1) & sp.st < timeNoise(trial,end);
% %             s = hist(sp.st(ind), timeNoise(trial,:));
% %             Y = fft(s - mean(s));
% %             P = abs(Y/size(timeNoise,2));
% %             P = P(1:size(timeNoise,2)/2+1);
% %             P(2:end-1) = 2*P(2:end-1);
% %             pow(trial,:) = P;
% %         end
% %         pow = mean(pow,1);
% %         stimPowers(n) = pow(freqInds(1)) / mean(pow(freqInds(2:end)));
%         
% %         % noise delay
% %         [~, ~, noisePsths_eachPixel{n}] = sparseNoiseRF(sp.st(sp.clu == units(n)), ...
% %             stimTimes, stimPositions, paramsNoisePsth);
% %         noisePsths{n} = mean(noisePsths_eachPixel{n},1);
% %         noisePsths_smooth{n} = conv(noisePsths{n} - mean(noisePsths{n}), ...
% %             winHamm, 'same') + mean(noisePsths{n});
% %         m_spont_noise = mean(noisePsths_smooth{n}(bin0 + (-2:1)));
% %         respSquNoise = (noisePsths_smooth{n} - m_spont_noise) .^ 2;
% %         noiseDerivs{n} = medfilt1(diff(respSquNoise), 5);
% %         threshNoise = mean(noiseDerivs{n}(1:bin20)) + ...
% %             thresholdFactorNoise * std(noiseDerivs{n}(1:bin20));
% %         delayIndNoise = find(noiseDerivs{n}(bin0+1:end) > threshNoise, 1) + bin0;
% %         if ~isempty(delayIndNoise)
% %             delayNoise(n) = binsNoise(delayIndNoise);
% %             timeCourses_noise(n,:) = (noisePsths_smooth{n}-noisePsths{n}(bin0)) ./ ...
% %                 max(noisePsths_smooth{n}-noisePsths{n}(bin0));
% %         end
% 
%         % grating delay
%         gratingPsths{n} = NaN(parsGratings.npfilestimuli, length(binsGrating));
%         for s = 1:parsGratings.npfilestimuli
%             ind = find(stimSeq == s);
%             gratingPsths{n}(s,:) = psthAndBA(sp.st(sp.clu == units(n)), ...
%                 stimOn(ind), window, binSizeGrating);
%             gratingPsths{n}(s,:) = conv(gratingPsths{n}(s,:), winHamm, 'same');
%         end
%         m_spont_gratings = mean(gratingPsths{n}(blank,stimBins));
%         respSquGratings = (gratingPsths{n}(~blank,:) - m_spont_gratings).^2;
%         gratingDerivs{n} = medfilt1(diff(mean(respSquGratings,1)), 5);
%         threshGratings = max(minThresh, mean(gratingDerivs{n}(1:minBin)) + ...
%             thresholdFactorGratings * std(gratingDerivs{n}(1:minBin)));
%         delayIndGratings = find(gratingDerivs{n}(minBin:maxBin) > ...
%             threshGratings, 1) + minBin - 1;
%         if ~isempty(delayIndGratings)
%             delayGratings(n) = binsGrating(delayIndGratings);
%             timeCourses_gratings(n,:) = (mean(gratingPsths{n}(~blank,:)) - ...
%                 m_spont_gratings) ./ max(mean(gratingPsths{n}(~blank,:)) - ...
%                 m_spont_gratings);
%         end
%         
%         % depth
%         depths(n) = mean(sp.spikeDepths(sp.clu == units(n)));
%     end
%     fprintf('\n')
%     
%     gu = ismember(units, units_good);
%     mu = ismember(units, units_mua);
%     rest = ~(gu | mu); % not sorted
% %     
% %     cols = lines(5);
% %     ax = zeros(1,5);
% %     markSz = 15;
% %     
% %     figure('Position', [565 42 1356 1074])
% %     subplot(1,5,1)
% %     plot(meanRate(rest), depths(rest), 'k.')
% %     hold on
% %     plot(meanRate(gu)', depths(gu)', '.', 'Color', cols(1,:), 'MarkerSize', markSz)
% %     plot(meanRate(mu)', depths(mu)', 'o', 'Color', cols(1,:))
% %     title('Firing rate')
% %     xlabel('spikes/s')
% %     ylabel('Depth (um)')
% %     set(gca, 'XGrid', 'on', 'YGrid', 'on')
% %     ax(1) = gca;
% %     
% %     subplot(1,5,2)
% %     plot(stimPowers(rest), depths(rest), 'k.')
% %     hold on
% %     plot(stimPowers(gu), depths(gu), '.', 'Color', cols(2,:), 'MarkerSize', markSz)
% %     plot(stimPowers(mu), depths(mu), 'o', 'Color', cols(2,:))
% %     title('Visual drive')
% %     xlabel('stim. / neighb. freq.s')
% %     set(gca, 'XGrid', 'on', 'YGrid', 'on')
% %     ax(2) = gca;
% %     
% %     subplot(1,5,3)
% %     plot(amplitude(rest), depths(rest), 'k.')
% %     hold on
% %     plot(amplitude(gu), depths(gu), '.', 'Color', cols(3,:), 'MarkerSize', markSz)
% %     plot(amplitude(mu), depths(mu), 'o', 'Color', cols(3,:))
% %     title('Spike ampl.')
% %     xlabel('uV')
% %     set(gca, 'XGrid', 'on', 'YGrid', 'on')
% %     ax(3) = gca;
% %     
% %     subplot(1,5,4)
% %     plot(delayNoise(rest), depths(rest), 'k.')
% %     hold on
% %     plot(delayNoise(gu), depths(gu), '.', 'Color', cols(4,:), 'MarkerSize', markSz)
% %     plot(delayNoise(mu), depths(mu), 'o', 'Color', cols(4,:))
% %     title('Resp. delay (noise)')
% %     xlabel('s')
% %     set(gca, 'XGrid', 'on', 'YGrid', 'on', 'XScale', 'log')
% %     mini = min(delayNoise);
% %     maxi = max(delayNoise);
% %     range = maxi - mini;
% %     mini = mini - .05 * range;
% %     maxi = maxi + .05 * range;
% %     xlim([mini maxi])
% %     ax(4) = gca;
% %     
% %     subplot(1,5,5)
% %     plot(delayGratings(rest), depths(rest), 'k.')
% %     hold on
% %     plot(delayGratings(gu), depths(gu), '.', 'Color', cols(5,:), 'MarkerSize', markSz)
% %     plot(delayGratings(mu), depths(mu), 'o', 'Color', cols(5,:))
% %     title('Resp. delay (gratings)')
% %     xlabel('uV')
% %     set(gca, 'XGrid', 'on', 'YGrid', 'on', 'XScale', 'log')
% %     mini = min(delayGratings);
% %     maxi = max(delayGratings);
% %     range = maxi - mini;
% %     mini = mini - .05 * range;
% %     maxi = maxi + .05 * range;
% %     xlim([mini maxi])
% %     ax(5) = gca;
% %     
% %     linkaxes(ax, 'y')
% %     
% %     fig = gcf;
% %     fig.PaperPositionMode = 'auto';
% %     print(fullfile(folder, 'depthPlot.jpg'), '-djpeg','-r0')
% %     
% %     ylim([min(depths(gu))-20 max(depths(gu))+20])
% %     
% %     print(fullfile(folder, 'depthPlot_zoomIn.jpg'), '-djpeg','-r0')
% %     
% %     close(fig)
% end

%% (1) Psths in response to gratings and to visual noise

% for k = 1 %1:length(db)
%     alignDir = fullfile(subjectsFolder, db(k).subject, db(k).date, 'alignments');
%     %% Load spike data
%     sp = loadAllKsDir(db(k).subject, db(k).date);
%     [~, ~, ~, ~, ~, ~, hasTimeline] = ...
%         dat.whichExpNums(db(k).subject, db(k).date);
%     TLexp = find(hasTimeline);
%     TLexp = TLexp(end);
%     
%     
%     %% Load and prepare data for tuning curves
%     data = load(fullfile(protocolFolder, db(k).subject, db(k).date, ...
%         num2str(db(k).expOri), sprintf('%s_%d_%s_parameters.mat', ...
%         db(k).date, db(k).expOri, db(k).subject)));
%     parsGratings = data.parameters.Protocol;
%     
%     bTLtoMaster = readNPY(fullfile(alignDir, ...
%         sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, sp.name)));
%     stimOnTL = readNPY(fullfile(alignDir, ...
%         sprintf('mpep_%d_onsets_in_timeline_%d.npy', db(k).expOri, TLexp)));
%     stimOffTL = readNPY(fullfile(alignDir, ...
%         sprintf('mpep_%d_offsets_in_timeline_%d.npy', db(k).expOri, TLexp)));
%     stimOn = applyCorrection(stimOnTL, bTLtoMaster);
%     stimOff = applyCorrection(stimOffTL, bTLtoMaster);
%     stimDur = mean(stimOff - stimOn);
%     window = [-0.2 stimDur+0.2];
%     binBorders = window(1) : binSizeGrating : window(2);
%     binsGrating = binBorders(1:end-1) + binSizeGrating/2;
%     stimIDs = repmat((1:parsGratings.npfilestimuli)',parsGratings.nrepeats,1);
%     [~,order] = sort(parsGratings.seqnums(:));
%     stimSeq = stimIDs(order);
%     stimBins = binsGrating>0 & binsGrating<stimDur;
%     blank = parsGratings.pars(15,:) == 0;
%     directions = parsGratings.pars(6,:);
%     directions(blank) = [];
%     
%     %% Load and prepare data for receptive fields
%     data = load(fullfile(protocolFolder, db(k).subject, db(k).date, ...
%         num2str(db(k).expNoise), sprintf('%s_%d_%s_parameters.mat', ...
%         db(k).date, db(k).expNoise, db(k).subject)));
%     parsNoise = data.parameters.Protocol;
%     stimFile = str2func(strtok(parsNoise.xfile, '.'));
%     % load myScreenInfo
%     load(fullfile(hardwareFolder, db(k).subject, db(k).date, ...
%         num2str(db(k).expNoise), sprintf('%s_%d_%s_hardwareInfo.mat', ...
%         db(k).date, db(k).expNoise, db(k).subject)));
%     myScreenInfo.windowPtr = NaN;
%     
%     bTLtoMaster = readNPY(fullfile(alignDir, ...
%         sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, sp.name)));
%     stimOnTL = readNPY(fullfile(alignDir, ...
%         sprintf('mpep_%d_onsets_in_timeline_%d.npy', db(k).expNoise, TLexp)));
%     
%     noiseOn = applyCorrection(stimOnTL, bTLtoMaster);
%     
%     % call x-file to create stimuli
%     SS = stimFile(myScreenInfo, parsNoise.pars);
%     stimFrames = cat(3, SS.ImageTextures{:});
%     
%     framesPerImage = parsNoise.pars(6,1);
%     frameTimes = (0 : size(stimFrames, 3)-1) * framesPerImage / myScreenInfo.FrameRate;
%     noiseFrameTimes = (frameTimes + noiseOn)';
%     
%     rows = size(stimFrames,1);
%     cols = size(stimFrames,2);
%     
%     % white or black squares
%     ind = find(stimFrames ~= 0);
%     t = ceil(ind / (rows * cols));
%     ind = mod(ind, (rows * cols));
%     ind(ind == 0) = rows * cols;
%     x = ceil(ind / rows);
%     y = mod(ind, rows);
%     y(y == 0) = rows;
%     time = noiseFrameTimes(t,:);
%     stimPositions = [repmat(y, 3, 1) repmat(x, 3, 1)];
%     stimTimes = reshape(time, [], 1);
%     
%     %% Make plot for each unit
%     
%     folder = fullfile(plotFolder, 'plots_gratingPSTHs', ...
%         [db(k).subject '_' db(k).date]);
%     if ~isdir(folder)
%         mkdir(folder)
%     end
%     
%     units = unique(sp.clu);
%     col = lines(3);
%     
%     [~, stats] = sparseNoiseRF(sp.st(sp.clu == units(1)), ...
%         stimTimes, stimPositions, params);
%     binsNoise = stats.timeBins;
%     bin0 = find(binsNoise < 0, 1, 'last');
%     bin20 = find(binsNoise < 0.02, 1, 'last');
%     
%     maxBin = find(binsGrating>.4, 1);
%     minBin = find(stimBins, 1);
%     
%     % loop through all units to get measures
%     for n = 1:length(units)
%         % grating delay
%         psth = NaN(parsGratings.npfilestimuli, length(binsGrating));
%         for s = 1:parsGratings.npfilestimuli
%             ind = find(stimSeq == s);
%             [psth(s,:),~,~,~,~,spikeCount] = psthAndBA(sp.st(sp.clu == units(n)), ...
%                 stimOn(ind), window, binSizeGrating);
%             psth(s,:) = conv(psth(s,:), winHamm, 'same');
%             if s == find(blank)
%                 respEachBlank = spikeCount ./ binSizeGrating;
%             end
%         end
%         
%         m_spont = mean(psth(blank,stimBins));
%         se_spont = std(respEachBlank,0,1) ./ sqrt(size(respEachBlank,1));
%         se_spont = mean(se_spont);
%         respSqu = (psth(~blank,:) - m_spont).^2;
%         peak = max(mean(respSqu(:,minBin:maxBin),1));
%         deriv = medfilt1(diff(mean(respSqu,1)), 5);
%         thresh = max(minThresh, mean(deriv(1:minBin)) + ...
%             thresholdFactorGratings * std(deriv(1:minBin)));
%         delayInd = find(deriv(minBin:maxBin) > thresh, 1) + minBin - 1;
%         
%         ax = zeros(1,3);
% 
%         figure('Position', [3 446 1914 668])
%         subplot(3,3,1:2)
%         hold on
%         plot(binsGrating, psth(~blank,:), 'Color', [1 1 1].*.5)
%         plot(binsGrating, mean(psth(~blank,:),1), 'k', 'LineWidth', 2)
%         plot(binsGrating, psth(blank,:), ':', 'Color', col(2,:), 'LineWidth', 2)
%         h = plot(binsGrating([1 end]), [1 1].*m_spont, 'Color', col(2,:), 'LineWidth', 2);
%         plot(binsGrating([1 end]), [1 1].*(m_spont+1*se_spont), 'Color', col(2,:))
%         plot(binsGrating([1 end]), [1 1].*(m_spont-1*se_spont), 'Color', col(2,:))
%         data = [psth(:); m_spont+se_spont; m_spont-se_spont];
%         if ~isempty(delayInd)
%             plot([1 1].*binsGrating(delayInd), [min(data) max(data)], 'r', ...
%                 'LineWidth', 2)
%         end
%         legend(h, 'blank response')
%         ylabel('Firing rate')
%         title('Response to drifting gratings')
%         set(gca, 'box', 'off', 'XTick', -.15:.05:1)
%         axis tight
%         ax(1) = gca;
%         
%         subplot(3,3,4:5)
%         hold on
% %         plot(binsGrating, respSqu, 'Color', [1 1 1].*.5)
%         plot(binsGrating, mean(respSqu,1), 'k', 'LineWidth', 2)
%         h = plot(binsGrating([1 end]), [1 1].*peak/2, 'Color', col(1,:), 'LineWidth', 2);
%         if ~isempty(delayInd)
%             plot([1 1].*binsGrating(delayInd), [min(mean(respSqu,1)) ...
%                 max(mean(respSqu,1))], 'r', 'LineWidth', 2)
%         end
%         legend(h, 'half peak')
%         ylabel('Squared fir. rate')
%         set(gca, 'box', 'off', 'XTick', -.15:.05:1)
%         axis tight
%         ax(2) = gca;
%         
%         subplot(3,3,7:8)
%         hold on
%         plot(binsGrating(1:end-1), deriv, 'k')
%         plot(binsGrating([1 end-1]), [1 1].*thresh, 'r')
%         if ~isempty(delayInd)
%             h = plot([1 1].*binsGrating(delayInd), [min(deriv) max(deriv)], 'r', ...
%                 'LineWidth', 2);
%             legend(h, 'delay (thresh. crossed)')
%         end
%         ylabel('1st derivative')
%         xlabel('Time after stim. onset (s)')
%         set(gca, 'box', 'off', 'YGrid', 'on', 'XTick', -.15:.05:1)
%         axis tight
%         ax(3) = gca;
%         
%         linkaxes(ax, 'x')
%         xlim([binsGrating(1) 0.4])
%         
%         
%         % visual noise delay
%         [~, stats, psth] = sparseNoiseRF(sp.st(sp.clu == units(n)), ...
%             stimTimes, stimPositions, params);
%         respMean = mean(psth,1);
%         respMean = conv(respMean - mean(respMean), winHamm, 'same') + mean(respMean);
%         m_spont = mean(respMean(bin0 + (-2:1)));
%         se_spont = std(psth(:,bin0 + (-2:1)),0,1) ./ sqrt(size(psth,1));
%         se_spont = mean(se_spont);
%         respSqu = (respMean - m_spont) .^ 2;
%         peak = max(respSqu(bin0+1:end));
%         deriv = diff(respSqu);
%         derivFilt = medfilt1(deriv, 5);
%         thresh = mean(derivFilt(1:bin20)) + ...
%             thresholdFactorNoise * std(derivFilt(1:bin20));
%         delayInd = find(derivFilt(bin0+1:end) > thresh, 1) + bin0;
%         
%         SVDtime = stats.timeCourse ./ sum(abs(stats.timeCourse - ...
%             mean(stats.timeCourse))) .* sum(abs(respMean - mean(respMean)));
%         SVDtimeFilt = conv(SVDtime - mean(SVDtime), winHamm, 'same') + mean(SVDtime);
%         SVD_spont = mean(SVDtimeFilt(bin0 + (-2:1)));
%         SVDsqu = (SVDtimeFilt - SVD_spont) .^ 2;
%         SVDderiv = diff(SVDsqu);
%         SVDderivFilt = medfilt1(SVDderiv, 5);
%         
%         ax = zeros(1,3);
%         
%         subplot(3,3,3)
%         hold on
%         plot(binsNoise, SVDtime + mean(respMean), 'Color', [1 1 1].*.5)
%         plot(binsNoise, SVDtimeFilt + mean(respMean), 'Color', [1 1 1].*.5, 'LineWidth', 2)
%         plot(binsNoise,mean(psth),'k')
%         plot(binsNoise, respMean, 'k', 'LineWidth', 2)
%         h = plot(binsNoise([1 end]), [1 1].*m_spont, 'Color', col(2,:), 'LineWidth', 2);
%         plot(binsNoise([1 end]), [1 1].*(m_spont+se_spont), 'Color', col(2,:))
%         plot(binsNoise([1 end]), [1 1].*(m_spont-se_spont), 'Color', col(2,:))
%         data = [mean(psth)'; respMean(:); m_spont+se_spont; m_spont-se_spont; ...
%             SVDtime(:) + mean(respMean)];
%         if ~isempty(delayInd)
%             plot([1 1].*binsNoise(delayInd), [min(data) max(data)], 'r', ...
%                 'LineWidth', 2)
%         end
%         legend(h, 'spont. act.')
%         title('Response to visual noise')
%         set(gca, 'box', 'off')
%         axis tight
%         ax(1) = gca;
%         
%         subplot(3,3,6)
%         hold on
%         plot(binsNoise, SVDsqu, 'Color', [1 1 1].*.5, 'LineWidth', 2)
%         plot(binsNoise, respSqu, 'k', 'LineWidth', 2)
%         h = plot(binsNoise([1 end]), [1 1].*peak/2, 'Color', col(1,:), ...
%             'LineWidth', 2);
%         if ~isempty(delayInd)
%             plot([1 1].*binsNoise(delayInd), [min(respSqu) max(respSqu)], 'r', ...
%                 'LineWidth', 2)
%         end
%         legend(h, 'half peak')
%         set(gca, 'box', 'off')
%         axis tight
%         ax(2) = gca;
%         
%         subplot(3,3,9)
%         hold on
%         plot(binsNoise(1:end-1), SVDderiv, 'Color', [1 1 1].*.5)
%         plot(binsNoise(1:end-1), SVDderivFilt, 'Color', [1 1 1].*.5, 'LineWidth', 2)
%         plot(binsNoise(1:end-1), deriv, 'k')
%         plot(binsNoise(1:end-1), derivFilt, 'k', 'LineWidth', 2)
%         plot(binsNoise([1 end-1]), [1 1].*thresh, 'r')
%         if ~isempty(delayInd)
%             plot([1 1].*binsNoise(delayInd), [min(deriv) max(deriv)], 'r', ...
%                 'LineWidth', 2)
%         end
%         xlabel('Time after stim. onset (s)')
%         set(gca, 'box', 'off', 'YGrid', 'on')
%         axis tight
%         ax(3) = gca;
%         
%     
%         fig = gcf;
%         fig.PaperPositionMode = 'auto';
%         print(fullfile(folder, sprintf('neuron%04d.jpg', ...
%             units(n))), '-djpeg','-r0')
%         close(fig)
%     end
% end




% Load and prepare data for running correlation

%     rotary = double(tl.rawDAQData(:, strcmp({tl.hw.inputs.name}, ...
%         'rotaryEncoder')));
%     tlTime = applyCorrection(tl.rawDAQTimestamps, bTLtoMaster);
%     runSpeed = nonVis.getRunningSpeed_wheel(rotary, tlTime, runningSigma);
%     runningTime = runSpeed.t;
%     runSpeed = runSpeed.total;
%     cmPerUnit = 2*pi * 8.75 / (4 * 1024);
%     runSpeed = runSpeed * cmPerUnit;
%     binSizeRun = 0.001;
%     numD = size(db(k).darkTime,1);
%     binsRun = cell(1, numD);
%     time = cell(1, numD);
%     running = cell(1, numD);
%     for d = 1:numD
%         binsRun{d} = db(k).darkTime(d,1) : binSizeRun : db(k).darkTime(d,2);
%         time{d} = binsRun{d}(1:end-1) + binSizeRun/2;
%         running{d} = interp1(runningTime, runSpeed, time{d}, 'pchip');
%         time{d}(end+1) = NaN;
%         if d>1
%             time{d} = time{d} - (time{d}(1)-time{d-1}(end-1)) + 100;
%         end
%     end
%     time = cat(2, time{:});
%     time(end) = [];