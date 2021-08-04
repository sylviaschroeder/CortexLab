%% Parameters
samplingRate = 2500; % in ....imec.lf.meta file
numChans = 385;
stimWin1 = [0 .1];
stimWin2 = [0 .1];
baseWin = [-.05 0];

%% Folders
folderData = '\\zubjects.cortexlab.net\Subjects';
folderSave = 'C:\STORAGE\OneDrive - University of Sussex\Projects\2021_Joanna_competition in SC';
folderPlots = fullfile(folderSave, 'Results\LFP plots');

folderScript = 'C:\dev\workspace\CortexLab';
folderTools = 'C:\STORAGE\workspaces';

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')));
addpath(genpath(fullfile(folderTools, 'spikes')));
addpath(genpath(fullfile(folderScript)));

%% Define datasets

% subject = 'SS093';
% date = '2018-05-24';
% probe = 'K1';
db = db_ephys_task;
for k = 1:length(db)
    subject = db(k).subject;
    date = db(k).date;
    fprintf('Dataset %d of %d: %s %s\n', k, length(db), subject, date)
        
    %% Load stimulus data
    folderAlf = fullfile(folderData, subject, date, 'alf');
    noise = io.getVisNoiseInfo(folderAlf);
    noiseStim = (noise.frames(:,:,noise.stimOrder) + 1) ./ 2;
    noiseDeriv = diff(cat(3, zeros(size(noiseStim,1), size(noiseStim,2)), ...
        noiseStim), 1, 3); % 1 if switching to white from black, -1 if switching to black from white, 0, if value does not change
    noiseDeriv = reshape(noiseDeriv, [], size(noiseStim,3));
    
    for pr = 1:length(db(k).probes)
        probe = db(k).probes{pr};
        fprintf('  Probe %d (%s)\n', pr, probe)
        
        %% Get channel coordinates
        x = repmat([43; 11; 59; 27], 96, 1);
        y = reshape(repmat(20:20:3840, 2, 1), [], 1);
        channelCoords = [x, y];
        
        %% Load LFP data
        folderAlign = fullfile(folderData, subject, date, 'alignments');
        folderLFP = fullfile(folderData, subject, date, sprintf('ephys_%s', probe));
        fileLFP = dir(fullfile(folderLFP, '*.lf.bin'));
        
        lfp = memmapfile(fullfile(folderLFP, fileLFP.name), ...
            'Format',  {'int16', [numChans fileLFP.bytes/(numChans*2)], 'x'});
        t_lfp = (1:fileLFP.bytes/(numChans*2))' ./ samplingRate;
        
        % correct time to align with master if necessary
        fileAlign = dir(fullfile(folderAlign, sprintf('correct_ephys_%s_*.npy', probe)));
        if ~isempty(fileAlign)
            probeToMaster = readNPY(fullfile(folderAlign, fileAlign.name));
            t_lfp = applyCorrection(t_lfp, probeToMaster);
        end
        
        start = find(t_lfp > noise.times(1)-5, 1);
        stop = find(t_lfp > noise.times(end)+1, 1);
        
        lfp = lfp.Data.x(1:end-1,start:stop); % last channel is external input -> not LFP
        t_lfp = t_lfp(start:stop);
        
        % find reference/noise channels
        lfpStd = std(double(lfp(:,1:100000)), 0, 2);
        [stdSorted,order] = sort(lfpStd);
        [~,jumpInStd] = max(diff(stdSorted(1:100)));
        noiseChans = sort(order(1:jumpInStd));
        
        % ignore noise channels
        lfp(noiseChans,:) = [];
        channelCoords(noiseChans,:) = [];
        
        %% Get stimulus triggered LFP + most driving pixel
        winSamples = stimWin1(1)*samplingRate : stimWin1(2)*samplingRate;
        baseSamples = baseWin(1)*samplingRate : baseWin(2)*samplingRate;
        lfpEvoked = zeros(size(lfp,1), size(noiseDeriv,1), length(winSamples)); % [channels x pixels x time]
        
        for pix = 1:size(noiseDeriv,1)
            stimEv = noiseDeriv(pix,:) ~= 0;
            t = noise.times(stimEv);
            [~,~,t_ind] = histcounts(t, [0; t_lfp+1/(2*samplingRate)]);
            base = lfp(:,t_ind + baseSamples);
            base = reshape(base, size(lfp,1), length(t_ind), length(baseSamples));
            base = mean(double(base), 3);
            ev = lfp(:, t_ind + winSamples);
            ev = reshape(ev, size(lfp,1), length(t_ind), length(winSamples));
            lfpEvoked(:,pix,:) = mean(double(ev) - base, 2);
        end
        
        % find best driving pixel and time of max response
%         [~, maxChan] = max(max(max(abs(lfpEvoked), [], 3), [], 2), [], 1);
        [~, maxChan] = min(min(min(lfpEvoked, [], 3), [], 2), [], 1);
        [~, maxT] = min(min(lfpEvoked(maxChan,:,:), [], 2), [], 3);
        [~, maxPix] = min(lfpEvoked(maxChan,:,maxT), [], 2);
        maxRF = reshape(squeeze(lfpEvoked(maxChan,:,maxT)), size(noiseStim,1), size(noiseStim,2));
        maxTimeCourse = smooth(squeeze(lfpEvoked(maxChan, maxPix, :)), 20);
        [~, maxT] = min(maxTimeCourse);
        
        %% Find upper border of SC (according to: Ito, Feldheim, Litke, 2017)
        % determine amplitude profile for right and left channels
        chans_right = find(channelCoords(:,1) > 30);
        ampProfile_right = lfpEvoked(chans_right,maxPix,maxT);
        chans_left = find(channelCoords(:,1) < 30);
        ampProfile_left = lfpEvoked(chans_left,maxPix,maxT);
        
        % find borders on right and left channels
        ampProfSm_right = smooth(ampProfile_right, 5);
        [minR, ind] = min(ampProfSm_right);
        ind2 = find(ampProfSm_right(ind+1:end) > minR/20, 1);
        if isempty(ind2)
            fprintf('  Something wrong with LFP. Amplitudes above channel with minimum do not approach zero.\n')
            continue
        end
        border_right = interp1(ampProfSm_right(ind : ind+ind2), ...
            channelCoords(chans_right(ind : ind+ind2),2), minR/2);
        ampProfSm_left = smooth(ampProfile_left, 5);
        [minL, ind] = min(ampProfSm_left);
        ind2 = find(ampProfSm_left(ind+1:end) > minL/20, 1);
        if isempty(ind2)
            fprintf('    Something wrong with LFP. Amplitudes above channel with minimum do not approach zero.')
            continue
        end
        border_left = interp1(ampProfSm_left(ind : ind+ind2), ...
            channelCoords(chans_left(ind : ind+ind2),2), minL/2);
        upperBorder = mean([border_right border_left]);
        
        % save location of upper border to data
        f = fullfile(folderSave, 'Data', subject, date, probe);
        if ~isfolder(f)
            mkdir(f)
        end
        writeNPY(upperBorder, fullfile(f, '_ss_colliculusTop.depth.npy'));
        
        % plot amplitude profile
        h = [0 0];
        figure
        hold on
        plot(ampProfile_right, channelCoords(chans_right,2), 'b')
        h(1) = plot(ampProfSm_right, channelCoords(chans_right,2), ...
            'b', 'LineWidth', 2);
        plot(minR/2, border_right, 'ob', 'LineWidth', 2)
        plot(ampProfile_left, channelCoords(chans_left,2), 'r')
        h(2) = plot(ampProfSm_left, channelCoords(chans_left,2), ...
            'r', 'LineWidth', 2);
        plot(minL/2, border_left, 'or', 'LineWidth', 2)
        ylim([upperBorder - 1000, min(upperBorder + 500, max(channelCoords(:,2)))])
        legend(h, 'right channels','left channels', 'Location', 'NorthWest')
        xlabel('LFP amplitude')
        ylabel('Depth on probe (\mum)')
        title(sprintf('Top edge of SC (%s %s %s)', subject, date, probe))
        saveas(gcf, fullfile(folderPlots, ...
            sprintf('%s_%s_%s_amplitudeProfile.jpg', subject, date, probe)))
        close gcf
        
        %% Plot evoked LFP
        % merge signal from same y-position
        yPos = unique(channelCoords(:,2));
        lfpMerged = zeros(length(yPos), length(winSamples));
        for y = 1:length(yPos)
            ind = channelCoords(:,2) == yPos(y);
            lfpMerged(y,:) = mean(lfpEvoked(ind,maxPix,:),1);
        end
        
        % interpolate evoked LFP to equidistance depths on probe
        winTime = winSamples ./ samplingRate;
        newYPos = yPos(1):20:yPos(end);
        if any(yPos ~= newYPos)
            [X,Y] = meshgrid(winTime, yPos);
            [X2, depths] = meshgrid(winTime, yPos(1):20:yPos(end));
            lfpMerged = interp2(X, Y, lfpMerged, X2, depths);
            yPos = newYPos;
        end
        
        % plot evoked LFP
        maxi = max(abs(lfpMerged(:)));
        figure('Position', [950 110 560 870])
        hold on
        imagesc(winTime([1 end]), yPos([1 end]), lfpMerged, [-maxi maxi])
        set(gca, 'YDir', 'normal')
        plot(winTime([1 end]), [1 1] .* upperBorder, 'k', 'LineWidth', 1)
        xlim(winTime([1 end]))
        ylim([upperBorder - 500 min(upperBorder + 200, max(channelCoords(:,2)))])
        colormap(colmaps.getBlueWhiteRedMap(201))
        c = colorbar;
        c.Label.String = 'mV';
        xlabel('Time from frame onset (s)')
        ylabel('Probe depth (\mum)')
        title(sprintf('Evoked LFP (%s %s %s)', subject, date, probe))
        saveas(gcf, fullfile(folderPlots, ...
            sprintf('%s_%s_%s_evokedLFP.jpg', subject, date, probe)))
        close gcf
        
        %% Current-source density analysis to determine border between intermediate and deep SC
        % according to Stitt, ..., Engel (2013), and Lee, Tran, Turan, Meister
        % (2020)
        winSamples = stimWin2(1)*samplingRate : stimWin2(2)*samplingRate;
        winTime = winSamples ./ samplingRate;
        stimEv = noiseDeriv(maxPix,:) ~= 0;
        t = noise.times(stimEv);
        [~,~,t_ind] = histcounts(t, [0; t_lfp+1/(2*samplingRate)]);
        base = lfp(:,t_ind + baseSamples);
        base = reshape(base, size(lfp,1), length(t_ind), length(baseSamples));
        base = mean(double(base), 3);
        ev = lfp(:, t_ind + winSamples);
        ev = double(reshape(ev, size(lfp,1), length(t_ind), length(winSamples))) - base;
        ev = flip(ev,1);
        
        % do CSD for each x-position separately
        xPos = unique(channelCoords(:,1));
        CSD = cell(1, length(xPos));
        for x = 1:length(xPos)
            ind = channelCoords(:,1) == xPos(x);
            yPos = channelCoords(ind,2);
            ySteps1 = diff(yPos);
            yPos = yPos(1:end-1) + 0.5 .* ySteps1;
            ySteps2 = diff(yPos);
            yPos = flip(yPos(1:end-1) + 0.5 .* ySteps2);
            ySteps1 = ySteps1 ./ 1000; % in mm
            ySteps2 = ySteps2 ./ 1000; % in mm
            csdX = squeeze(mean(diff(diff(ev(ind,:,:) ./ 1000, 1,1) ./ ySteps1, ...
                1,1) ./ ySteps2, 2)); % dividing by spatial distance (twice) is important
            % as distances vary (if channel is missing)
            % unit: mV / mm^2
            CSD{x} = [yPos csdX];
        end
        
        % plot CSDs for each x-position
        figure('Position', [680 220 1220 760])
        hold on
        for x = 1:length(xPos)
            plot(winTime + (stimWin2(2)+.1).*(x-1), ...
                CSD{x}(:,2:end) + CSD{x}(:,1))
        end
        ylim([upperBorder - 500 min(upperBorder + 200, max(channelCoords(:,2)))])
        set(gca, 'XTick', stimWin2)
        xlabel('Time from frame onset (s)')
        ylabel('Probe depth (\mum)')
        title(sprintf('CSD for chans at each x-position along probe (%s %s %s)', subject, date, probe))
        saveas(gcf, fullfile(folderPlots, ...
            sprintf('%s_%s_%s_CSD-perX.jpg', subject, date, probe)))
        close gcf
        
        % merge CSDs
        merged = cat(1,CSD{:});
        yPos = flip(unique(merged(:,1)));
        CSDmerged = zeros(length(yPos), length(winSamples));
        for y = 1:length(yPos)
            ind = merged(:,1) == yPos(y);
            CSDmerged(y,:) = mean(merged(ind,2:end),1);
        end
        
        % smooth and interpolate CSD to equidistance depths on probe
        [X,Y] = meshgrid(winTime, yPos);
        [X2, depths] = meshgrid(winTime, flip(yPos(end):10:yPos(1)));
        CSDinterp = interp2(X, Y, CSDmerged, X2, depths);
        win = normpdf(-3:3, 0, 1);
        c = convmtx(win' , size(CSDinterp,1));
        CSDinterp = c * CSDinterp;
        CSDinterp([1:3 end-2:end],:) = [];
        CSDinterp(2:2:end,:) = [];
        depths = depths(:,1);
        depths(2:2:end) = [];
        
        % find border between SGS (superficial cell layer) and SO (stratum opticum)
        avgCSDs = mean(CSDinterp, 2);
        [maxi, indMax] = max(avgCSDs);
        [mini, indMin] = min(avgCSDs(indMax+1:end));
        SGS_SO_border = depths(indMax) + (maxi / (maxi - mini)) * ...
            (depths(indMax + indMin) - depths(indMax));
        
        % plot interpolated CSD
        maxi = max(abs(CSDinterp(:)));
        figure('Position', [950 110 560 870])
        hold on
        imagesc(winTime([1 end]), depths([1 end]), CSDinterp, [-maxi maxi])
        set(gca, 'YDir', 'normal')
        plot(winTime([1 end]), [1 1] .* upperBorder, 'k', 'LineWidth', 1)
        plot(winTime([1 end]), [1 1] .* SGS_SO_border, 'k', 'LineWidth', 1)
        plot(0.01, depths(indMax), 'or', 'LineWidth', 2)
        plot(0.01, depths(indMax + indMin), 'ob', 'LineWidth', 2)
        xlim(winTime([1 end]))
        ylim([upperBorder - 500 min(upperBorder + 200, max(channelCoords(:,2)))])
        colormap(colmaps.getBlueWhiteRedMap(201))
        c = colorbar;
        c.Label.String = 'mV/mm^2';
        xlabel('Time from frame onset (s)')
        ylabel('Probe depth (\mum)')
        title(sprintf('CSD merged and interpolated (%s %s %s)', subject, date, probe))
        saveas(gcf, fullfile(folderPlots, ...
            sprintf('%s_%s_%s_CSD.jpg', subject, date, probe)))
        close gcf
    end
end