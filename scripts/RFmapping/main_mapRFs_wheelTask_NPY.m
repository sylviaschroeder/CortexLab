%% Folders
folderBase = 'C:\Users\Sylvia\OneDrive - University of Sussex\Lab\DATA\NPY\task_2p';
folderTools = 'C:\dev\toolboxes';
folderThisRepo = 'C:\dev\workspaces\CortexLab';
folderResults = 'C:\Users\Sylvia\OneDrive - University of Sussex\Lab\RESULTS\wheelTask\RFs_new';

%% Parameters
% for correcting baseline drifts of calcium traces at start of experiments
driftWin = 20; % in s, window to test whether baseline is higher than normal
driftThresh = 1.5; % in std, threshold for drift
correctWin = 150; % in s, window to fit exponential

% for receptive field estimates
% used for fitting 2 RFs (ON and OFF simultaneously), and fitting running
% kernels and RFs simultaneously
% lambdasStim = logspace(-5, 0, 6);
lambdasStim = logspace(-5, -2, 4);
RFlimits = [0.2 0.4];
% RFlimits = [0.2 1.0];
crossFolds = 10;

% for evaluation of receptive fields (significance/goodness)
minExplainedVarianceStim = 0.01;
minPVal = 0.05;

% colormaps
red = [1 0 .5];
blue = [0 .5 1];
black = [1 1 1].*0.5;
grad = linspace(0,1,100)';
reds = red.*flip(grad) + [1 1 1].*grad;
blacks = black.*flip(grad) + [1 1 1].*grad;
cm_ON = [blacks; flip(reds(1:end-1,:),1)];
blues = blue.*flip(grad) + [1 1 1].*grad;
cm_OFF = [blacks; flip(blues(1:end-1,:),1)];
colormaps = {cm_ON, cm_OFF};

titles = {'ON field','OFF field'};

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(fullfile(folderThisRepo))

%% Fit RFs and get cross-validated explained variance

subjects = dir(fullfile(folderBase, 'SS*'));
for subj = 1:length(subjects)
    name = subjects(subj).name;
    dates = dir(fullfile(folderBase, name, '2*'));
    for dt = 1:length(dates)
        date = dates(dt).name;
        folder = fullfile(folderBase, name, date);

        if ~isfile(fullfile(folder, '_ss_sparseNoise.times.npy'))
            continue
        end
        
        fResults = fullfile(folderResults, name, date);
        if ~isfolder(fResults)
            mkdir(fResults)
        end
        
        % load data
        traces = readNPY(fullfile(folder, '_ss_2pCalcium.dff.npy'));
        time = readNPY(fullfile(folder, '_ss_2pCalcium.timestamps.npy'));
        delays = readNPY(fullfile(folder, '_ss_2pPlanes.delay.npy'));
        planes = readNPY(fullfile(folder, '_ss_2pRois._ss_2pPlanes.npy'));
        stimTimes = readNPY(fullfile(folder, '_ss_sparseNoise.times.npy'));
        stimMaps = readNPY(fullfile(folder, '_ss_sparseNoiseID.map.npy'));
        stimSeq = readNPY(fullfile(folder, '_ss_sparseNoise._ss_sparseNoiseID.npy'));
        stimPos = readNPY(fullfile(folder, '_ss_sparseNoiseArea.edges.npy')); % [left right top bottom]
        % flip sign of stimulus borders along y-axis -> positive numbers
        % are above horizon/monitor centre
        stimPos(3:4) = -stimPos(3:4);

        % if stimulus covers both hemifields (SS052, SS057, SS058), only
        % consider right hemifield (affected datasets all recorded in left
        % hemisphere)
        squW = diff(stimPos(1:2)) / size(stimMaps,3);
        squH = -diff(stimPos(3:4)) / size(stimMaps,2);
        if stimPos(1) * stimPos(2) < 0
            % determine left edge of all pixel columns
            leftEdges = stimPos(1) + (0:size(stimMaps,3)-1) .* squW;
            validPix = leftEdges >= 0;
            stimMaps = stimMaps(:,:,validPix);
            stimPos(1) = leftEdges(find(validPix,1));
        end
        
        % interpolate calcium traces to align all to same time
        t_ind = time > stimTimes(1) - 10 & time < stimTimes(end) + 10;
        traces = traces(t_ind,:);
        time = time(t_ind);
        timeBin = median(diff(time));
        delta_t = median(diff(delays));
        upsample = round(timeBin / delta_t);
        timeBin = timeBin / upsample;
        t = reshape((time + (0:upsample-1) * timeBin)', [], 1);
        tr_up = NaN(length(t), size(traces,2));
        for d = 1:length(delays)
            indUnits = find(planes == d);
            for n = indUnits'
                if all(isnan(traces(:,n)))
                    continue
                end
                nanInd1 = isnan(traces(:,n));
                tr_up(:,n) = interp1(time(~nanInd1) + delays(d), ...
                    traces(~nanInd1,n), t, 'pchip');
                nanInd2 = reshape(repmat(nanInd1, 1, upsample), [], 1);
                tr_up(nanInd2,n) = NaN;
            end
        end
        
        % remove strong baseline decay at start of experiment in cells that
        % show it
        indUnits = find(nanmean(tr_up(1:round(driftWin / timeBin),:),1) > ...
            nanmean(tr_up,1) + driftThresh .* nanstd(tr_up,0,1));
        ind = round(correctWin / timeBin);
        for iUnit = 1:length(indUnits)
            y = tr_up(:, indUnits(iUnit));
            y = fillmissing(y, 'linear');
            % fit double exponential to start of trace
            f = fit((1:length(y))', y, ...
                @(a,b,c,d,e,x) a + b .* exp(-x ./ c) + d .* exp(-x ./ e), ...
                'Lower', [0 0 0 0 0], ...
                'Upper', [max(y) max(y) 500 max(y) 500], ...
                'StartPoint', [min(y) mean(y) 50 mean(y) 5]);
            % remove fit
            tr_up(:, indUnits(iUnit)) = y - f(1 : size(tr_up,1)) + f.a;
        end
        
        % more stimulus information
        stimFrames = stimMaps(stimSeq,:,:);
        stimFrameDur = median(diff(stimTimes));
        RFtimesInFrames = floor(RFlimits(1) / stimFrameDur) : ...
            ceil(RFlimits(2) / stimFrameDur);
        clear stimMaps stimSeq traces

        % map RF
        [rFields, ev, pred, t_pred] = ...
            whiteNoise.getReceptiveField( ...
            tr_up, t, stimFrames, stimTimes, ...
            RFtimesInFrames, lambdasStim, crossFolds);
        
        v = mean(ev,3); % [neuron x lamStim], average across cross-folds
        [maxEV, maxStimLam] = max(v,[],2);

        % test signficance of each RF
        [ev, ev_shift] = ...
            whiteNoise.receptiveFieldShiftTest( ...
            tr_up, t, stimFrames, stimTimes, ...
            RFtimesInFrames, rFields, maxStimLam, 500);
        pvals = sum(ev_shift > ev, 2) ./ size(ev_shift,2);
        pvals(isnan(ev)) = NaN;

%         % ONLY FOR REFITTING OF PARAMETERS
%         rFields = readNPY(fullfile(folder, '_ss_rf.maps.npy'));
%         rFields = permute(rFields, [2 3 4 5 1]);
%         maxEV = readNPY(fullfile(folder, '_ss_rf.explVars.npy'));
%         pvals = readNPY(fullfile(folder, '_ss_rf.pValues.npy'));
%         maxStimLam = readNPY(fullfile(folder, '_ss_rf.lambdasStim.npy'));

        % fit Gaussian
        fitPars = NaN(size(tr_up,2), 6); % each row: [amplitude, xCenter, 
                                         % xStd, yCenter, yStd, rotation]
        for iUnit = 1:size(fitPars,1)
            if maxEV(iUnit) < minExplainedVarianceStim || pvals(iUnit) > minPVal
                continue
            end
            rf = rFields(:,:,:,:,iUnit); % [rows x cols x t x ON/OFF]
            rf(:,:,:,2) = -rf(:,:,:,2);
            [~,mxTime] = max(max(reshape(permute(abs(rf), [1 2 4 3]), ...
                [], size(rf,3)), [], 1));
            rf = squeeze(rf(:,:,mxTime,:));
            r = reshape(rf,[],2);
            [~,mxPos] = max(max(abs(r),[],2));
            if r(mxPos,1) < 0
                rf(:,:,1) = -rf(:,:,1);
            end
            if r(mxPos,2) < 0
                rf(:,:,2) = -rf(:,:,2);
            end
            rf = mean(rf,3);
            % interpolate RF so that pixels are square with edge length of 
            % 1 visual degree
            [x0, y0] = meshgrid(1:size(rf,2), 1:size(rf,1));
            x1 = linspace(0.5, size(rf,2)+0.5, diff(stimPos(1:2)));
            x1(x1<1 | x1>size(rf,2)) = [];
            y1 = linspace(0.5, size(rf,1)+0.5, -diff(stimPos(3:4)));
            y1(y1<1 | y1>size(rf,1)) = [];
            [x1, y1] = meshgrid(x1, y1);
            rf = interp2(x0, y0, rf, x1, y1);
            % fit Gaussian
            fitPars(iUnit,:) = whiteNoise.fit2dGaussRF(rf, false);
        end
        % transform RF position to absolute values (relative to screen)
        fitPars(:,2) = stimPos(1) + squW/2 + fitPars(:,2);
        fitPars(:,4) = stimPos(3) - squH/2 - fitPars(:,4);
        % mirror RF orientation to account for flipped y-axis direction
        % (top is positive)
        fitPars(:,6) = -fitPars(:,6);
        
        % save results
        writeNPY(permute(rFields, [5 1 2 3 4]), fullfile(folder, '_ss_rf.maps.npy'))
        writeNPY(maxEV, fullfile(folder, '_ss_rf.explVars.npy'))
        writeNPY(maxStimLam, fullfile(folder, '_ss_rf.lambdasStim.npy'))
        writeNPY(stimPos, fullfile(folder, '_ss_rfDescr.edges.npy'))
        writeNPY(RFtimesInFrames * stimFrameDur, fullfile(folder, '_ss_rfDescr.timestamps.npy'))        
        writeNPY(pvals, fullfile(folder, '_ss_rf.pValues.npy'))
        writeNPY(fitPars, fullfile(folder, '_ss_rf.gaussFitPars.npy'))

        % plot results
        zTraces = (tr_up - nanmean(tr_up,1)) ./ nanstd(tr_up,0,1);
        [x0, y0] = meshgrid(linspace(stimPos(1), stimPos(2), 100), ...
            flip(linspace(stimPos(4), stimPos(3), 100)));
        for iUnit = 1:size(tr_up,2)
            rf = rFields(:,:,:,:,iUnit); % [rows x cols x t x ON/OFF]
            rf(:,:,:,2) = -rf(:,:,:,2);
            [mx,mxTime] = max(max(reshape(permute(abs(rf), [1 2 4 3]), ...
                [], size(rf,3)), [], 1));
            if ~isnan(fitPars(iUnit,1))
                fitRF = whiteNoise.D2GaussFunctionRot(fitPars(iUnit,:), cat(3, x0, y0));
            end
            figure('Position', [185 195 1360 780])
            for f = 1:2
                subplot(3,2,[0 2]+f)
                imagesc([stimPos(1)+squW/2 stimPos(2)-squW/2], ...
                    [stimPos(3)-squH/2 stimPos(4)+squH/2], ...
                    rf(:,:,mxTime,f),[-mx mx])
                if ~isnan(fitPars(iUnit,1))
                    hold on
                    contour(x0, y0, fitRF, [1 1].*fitPars(iUnit,1)/2, 'k', 'LineWidth', 1)
                end
                axis image
                set(gca, 'box', 'off', 'XTick', stimPos(1:2), ...
                    'YTick', [stimPos(4) 0 stimPos(3)], 'YDir', 'normal', ...
                    'YTickLabel', [stimPos(4) 0 stimPos(3)])
                xlim([stimPos(1) stimPos(2)])
                colormap(gca, colormaps{f})
                title(titles{f})
                colorbar
                if stimPos(1) * stimPos(2) < 0 % RF map crosses vertical midline
                    set(gca, 'XTick', [stimPos(1) 0 stimPos(2)])
                end
            end
            subplot(3,2,[5 6])
            hold on
            plot(t, zTraces(:,iUnit), 'k')
            plot(t_pred, pred(:,iUnit), 'r')
            legend('data', 'prediction')
            xlabel('Time (s)')
            ylabel('Delta-F (z-scored)')
            xlim(t([1 end]))
            sgtitle(sprintf('ROI %d (EV: %.1f%%, p: %.3f, lambda: %.0e', ...
                iUnit, maxEV(iUnit)*100, pvals(iUnit), lambdasStim(maxStimLam(iUnit))))

            saveas(gcf, fullfile(fResults, sprintf('unit%03d.png', iUnit)))
            close gcf
        end
    end
end