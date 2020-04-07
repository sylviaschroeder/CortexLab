%% Parameters
% label = 'neurons';
label = 'boutons';

smoothStd = 1; %in sec

if strcmp(label, 'neurons')
    nonvis = 'pupil';
    stim = 'grayScreen';
else
    nonvis = 'running';
    stim = 'dark';
end

%% Folders
% file location depending on PCs
folderPC = 'C:\STORAGE\OneDrive - University College London';
% folderPC = 'C:\Users\skgtchr\OneDrive - University College London'; % ZEN notebook
% folderPC = 'C:\Storage\OneDrive - University College London'; % Lenovo laptop
% data folders
folderROIData = fullfile(folderPC, 'Lab\DATA\InfoStructs');
if strcmp(label, 'neurons')
    folderResults = fullfile(folderPC, 'Lab\RESULTS\nonvisualEffects\modelGratingResp');
    db_driftingGratings_blank
else
    folderResults = fullfile(folderPC, 'Lab\RESULTS\boutons\nonvisualEffects');
    db_boutons_driftingGratings_blanks
end
folderPlots = fullfile(folderResults, 'plots_pupil', 'stimulusTraces');

%% Load data
data = load(fullfile(folderResults, 'kernelFit', 'results.mat'));
kernels = data.results;
data = load(fullfile(folderResults, 'pupil', 'tuning_prefDirSigmaDIFixed_lineFit.mat'));
tuning = data.tuning;
degrees = data.x;
data = load(fullfile(folderResults, 'pupil', 'nullTuning_prefDirSigmaDIFixed.mat'));
null = data.null;
data = load(fullfile(folderResults, 'pupil', 'corrsDuringGratingsAndGrayScreen_sigma1.00.mat'));
corrs = data.corrs;
if strcmp(label, 'neurons')
    corrsSpont = corrs;
else
    data = load(fullfile(folderResults, nonvis, 'corrsDuringGratingsAndGrayScreen_sigma1.00.mat'));
    corrsSpont = data.corrs;
end

%% Plot traces for each neuron
colors = [0 0 0; 1 0 0];
offsets = [5 5]; % in s, time before and after stimulus to be plotted

for iExp = 1:length(tuning)
    folder = fullfile(folderROIData, tuning(iExp).subject, ...
        tuning(iExp).date);
    file = [tuning(iExp).date '_%d_' tuning(iExp).subject '_2P_plane%03d_ROI.mat'];
    for iPlane = 1:length(tuning(iExp).plane)
        % load raw plane data for gratings
        data = load(fullfile(folder, num2str(tuning(iExp).exp), ...
            sprintf(file,tuning(iExp).exp, tuning(iExp).planes(iPlane))));
        meta = data.meta;
        meta.folderTL = strrep(meta.folderTL, 'ioo.ucl.ac.uk', ...
            'cortexlab.net');
        traceTimes = ppbox.getFrameTimes(meta);
        traceBin = median(diff(traceTimes));
        offSamples = round(offsets ./ traceBin);
        stdSamples = round(smoothStd / traceBin);
        convWindow = normpdf(-4*stdSamples:4*stdSamples, 0, stdSamples);
        
        if strcmp(label, 'neurons')
            % load raw plane data for gray screens
            if ~isempty(db(iExp).expGrayScreen)
                data = load(fullfile(folder, num2str(db(iExp).expGrayScreen), ...
                    sprintf(file,db(iExp).expGrayScreen, tuning(iExp).planes(iPlane))));
                metaSpont = data.meta;
                metaSpont.folderTL = strrep(metaSpont.folderTL, 'ioo.ucl.ac.uk', ...
                    'cortexlab.net');
                traceTimesSpont = ppbox.getFrameTimes(metaSpont);
                
                if iPlane == 1
                    % get pupil data
                    [nonvisSpont, nonvisTimeSpont] = nonVis.loadPupilData(metaSpont);
                    if ~isempty(nonvisSpont)
                        nonvisTimeSpont(length(nonvisSpont.x)+1:end) = [];
                        nonvisSpont = nonVis.getPupilDiam(nonvisSpont);
                    end
                end
            else
                metaSpont = [];
            end
        else
            % load raw plane data for darkness
            if ~isempty(db(iExp).expDark)
                data = load(fullfile(folder, num2str(db(iExp).expDark), ...
                    sprintf(file,db(iExp).expDark, tuning(iExp).planes(iPlane))));
                metaSpont = data.meta;
                metaSpont.folderTL = strrep(metaSpont.folderTL, 'ioo.ucl.ac.uk', ...
                    'cortexlab.net');
                traceTimesSpont = ppbox.getFrameTimes(metaSpont);
                
                if iPlane == 1
                    % get running data
                    ballData = nonVis.getRunningSpeed(metaSpont);
                    if ~isempty(ballData)
                        nonvisSpont = ballData.total / median(diff(ballData.t)) / 53;
                        nonvisTimeSpont = ballData.t;
                    end
                end
            else
                metaSpont = [];
            end
        end
        
        if iPlane == 1
            % get stimulus information
            [stimTimes, stimSeq, stimMatrix, ~, samplingRate] = ...
                ssLocal.getStimulusResponseInfo(meta);
            stimDuration = mean(stimTimes.offset - stimTimes.onset);
            [orientations, blankStim] = gratings.getOrientations(stimSeq);
            plotRows = floor(sqrt(size(orientations,1)));
            plotCols = ceil(size(orientations,1) / plotRows);
            
            % get pupil data
            [pupil, pupilTime] = nonVis.loadPupilData(meta);
            if ~isempty(pupil)
                pupilTime(length(pupil.x)+1:end) = [];
                pupil = nonVis.getPupilDiam(pupil);
            end
            
            % get condition of trials according to pupil size (large or
            % small)
            validPlanes = find(~cellfun(@isempty, {tuning(iExp).plane.cellIDs}));
            validROI = [];
            for p = 1:length(validPlanes)
                validROI = find(~isnan(tuning(iExp).plane(validPlanes(p)).isTuned), 1);
                if ~isempty(validROI)
                    break
                end
            end
            if isempty(validROI)
                continue
            end
            conditions = NaN(size(tuning(iExp).plane(validPlanes(p)) ...
                .cond(1).cell(validROI).responses));
            for c = 1:2
                conditions(~isnan(tuning(iExp).plane(validPlanes(p)).cond(c) ...
                    .cell(validROI).responses)) = c;
            end
            conditions(end+1,:) = NaN;
            seq = reshape(stimSeq.seq, size(stimMatrix,1), []);
            seq = reshape(seq + (0:size(seq,2)-1).*size(seq,1), [], 1);
            condsSeq = reshape(conditions,[],1);
            condsSeq = condsSeq(seq);
        end
        
        if ~isempty(pupil)
            ind = isnan(pupil);
            indInterp = hist(pupilTime(ind), traceTimes) > 0;
            pupilFilt = interp1(pupilTime(~ind), pupil(~ind), traceTimes, 'pchip')';
            pupilFilt = conv([ones(ceil((length(convWindow)-1)/2),1) .* ...
                mean(pupilFilt(1:round(length(convWindow)/2))); ...
                pupilFilt; ones(floor((length(convWindow)-1)/2),1) .* ...
                mean(pupilFilt(end-round(length(convWindow)/2):end))], ...
                convWindow, 'valid');
            pupilFilt(indInterp) = NaN;
        end
        if ~isempty(nonvisSpont)
            ind = isnan(nonvisSpont);
            indInterp = hist(nonvisTimeSpont(ind), traceTimesSpont) > 0;
            nonvisSpont = interp1(nonvisTimeSpont(~ind), nonvisSpont(~ind), ...
                traceTimesSpont, 'pchip')';
            nonvisSpont = conv([ones(ceil((length(convWindow)-1)/2),1) .* ...
                mean(nonvisSpont(1:round(length(convWindow)/2))); ...
                nonvisSpont; ones(floor((length(convWindow)-1)/2),1) .* ...
                mean(nonvisSpont(end-round(length(convWindow)/2):end))], ...
                convWindow, 'valid');
            nonvisSpont(indInterp) = NaN;
        end
        
        responses = ssLocal.getTracesPerStimulus(meta.F_final, ...
            stimMatrix, offSamples);
        respTime = ((1:size(responses,4))-1-offSamples(1)) .* traceBin;
        
        plotF = fullfile(folderPlots, sprintf('%s_%s_%d', ...
            tuning(iExp).subject, tuning(iExp).date, tuning(iExp).exp), ...
            num2str(tuning(iExp).planes(iPlane)));
        if ~isfolder(plotF)
            mkdir(plotF)
        end
        
        for iCell = 1:length(tuning(iExp).plane(iPlane).cellIDs)
            ID = tuning(iExp).plane(iPlane).cellIDs(iCell);
            krnID = kernels(iExp).plane(iPlane).cellIDs == ID;
            if isempty(kernels(iExp).plane(iPlane).kernelFit(krnID).kernel)
                continue
            end
            
            minis = [NaN NaN];
            maxis = [NaN NaN];
            stimMeans = [NaN NaN];
            nullMinis = [];
            nullMaxis = [];
            nullStimMeans = [];
            for c = 1:2
                pars = tuning(iExp).plane(iPlane).cond(c).cell(iCell).parameters;
                curve = tuning(iExp).plane(iPlane).cond(c).cell(iCell).curve;
                if length(pars) == 1 % not tuned
                    stimMeans(c) = pars;
                else
                    stimMeans(c) = mean(curve);
                    oris = mod(pars(1) + [0 90 180], 360);
                    sr = gratings.orituneWrappedConditions(pars, oris);
                    maxis(c) = sr(1);
                    if sr(1)-sr(2)>0
                        [minis(c), ind] = min(sr(2:3));
                    else % suppressed neurons
                        [minis(c), ind] = max(sr(2:3));
                    end
                    ind = ind+1;
                end
                
                pars = null(iExp).plane(iPlane).cond(c).cell(iCell).parameters;
                if size(pars,2) == 1 % not tuned
                    nullStimMeans = [nullStimMeans, pars];
                else
                    sr = NaN(size(pars,1), 3);
                    curves = NaN(size(pars,1), length(degrees));
                    for p = 1:size(pars,1)
                        oris = mod(pars(p,1) + [0 90 180], 360);
                        sr(p,:) = gratings.orituneWrappedConditions(pars(p,:), oris);
                        curves(p,:) = gratings.orituneWrappedConditions(pars(p,:), degrees);
                    end
                    nullStimMeans = [nullStimMeans, mean(curves,2)];
                    nullMaxis = [nullMaxis, sr(:,1)];
                    nullMinis = [nullMinis, sr(:,ind)];
                end
            end
            mods = abs(maxis - minis);
            nullMods = abs(nullMaxis - nullMinis);
            modFun = @(a,b) (b-a)./(abs(a)+abs(b));
            DImaxis = modFun(maxis(1), maxis(2));
            DImeans = modFun(stimMeans(1), stimMeans(2));
            DImods = modFun(mods(1), mods(2));
            nullDImeans = prctile(modFun(nullStimMeans(:,1), nullStimMeans(:,2)), [2.5 97.5]);
            if ~isempty(nullMaxis)
                nullDImaxis = prctile(modFun(nullMaxis(:,1), nullMaxis(:,2)), [2.5 97.5]);
                nullDImods = prctile(modFun(nullMods(:,1), nullMods(:,2)), [2.5 97.5]);
            else
                nullDImaxis = [NaN NaN];
                nullDImods = [NaN NaN];
            end
            
            % Plot continuous pupil and calcium traces, and fitted kernel
            figure('Position', [1921 1 1920 1123])
            annotation('textbox', [0 .95 1 .03], 'String', ...
                sprintf('ROI %d: DI mean %.2f [%.2f %.2f], DI pref %.2f [%.2f %.2f], DI depth %.2f [%.2f %.2f]', ...
                ID, DImeans, nullDImeans, DImaxis, nullDImaxis, DImods, nullDImods), ...
                'HorizontalAlignment', 'Center', 'FontSize', 14, ...
                'FontWeight', 'bold', 'LineStyle', 'none')
            ax = [0 0];
            ps = medfilt1(pupil, 21);
            tr = medfilt1(meta.F_final(:,ID), 11);
            pred = kernels(iExp).plane(iPlane).kernelFit(krnID).prediction;
            minis = [min(ps), min([tr;pred])];
            maxis = [max(ps), max([tr;pred])];
            crr = corrs(iExp).plane(iPlane).gratings.rhos(ID);
            nl = prctile(corrs(iExp).plane(iPlane).gratings.nullRhos( ...
                ID,:), [2.5 97.5]);
            t = [stimTimes.onset, stimTimes.offset];
            h = zeros(1,4);
            for p = 1:2
                subplot(3,3,(p-1)*3+(1:3))
                hold on
                for c = 1:2
                    tt = t(condsSeq == c,:);
                    tt = [tt, flip(tt,2)]';
                    hh = fill(tt, [[1 1]'.*minis(p); [1 1]'.*maxis(p)], 'k', 'FaceColor', ...
                        colors(c,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
                    h((p-1)*2+c) = hh(1);
                end
            end
            % for gratings
            subplot(3,3,1:3)
            ax(1) = gca;
            plot(pupilTime, ps, 'k', 'LineWidth', 1)
            ylim([minis(1) maxis(1)])
            legend(h(1:2), {'small pupil','large pupil'})
            ylabel('Pupil size')
            title(sprintf('During gratings, correlation with pupil %.3f (null: [%.3f %.3f])', ...
                crr, nl(1), nl(2)))
            subplot(3,3,4:6)
            h = [0 0];
            ax(2) = gca;
            mn = min(length(traceTimes), length(pred));
            h(1) = plot(traceTimes, tr, 'k', 'LineWidth', 1);
            h(2) = plot(traceTimes(1:mn), pred(1:mn), 'LineWidth', 1);
            ylim([minis(2) maxis(2)])
            legend(h, {'data', 'prediction'});
            linkaxes(ax, 'x')
            xlim(traceTimes([1 end]))
            xlabel('Time (in s)')
            ylabel('\DeltaF/F')
            subplot(3,3,7)
            krnl = kernels(iExp).plane(iPlane).kernelFit(krnID).kernel;
            hold on
            plot([1 1].*stimDuration, [min(krnl) max(krnl)], 'k:')
            plot(kernels(iExp).plane(iPlane).kernelTime, krnl, 'k', ...
                'LineWidth', 2)
            axis tight
            xlabel('Time from stimulus onset (s)')
            ylabel('Norm. resp.')
            title('Fitted kernel')
            
            % Plot continuous traces, filtered in same way as for
            % determining correlation values
            figure('Position', [1921 1 1920 1123])
            annotation('textbox', [0 .95 1 .03], 'String', ...
                sprintf('ROI %d: DI mean %.2f [%.2f %.2f], DI pref %.2f [%.2f %.2f], DI depth %.2f [%.2f %.2f]', ...
                ID, DImeans, nullDImeans, DImaxis, nullDImaxis, DImods, nullDImods), ...
                'HorizontalAlignment', 'Center', 'FontSize', 14, ...
                'FontWeight', 'bold', 'LineStyle', 'none')
            % for gratings
            crr = corrs(iExp).plane(iPlane).gratings.rhos(ID);
            nl = prctile(corrs(iExp).plane(iPlane).gratings.nullRhos( ...
                ID,:), [2.5 97.5]);
            subplot(4,1,1)
            plot(traceTimes, pupilFilt, 'k', 'LineWidth', 1)
            axis tight
            set(gca,'box','off')
            ylabel('Pupil size')
            title(sprintf('During gratings, correlation with pupil %.3f (null: [%.3f %.3f])', ...
                crr, nl(1), nl(2)))
            subplot(4,1,2)
            tr = meta.F_final(:,ID);
            tr = conv([ones(ceil((length(convWindow)-1)/2),1) .* ...
                mean(tr(1:round(length(convWindow)/2))); ...
                tr; ...
                ones(floor((length(convWindow)-1)/2),1) .* ...
                mean(tr(end-round(length(convWindow)/2)))], ...
                convWindow, 'valid');
            plot(traceTimes, tr, 'k', 'LineWidth', 1)
            axis tight
            set(gca,'box','off')
            xlabel('Time (s)')
            ylabel('\DeltaF/F')            
            % for spontaneous activity
            if ~isempty(metaSpont)
                crr = corrsSpont(iExp).plane(iPlane).(stim).rhos(ID);
                nl = prctile(corrsSpont(iExp).plane(iPlane).(stim).nullRhos( ...
                    ID,:), [2.5 97.5]);
                
                subplot(4,1,3)
                ax(1) = gca;
                plot(traceTimesSpont, nonvisSpont, 'k', 'LineWidth', 1)
                axis tight
                set(gca,'box','off')
                ylabel(nonvis)
                title(sprintf('During %s, correlation with %s %.3f (null: [%.3f %.3f])', ...
                    stim, nonvis, crr, nl(1), nl(2)))
                subplot(4,1,4)
                ax(2) = gca;
                trSpont = metaSpont.F_final(:,ID);
                trSpont = conv([ones(ceil((length(convWindow)-1)/2),1) .* ...
                    mean(trSpont(1:round(length(convWindow)/2))); ...
                    trSpont; ...
                    ones(floor((length(convWindow)-1)/2),1) .* ...
                    mean(trSpont(end-round(length(convWindow)/2)))], ...
                    convWindow, 'valid');
                plot(traceTimesSpont, trSpont, 'k', 'LineWidth', 1);
                linkaxes(ax, 'x')
                axis tight
                set(gca,'box','off')
                xlim(traceTimesSpont([1 end]))
                xlabel('Time (in s)')
                ylabel('\DeltaF/F')
            end
            
            figure
            plot(pupilFilt, tr, 'k.')
            xlabel('Pupil')
            ylabel('\DeltaF/F')
            
%             savefig(fullfile(plotF, sprintf('ROI%03d_contTrace.fig', ...
%                 ID)))
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(fullfile(plotF, sprintf('ROI%03d_contTrace.jpg', ...
                ID)),'-djpeg','-r0')
            close gcf
            
            % Plot responses aligned to stimulus onset
            figure('Position', [1921 1 1920 1123])
            r = permute(responses(ID, :, :, :), [4 3 2 1]); % [time x trials x ori]
            yLimits = [min(r(:)) max(r(:))];
            h = [0 0];
            for ori = 1:size(orientations,1)
                subplot(plotRows, plotCols, ori)
                hold on
                plot([0 0], yLimits, 'k')
                plot([1 1].*stimDuration, yLimits, 'k')
                for c = 1:2
                    plot(respTime, r(:,conditions(ori,:)==c, ori), 'Color', ...
                        colors(c,:).*0.3 + [1 1 1].*0.7)
                end
                for c = 1:2
                    h(c) = plot(respTime, nanmean(r(:,conditions(ori,:)==c, ori),2), ...
                        'Color', colors(c,:), 'LineWidth', 2);
                end
                xlim(respTime([1 end]))
                ylim(yLimits)
                title(sprintf('Direction: %d deg', orientations(ori,1)))
                if ori == (plotCols-1)*plotRows + 1
                    xlabel('Time from stimulus onset (s)')
                    ylabel('\DeltaF/F')
                    legend(h, {'small pupil','large pupil'})
                end
            end
            annotation('textbox', [0 .95 1 .03], 'String', ...
                sprintf('ROI %d: DI mean %.2f [%.2f %.2f], DI pref %.2f [%.2f %.2f], DI depth %.2f [%.2f %.2f]', ...
                ID, DImeans, nullDImeans, DImaxis, nullDImaxis, DImods, nullDImods), ...
                'HorizontalAlignment', 'Center', 'FontSize', 14, ...
                'FontWeight', 'bold', 'LineStyle', 'none')
%             savefig(fullfile(plotF, sprintf('ROI%03d_stimTraces.fig', ...
%                 ID)))
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(fullfile(plotF, sprintf('ROI%03d_stimTraces.jpg', ...
                ID)),'-djpeg','-r0')
            close gcf
        end
    end
end