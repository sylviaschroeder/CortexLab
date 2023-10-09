%% Folders
folderBase = 'Z:\UCLData\2P_Task';
folderTools = 'C:\dev\toolboxes';
folderThisRepo = 'C:\dev\workspaces\CortexLab';
folderResults = 'C:\Users\Sylvia\OneDrive - University of Sussex\Lab\RESULTS\wheelTask\preprocessing';
folderSaveData = 'C:\Users\Sylvia\OneDrive - University of Sussex\Lab\DATA\NPY\task_2p';

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(fullfile(folderThisRepo))

%% Build db
make_db_wheelTask2023;

%% Parameters
timeGap = 600; % gap (in sec) introduced between experiments

%% Preprocess traces across experiments!
% parameters chosen on 19.12.2016
timeWindow = 180; % in sec (for neuron imaging)
prctileNeuron = 8;
prctileNpil = 40;
windowRed = 10; % in sec
for iSet = 1:length(db)
    fprintf('Dataset %d of %d\n', iSet, length(db))
    subject = db(iSet).mouse_name;
    date = db(iSet).date;
    nExps = length(db(iSet).expts);

    folder = fullfile(folderBase, subject, date);
    planes = dir(fullfile(folder, 'suite2p', 'plane*'));
    goodPlanes = ~endsWith({planes.name}, '0');
    nPlanes = length(planes);
    F_final = cell(nPlanes, nExps);
    frameTimes = cell(nPlanes, nExps);
    cellIDs = cell(nPlanes, 1);
    for iExp = 1:nExps
        expID = db(iSet).expts(iExp);
        file = dir(fullfile(folder, num2str(expID), sprintf('%s_*_Timeline.mat', date)));
        if length(file) > 1
            disp('WARNING! There are several Timeline files. Loading first.')
        end
        d = load(fullfile(folder, num2str(expID), file(1).name));
        timeLine = d.Timeline;
        clear d
        ft = ppbox.getFrameTimes2023(timeLine);
        for iP = 1:nPlanes
            frameTimes{iP,iExp} = ft(iP:nPlanes:end);
        end
    end
    samplingRate = 1/median(diff(frameTimes{1,1}));

    for iPlane = 1:nPlanes
        if ~goodPlanes(iPlane)
            continue
        end
        fprintf('  %s of %d\n', planes(iPlane).name, nPlanes)

        folderPlane = fullfile(folder, 'suite2p', planes(iPlane).name);
        isROI = readNPY(fullfile(folderPlane, 'iscell.npy'));
        isROI = logical(isROI(:,1));
        nROIs = sum(isROI);
        f = readNPY(fullfile(folderPlane, 'F.npy'))';
        n = readNPY(fullfile(folderPlane, 'Fneu.npy'))';
        f2 = readNPY(fullfile(folderPlane, 'F_chan2.npy'))';
        n2 = readNPY(fullfile(folderPlane, 'Fneu_chan2.npy'))';
        nframes = readNPY(fullfile(folderPlane, 'frames_per_folder.npy'));
        badframes = readNPY(fullfile(folderPlane, 'badframes.npy'));

        cellIDs{iPlane} = find(isROI);
        F = cell(1, nExps);
        N = cell(1, nExps);
        F_red = cell(1, nExps);
        N_red = cell(1, nExps);
        bad = cell(1, nExps);
        starts = cumsum([1; nframes(1:end-1)]);
        ends = cumsum(nframes);
        for iExp = 1:nExps
            F{iExp} = f(starts(iExp):ends(iExp), isROI);
            N{iExp} = n(starts(iExp):ends(iExp), isROI);
            F_red{iExp} = f2(starts(iExp):ends(iExp), isROI);
            N_red{iExp} = n2(starts(iExp):ends(iExp), isROI);
            bad{iExp} = badframes(starts(iExp):ends(iExp));
            frameTimes{iPlane,iExp} = frameTimes{iPlane,iExp}(1:nframes(iExp));
        end
        clear f n f2 n2

        % Correct for neuropil contamination
        % (1) high-pass filter traces of each exp. separately
        F_filt_all = [];
        F0_all = [];
        Npil_filt_all = [];
        Npil0_all = [];
        opt.minNp = 5;
        opt.maxNp = 95;
        opt.constrainedFit = 1;
        opt.verbose = 0;
        for iExp = 1:nExps
            [F_filt, F0_tmp] = preproc.removeSlowDrift(F{iExp}, samplingRate, ...
                timeWindow, prctileNeuron);
            F_filt_all = [F_filt_all; F_filt];
            F0_all = [F0_all; F0_tmp];
            [Npil_filt, Npil0] = preproc.removeSlowDrift(N{iExp}, ...
                samplingRate, timeWindow, prctileNpil);
            Npil_filt_all = [Npil_filt_all; Npil_filt];
            Npil0_all = [Npil0_all; Npil0];
        end
        % (2) find common neuropil correction parameters across experiments
        % and perform correction on each experiment
        [~,t] = min((F0_all-Npil0_all) ./ Npil0_all, [], 1);
        ind = sub2ind(size(F0_all), t, 1:size(F0_all,2));
        [~, parameters] = preproc.estimateNeuropil(bsxfun(@plus, F_filt_all, F0_all(ind))', ...
            bsxfun(@plus, Npil_filt_all, Npil0_all(ind))', opt);
        F_woNpil = cell(1,nExps);
        for iExp = 1:nExps
            F_woNpil{iExp} = F{iExp} - bsxfun(@times, parameters.corrFactor(:,2)', ...
                N{iExp});
        end
        % (3) remove slow drift from neuropil corrected traces in each
        % experiment
        F0_all = [];
        F0 = cell(1,nExps);
        deltaF = cell(1,nExps);
        for iExp = 1:nExps
            [deltaF{iExp}, F0{iExp}] = preproc.removeSlowDrift( ...
                F_woNpil{iExp}, samplingRate, timeWindow, prctileNeuron);
            F0_all = [F0_all; F0{iExp}];
        end
        % (4) calculate delta-F-over-F with same denominator across
        % experiments
        F_delta = cell(1,nExps);
        for iExp = 1:nExps
            F_delta{iExp} = bsxfun(@rdivide, deltaF{iExp}, ...
                max(1, mean(F0_all, 1)));
        end

        % (5) remove slow drift of red traces for each experiment
        R_all = [];
        F_all = [];
        R0 = cell(1,nExps);
        R_final = cell(1,nExps);
        nRF = round(windowRed * samplingRate);
        for iExp = 1:nExps
            % remove exponential decay
            R = NaN(size(F_red{iExp}));
            ind = round(150 * samplingRate);
            for iCell = 1:nROIs
                y = double(F_red{iExp}(1:ind,iCell));
                f = fit((1:ind)', y, ...
                    @(a,b,c,d,e,x) a+b.*exp(-x./c)+d.*exp(-x./e), ...
                    'Lower', [0 0 0 0 0], ...
                    'Upper', [max(y) max(y) 500 max(y) 500], ...
                    'StartPoint', [min(y) mean(y) 50 mean(y) 5]);
                R(:,iCell) = F_red{iExp}(:,iCell) - f(1:size(R,1)) + f.a;
            end
            F_red{iExp} = R;

            [R_noDrift, R0{iExp}] = preproc.removeSlowDrift(F_red{iExp}, ...
                samplingRate, timeWindow, prctileNeuron);
            R_final{iExp} = medfilt1(R_noDrift, nRF, [], 1, 'includenan', 'truncate');
            R_all = [R_all; R_final{iExp}];
            F_all = [F_all; F_delta{iExp}];
        end
        % (6) fit red traces to green traces across all experiments at once
        slopes = zeros(1, nROIs);
        for c = 1:nROIs
            if all(any(isnan([R_all(:,c) F_all(:,c)]),2))
                continue
            end
            mdl = fitlm(R_all(:,c), F_all(:,c), 'RobustOpts', 'on');
            slopes(c) = mdl.Coefficients.Estimate(2);
        end
        % (7) subtract weighted red traces from green traces
        for iExp = 1:nExps
            F_final{iPlane,iExp} = F_delta{iExp} - slopes .* R_final{iExp};
            F_final{iPlane,iExp}(bad{iExp},:) = NaN;
        end

        % Plot various traces
        numPl = 4;
        fPlots = fullfile(folderResults, sprintf('%s_%s', subject, date), ...
            planes(iPlane).name);
        if ~isfolder(fPlots)
            mkdir(fPlots);
        end
        for c = 1:nROIs
            ax = zeros(1,4);

            t0 = 0;
            figure('Position', [10 50 1900 930])
            for iExp = 1:nExps
                t = frameTimes{iPlane,iExp} + t0;

                subplot(numPl,1,1)
                hold on
                plot(t, F{iExp}(:,c), 'k')
                plot(t, F_woNpil{iExp}(:,c), 'Color', [0 .5 .5])
                plot(t, F0{iExp}(:,c), 'Color', [0 1 1])
                legend('raw','npil. corr.','F0')
                axis tight
                ax(1) = gca;

                subplot(numPl,1,2)
                hold on
                plot(t, F_delta{iExp}(:,c), 'Color', [0 0 .5])
                plot(t, F_final{iPlane,iExp}(:,c), 'Color', [0 .5 0])
                legend('\DeltaF/F','F final')
                axis tight
                ax(2) = gca;

                subplot(numPl,1,3)
                hold on
                plot(t, R_final{iExp}(:,c), 'Color', [.5 0 0])
                legend('R final')
                axis tight
                ax(3) = gca;

                subplot(numPl,1,4)
                hold on
                plot(t, F_red{iExp}(:,c), 'Color', [.5 0 .5])
                plot(t, R0{iExp}(:,c), 'Color', [1 0 1])
                legend('raw R','R0')
                axis tight
                ax(4) = gca;

                t0 = t(end) + 200;
            end
            xlabel('Time (sec)')

            linkaxes(ax, 'x')

            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(fullfile(fPlots, sprintf('preprocessing_%03d.jpg', c)), ...
                '-djpeg','-r0')

            close gcf
        end
    end

    % Save data
    f = fullfile(folderSaveData, subject, date);
    if ~isfolder(f)
        mkdir(f);
    end

    goodPlanes = find(goodPlanes);
    [sz,~] = cellfun(@size, F_final(goodPlanes,:));
    sz = min(sz, [], 1);
    time = [];
    dff = [];
    t0 = 0;
    for iExp = 1:nExps
        time = [time; frameTimes{goodPlanes(1),iExp}(1:sz(iExp)) + t0];
        dff_exp = [];
        for iPlane = goodPlanes
            dff_exp = [dff_exp, F_final{iPlane,iExp}(1:sz(iExp),:)];
        end
        dff = [dff; dff_exp];
        t0 = time(end) + timeGap;
    end
    delays = [];
    pl = [];
    for iPlane = 1:length(goodPlanes)
        delays = [delays; median(frameTimes{goodPlanes(iPlane),1}(1:sz(iExp)) - ...
            frameTimes{goodPlanes(1),1}(1:sz(iExp)))];
        pl = [pl; ones(size(F_final{goodPlanes(iPlane),1}, 2), 1) .* iPlane];
    end
    
    writeNPY(time, fullfile(f, '_ss_2pCalcium.timestamps.npy'))
    writeNPY(dff, fullfile(f, '_ss_2pCalcium.dff.npy'))
    writeNPY(delays, fullfile(f, '_ss_2pPlanes.delay.npy'))
    writeNPY(pl, fullfile(f, '_ss_2pRois._ss_2pPlanes.npy'))
    writeNPY(cat(1, cellIDs{:}), fullfile(f, '_ss_2pRois.ids.npy'))
end

%% Classify cells into GAD-positive and -negative neurons
% clear opt
% opt.classThresholds = [55 85];
% opt.bloodThreshold = 20;
% opt.width = 20;
% for iSet = 1:length(db)
%     folder = fullfile(ops.InfoStorage, db(iSet).mouse_name, db(iSet).date);
%     files = [db(iSet).date '_%d_' db(iSet).mouse_name '_2P_plane%03d_ROI.mat'];
%     planes = db(iSet).planesToProcess;
%     expts = 1:length(db(iSet).expts);
%     for iPlane = planes
%         iExp = expts(1);
%         file = fullfile(folder, num2str(iExp), sprintf(files, iExp, iPlane));
%         data = load(file);
%         meta = data.meta;
%         data = load(fullfile(meta.folderDat, meta.filenameDat));
%         stat = data.dat.stat([data.dat.stat.iscell] == 1);
%         opt.zoom = meta.zoomFactor;
%         opt.scope = meta.microID;
%         opt.xrange = data.dat.ops.xrange;
%         opt.yrange = data.dat.ops.yrange;
%         [classes, opt2] = classifyCells(meta.targetFrame, ...
%             meta.targetFrameROI, opt, stat);
%         meta.ROI.isGad = classes;
%         meta.opsGAD = opt2;
%         save(file, 'meta');
%         for iExp = expts(2:end)
%             file = fullfile(folder, num2str(iExp), sprintf(files, iExp, iPlane));
%             data = load(file);
%             meta = data.meta;
%             meta.ROI.isGad = classes;
%             meta.opsGAD = opt2;
%             save(file, 'meta');
%         end
%     end
% end