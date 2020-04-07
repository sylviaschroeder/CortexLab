%% Define data
label = 'boutons';

%% Folders
% file location depending on PCs
folderPC = 'C:\STORAGE\OneDrive - University College London';
folderROIData = fullfile(folderPC, 'Lab\DATA\InfoStructs');

if strcmp(label, 'boutons')
    make_db_boutons;
end

% need to determine N_0 for each experiment
% collect F_0 and N_0 across all experiments, then average across
% experiments

%% Preprocess traces across experiments!

corrections = db;
[corrections.subject] = deal(db.mouse_name);
corrections = rmfield(corrections, 'mouse_name');
corrections = orderfields(corrections, [6 1:5]);

doCorrect = @(a,b,F)(a.*F + b);

for iSet = 1:length(db)
    fprintf('Dataset %d of %d\n', iSet, length(db))
    if iSet < 7
        expIndices = db(iSet).expts;
    else
        expIndices = 1:length(db(iSet).expts);
    end
    corrections(iSet).expts = expIndices;
    for iPlane = 1:length(db(iSet).planesToProcess)
        fprintf('  plane %d of %d\n', iPlane, length(db(iSet).planesToProcess))
        clear meta
        for iExp = 1:length(db(iSet).expts)
            % load meta (with green and red traces)
            data = load(fullfile(folderROIData, db(iSet).mouse_name, ...
                db(iSet).date, num2str(expIndices(iExp)), sprintf( ...
                '%s_%d_%s_2P_plane%03d_ROI.mat', db(iSet).date, ...
                expIndices(iExp), db(iSet).mouse_name, ...
                db(iSet).planesToProcess(iPlane))));
            meta(iExp) = data.meta;
        end
        
        frameTimes = ppbox.getFrameTimes(meta(1));
        samplingRate = 1/median(diff(frameTimes));
        
        timeWindow = meta(1).opsNpil.driftWindow;
        prctileNeuron = meta(1).opsNpil.prctileNeuron;
        prctileNpil = meta(1).opsNpil.prctileNpil;
        windowRed = NaN;
        if isfield(meta, 'redFiltWindow')
            windowRed = meta(1).redFiltWindow;
        end
            
        % Determine parameters needed to calculate correction
        F = cell(1, length(db(iSet).expts));
        F0_corr = cell(1, length(db(iSet).expts));
        F_filt = cell(1, length(db(iSet).expts));
        F0 = cell(1, length(db(iSet).expts));
        Npil_filt = cell(1, length(db(iSet).expts));
        Npil0 = cell(1, length(db(iSet).expts));
        NpilSlopes = cell(1, length(db(iSet).expts));
        NpilIntercepts = cell(1, length(db(iSet).expts));
        opt.minNp = meta(1).opsNpil.minNp;
        opt.maxNp = meta(1).opsNpil.maxNp;
        opt.constrainedFit = meta(1).opsNpil.constrainedFit;
        opt.verbose = 0;
       
        for iExp = 1:length(meta)
            fprintf('    experiment %d of %d\n', iExp, length(db(iSet).expts))
            [f, f0] = preproc.removeSlowDrift(meta(iExp).F, samplingRate, ...
                timeWindow, prctileNeuron);
            F_filt{iExp} = f;
            F0{iExp} = f0;
            [n, n0] = preproc.removeSlowDrift(meta(iExp).Npil, ...
                samplingRate, timeWindow, prctileNpil);
            Npil_filt{iExp} = n;
            Npil0{iExp} = n0;
            NpilSlopes{iExp} = meta(iExp).NpilSlopes;
            NpilIntercepts{iExp} = meta(iExp).NpilIntercepts;
            
            F0_corr{iExp} = meta(iExp).F0;
            F{iExp} = meta(iExp).F;
        end
        F_filt = cat(1, F_filt{:});
        F0_all = cat(1, F0{:});
        Npil_filt = cat(1, Npil_filt{:});
        Npil0_all = cat(1, Npil0{:});
        F0_corr_all = cat(1, F0_corr{:});
        F = cat(1, F{:});
        
        % for datasets where magnitude (slope) for neuropil correction
        % differs across experiments, find the magnitude across all
        % experiments
        if ~all(all(cat(1,NpilSlopes{:}) == NpilSlopes{1}, 1)) && ...
                ~all(all(cat(1,NpilIntercepts{:}) == NpilIntercepts{1}, 1))
            separate = true;
            % (2) find common neuropil correction parameters across experiments
            % and perform correction on each experiment
            [~,t] = min((F0_all-Npil0_all) ./ Npil0_all, [], 1);
            ind = sub2ind(size(F0_all), t, 1:size(F0_all,2));
            goodROIs = ~all(isnan(F_filt),1) & ~all(isnan(Npil_filt),1);
            [~, parameters] = preproc.estimateNeuropil(bsxfun(@plus, ...
                F_filt(:,goodROIs), F0_all(ind(goodROIs)))', ...
                bsxfun(@plus, Npil_filt(:,goodROIs), ...
                Npil0_all(ind(goodROIs)))', opt);
            newNpilSlopes = NaN(1, size(F,2));
            newNpilSlopes(goodROIs) = parameters.corrFactor(:,2)';
            indChange = true(1, length(newNpilSlopes));
        else
            separate = false;
            newNpilSlopes = NpilSlopes{1};
            indChange = false(1, length(newNpilSlopes));
        end
        
        % the target slope is between [0 2], adjust new slopes accordingly
        indChange = indChange | newNpilSlopes < 0 | newNpilSlopes > 2;
        newNpilSlopes(newNpilSlopes < 0) = 0;
        newNpilSlopes(newNpilSlopes > 2) = 2;
        % make sure that new slope <= 0.9 * prctile(F,prctileNeuron)/prctile(Npil,prctileNpil) so
        % delta-F = F - slope*Npil >= F - 0.9*prctile(F,prctileNeuron)
        ratios = prctile(F, prctileNeuron) ./ prctile(Npil0_all, prctileNpil);
        ind = newNpilSlopes > 0.9 .* ratios;
        newNpilSlopes(ind) = 0.9 .* ratios(ind);
        indChange = indChange | ind;
        indChange = find(indChange);
        
        %  perform neuropil correction, remove slow drift from neuropil 
        % corrected traces, remove slow drift of red traces (only for ROIs
        % that change)
        deltaF2 = cell(1, length(db(iSet).expts));
        F2_0 = cell(1, length(db(iSet).expts));
        R = cell(1, length(db(iSet).expts));
        for iExp = 1:length(meta)
            F2_woNpil = meta(iExp).F(:,indChange) - ...
                bsxfun(@times, newNpilSlopes(indChange), ...
                meta(iExp).Npil(:,indChange));
            [deltaF2{iExp}, F2_0{iExp}] = preproc.removeSlowDrift(F2_woNpil, ...
                samplingRate, timeWindow, prctileNeuron);
            
            if ~isfield(meta, 'R')
                R{iExp} = zeros(size(meta(iExp).F,1), length(indChange));
                continue
            end
            if strcmp(label, 'boutons')
                r = NaN(size(meta(iExp).R,1), length(indChange));
                ind = round(150 * samplingRate);
                for iCell = 1:length(indChange)
                    y = meta(iExp).R(1:ind,indChange(iCell));
                    f = fit((1:ind)', y, ...
                        @(a,b,c,d,e,x) a+b.*exp(-x./c)+d.*exp(-x./e), ...
                        'Lower', [0 0 0 0 0], ...
                        'Upper', [max(y) max(y) 500 max(y) 500], ...
                        'StartPoint', [min(y) mean(y) 50 mean(y) 5]);
                    r(:,iCell) = meta(iExp).R(:,indChange(iCell)) - f(1:size(R,1)) + f.a;
                end
            else
                r = meta(iExp).R(:,indChange);
            end
            % remove slow drift and filter red traces
            R_noDrift = preproc.removeSlowDrift(r, ...
                samplingRate, timeWindow, prctileNeuron);
            nRF = round(windowRed * samplingRate);
            R{iExp} = medfilt1(R_noDrift, nRF, [], 1, 'includenan', 'truncate');
        end
        R_all = cat(1, R{:});
        
        % calculate delta-F-over-F with same denominator across experiments
        % (only for ROIs that change)
        F2_delta = cellfun(@rdivide, deltaF2, ...
            repmat({max(1, nanmean(cat(1, F2_0{:})))}, 1, length(deltaF2)), ...
            'UniformOutput', false);
        F2_delta_all = cat(1, F2_delta{:});
        
        % fit red traces to green traces across all experiments (only for 
        % ROIs that change)
        slopes = zeros(1, length(indChange));
        if isfield(meta, 'R')
            for c = 1:length(indChange)
                if all(any(isnan([R_all(:,c) F2_delta_all(:,c)]),2))
                    continue
                end
                mdl = fitlm(R_all(:,c), F2_delta_all(:,c), 'RobustOpts', 'on');
                slopes(c) = mdl.Coefficients.Estimate(2);
            end
        end
        
        % subtract weighted red traces from green traces (for all ROIs!)
        F2_final = cell(1, max(expIndices));
        for iExp = 1:length(meta)
            F2_final{expIndices(iExp)} = meta(iExp).F_final;
            F2_final{expIndices(iExp)}(:,indChange) = ...
                F2_delta{iExp} - slopes .* R{iExp};
        end
        
        a = cell(1, max(expIndices));
        b = cell(1, max(expIndices));
        for iExp = 1:length(meta)
            if ~separate && iExp > 1
                a(expIndices(2:end)) = a(expIndices(1));
                b(expIndices(2:end)) = b(expIndices(1));
                break
            end
            a{expIndices(iExp)} = ones(1, size(meta(iExp).F,2));
            b{expIndices(iExp)} = zeros(1, size(meta(iExp).F,2));
            for iCell = 1:length(indChange)
                ind = ~any(isnan([meta(iExp).F_final(:,indChange(iCell)), ...
                    F2_final{expIndices(iExp)}(:,indChange(iCell))]), 2);
                B = [meta(iExp).F_final(ind,indChange(iCell)), ...
                    ones(sum(ind),1)] \ ...
                    F2_final{expIndices(iExp)}(ind,indChange(iCell));
                a{expIndices(iExp)}(indChange(iCell)) = B(1);
                b{expIndices(iExp)}(indChange(iCell)) = B(2);
            end
        end
        
        corrections(iSet).plane(iPlane).a = a;
        corrections(iSet).plane(iPlane).b = b;
        corrections(iSet).plane(iPlane).F_final = F2_final;
    end
end

save(fullfile(folderROIData, 'corrections.mat'), 'corrections', 'doCorrect')

%% Change assignment of experiments for datasets 7-18
for iSet = 7:18
    oldInds = corrections(iSet).expts;
    newInds = 1:length(corrections(iSet).expts);
    for iPlane = 1:length(corrections(iSet).plane)
        corrections(iSet).plane(iPlane).a(newInds(end)+1:end) = [];
        corrections(iSet).plane(iPlane).b(newInds(end)+1:end) = [];
        corrections(iSet).plane(iPlane).F_final(newInds(end)+1:end) = [];
    end
    corrections(iSet).expts = newInds;
end

%% Plot new and old traces
F_final = cell(1, length(meta));
for iExp = 1:length(meta)
    F_final{iExp} = meta(iExp).F_final;
end
F_final_all = cat(1, F_final{:});
F2_final_all = cat(1, F2_final{:});
for k = 1:20%length(indChange)
    figure('position',[8 494 1908 616])
    ax = [0 0];
    subplot(2,1,1)
    plot(F_final_all(:,indChange(k)))
    ax(1) = gca;
    subplot(2,1,2)
    plot(F2_final_all(:,indChange(k)))
    hold on
    f = [];
    for iExp=1:length(meta)
        f = [f; doCorrect(a{expIndices(iExp)}(indChange(k)), ...
            b{expIndices(iExp)}(indChange(k)), meta(iExp).F_final(:,indChange(k)))];
    end
    plot(f)
    ax(2) = gca;
    linkaxes(ax,'x')
end