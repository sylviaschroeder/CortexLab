%% Folders
% file location depending on PCs
folderPC = 'C:\STORAGE\OneDrive - University College London';
folderROIData = fullfile(folderPC, 'Lab\DATA\InfoStructs');

make_db_boutons;

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
    for iPlane = 1:length(db(iSet).planesToProcess)
        fprintf('  plane %d of %d\n', iPlane, length(db(iSet).planesToProcess))
        clear meta
        expIndices = zeros(1, length(db(iSet).expts));
        for iExp = 1:length(db(iSet).expts)
            expIndices(iExp) = db(iSet).expts(iExp);
%             expIndices(iExp) = iExp;
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
        else
            separate = false;
            newNpilSlopes = NpilSlopes{1};
        end
        
        % the target slope is between [0 2], adjust new slopes accordingly
        newNpilSlopes(newNpilSlopes < 0) = 0;
        newNpilSlopes(newNpilSlopes > 2) = 2;
        % make sure that new slope <= 0.9 * prctile(F,prctileNeuron)/prctile(Npil,prctileNpil) so
        % delta-F = F - slope*Npil >= F - 0.9*prctile(F,prctileNeuron)
        ratios = prctile(F, prctileNeuron) ./ prctile(Npil0_all, prctileNpil);
        ind = newNpilSlopes > 0.9 .* ratios;
        newNpilSlopes(ind) = 0.9 .* ratios(ind);
        
        m1 = cell(1, max(db(iSet).expts));
        m2 = cell(1, max(db(iSet).expts));
        m3 = cell(1, max(db(iSet).expts)); 
        if separate
            for iExp = 1:length(meta)
                m1{expIndices(iExp)} = max(1, nanmean(F0_corr{iExp}));
                m2{expIndices(iExp)} = max(1, nanmean(F0_corr{iExp} - ...
                    (newNpilSlopes - NpilSlopes{iExp}) .* Npil0{iExp}));
                R = 0;
                if isfield(meta, 'R_final')
                    R = meta(iExp).redFiltSlopes .* nanmean(meta(iExp).R_final);
                end
                m3{expIndices(iExp)} = (m1{expIndices(iExp)} - m2{expIndices(iExp)}) ...
                    .* R - (newNpilSlopes - NpilSlopes{iExp}) .* ...
                    nanmean(meta(iExp).Npil - Npil0{iExp});
            end
        else
            m1{expIndices} = max(1, nanmean(F0_corr_all,1));
            m2{expIndices} = max(1, nanmean(F0_corr_all - ...
                (newNpilSlopes - NpilSlopes{1}) .* Npil0_all));
            for iExp = 1:length(meta)
                R = 0;
                if isfield(meta, 'R_final')
                    R = meta(iExp).redFiltSlopes .* nanmean(meta(iExp).R_final);
                end
                m3{expIndices(iExp)} = (m1{expIndices(iExp)} - m2{expIndices(iExp)}) ...
                    .*R - (newNpilSlopes - NpilSlopes{1}) .* ...
                    nanmean(meta(iExp).Npil - Npil0{iExp});
            end
        end
        a = cellfun(@rdivide, m1, m2, 'UniformOutput', false);
        b = cellfun(@rdivide, m3, m2, 'UniformOutput', false);
        
        corrections(iSet).plane(iPlane).a = a;
        corrections(iSet).plane(iPlane).b = b;
    end
end

save(fullfile(folderROIData, 'corrections.mat'), 'corrections', 'doCorrect')