ops.ResultsSavePath     = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\F';
% ops.RootStorage         = '//zserver/Data/2P';
% ops.RootStorage         = '//zserver4/Data/2P';
% ops.RootStorage         = '//zserver.cortexlab.net/Data/Subjects';
% ops.RootStorage         = '//zubjects.cortexlab.net/Subjects';
ops.RootStorage         = 'J:\DATA\raw';
ops.InfoStorage         = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\InfoStructs';

%% Build db
% cd('C:\STORAGE\Dropbox\Code\scripts\preprocessing') % start this code in the directory with make_db
% make_db_nonvisual;
make_db_boutons;
% make_db_wheelTask;

%% Change filename of bin files in ops-struct of F_... file
% binFolder = 'J:\DATA\tempreg';
% for iSet = 1:length(db)
%     expDir = sprintf('%d_',db(iSet).expts);
%     expDir(end) = [];
%     regfile = fullfile(ops.ResultsSavePath, db(iSet).mouse_name, ...
%         db(iSet).date, expDir, sprintf('regops_%s_%s.mat', ...
%         db(iSet).mouse_name, db(iSet).date));
%     regdata = load(regfile);
%     for iPlane = union(1,db(iSet).planesToProcess)
%         Ffile = fullfile(ops.ResultsSavePath, db(iSet).mouse_name, ...
%             db(iSet).date, expDir, sprintf('F_%s_%s_plane%d_proc.mat', ...
%             db(iSet).mouse_name, db(iSet).date, iPlane));
%         filename = regdata.ops1{iPlane}.RegFile;
%         [~,name,ext] = fileparts(filename);
%         filename = fullfile(binFolder, [name, ext]);
%         if exist(Ffile, 'file')
%             data = load(Ffile);
%             data.dat.ops.RegFile = filename;
%             save(Ffile, '-struct', 'data');
%         end
%         regdata.ops1{iPlane}.RegFile = filename;
%     end
%     save(regfile, '-struct', 'regdata');
% end

%% Redo registration (green channel)
iCh = 1;
for iSet = 1 %1:length(db)
    expDir = sprintf('%d_',db(iSet).expts);
    expDir(end) = [];
    setDir = fullfile(ops.ResultsSavePath, db(iSet).mouse_name, ...
        db(iSet).date, expDir);
    
    list = dir(fullfile(setDir, sprintf('regops_%s_%s.mat', ...
        db(iSet).mouse_name, db(iSet).date)));
    list = {list.name};
    data = load(fullfile(setDir, list{1}));
    ops1 = data.ops1;
    redoRegistration(ops1, iCh, db(iSet).planesToProcess);
end

%% Get signal and neuropil
for iSet = 1:length(db)
    ops0 = build_ops3(db(iSet), ops);
    
    % Calculate neural signal and neuropil from movie or SVD components
    ops1 = ops0;
    ops1.inNeurop   = 1; % inner diameter of neuropil mask (3 for neurons, 1 for boutons)
    ops1.outNeurop  = 30; % outer diameter of neuropil mask (30 for neurons and boutons)
    ops1.microID    = db(iSet).microID; % ID of microscope (b: B-scope, b2: Bergamo2, m: MOM)
    ops1.processed  = 1;
    ops1.useSVD     = 0;
    ops1.getSignal  = 1;
    ops1.getNeuropil= 1;
    ops1.newFile    = 0;
    for iPlane = db(iSet).planesToProcess
        get_signals_and_neuropil(ops1, iPlane, 1);
    end
end

%% Register red channel and get traces
iCh = 2;
for iSet = 1:length(db)
    expDir = sprintf('%d_',db(iSet).expts);
    expDir(end) = [];
    setDir = fullfile(ops.ResultsSavePath, db(iSet).mouse_name, ...
        db(iSet).date, expDir);
    
    list = dir(fullfile(setDir, sprintf('regops_%s_%s.mat', ...
        db(iSet).mouse_name, db(iSet).date)));
    list = {list.name};
    data = load(fullfile(setDir, list{1}));
    ops1 = data.ops1;
    for iPlane = 1:length(ops1)
        % change RegFile to '..._red.bin'
        ops1{iPlane}.RegFile = [ops1{iPlane}.RegFile(1:end-4) '_red.bin'];
    end
    redoRegistration(ops1, iCh, db(iSet).planesToProcess);
    
%     % downsample movie if larger than 512 x 512
%     if data.ops.Lx == 1024
%         preproc.downsampleMovie(data.ops);
%         delete(data.ops.RegFile);
%         movefile([data.ops.RegFile(1:end-4) '_small.bin'], ...
%             data.ops.RegFile);
%     end
    
    for iPlane = db(iSet).planesToProcess
        data = load(fullfile(setDir, sprintf('F_%s_%s_plane%d_proc.mat', ...
            db(iSet).mouse_name, db(iSet).date, iPlane)));
        data.dat.ops.RegFile = [data.dat.ops.RegFile(1:end-4) '_red.bin'];
        path = strsplit(data.dat.ops.ResultsSavePath, '/');
        path(cellfun(@isempty,path)) = [];
        data.dat.ops.ResultsSavePath = fullfile(ops.ResultsSavePath, ...
            path{end-2:end}, 'red');
        data.dat.Fcell = {};
        data.dat.FcellNeu = {};
        if ~exist(data.dat.ops.ResultsSavePath, 'dir')
            mkdir(data.dat.ops.ResultsSavePath)
        end
        save(fullfile(data.dat.ops.ResultsSavePath, ...
            sprintf('F_%s_%s_plane%d_proc.mat', data.dat.ops.mouse_name, ...
            data.dat.ops.date, iPlane)), '-struct', 'data')
    end
    
    ops0 = build_ops3(db(iSet), ops);
    ops0.ResultsSavePath = fullfile(ops0.ResultsSavePath, 'red');
    
    % Calculate neural signal and neuropil from movie
    ops1 = ops0;
    ops1.inNeurop   = 1; % inner diameter of neuropil mask (1 for boutons, 3 for neurons)
    ops1.outNeurop  = 30; % outer diameter of neuropil mask
    ops1.microID    = db(iSet).microID; % ID of microscope (b: B-scope, b2: Bergamo2, m: MOM)
    ops1.processed  = 1;
    ops1.useSVD     = 0;
    ops1.getSignal  = 1;
    ops1.getNeuropil= 0;
    ops1.newFile    = 0;
    for iPlane = db(iSet).planesToProcess
        get_signals_and_neuropil(ops1, iPlane);
    end
end

%% Translate Suite2P output into meta structures and save
for iSet = 1:length(db)
    ops0 = build_ops3(db(iSet), ops);
    % Translate green and red into info structure
    ops2 = ops0;
    planesToProcess = db(iSet).planesToProcess;
    ops2.processed = 1;
    ops2.newFile   = 0;
    ops2.makeCompatibleWith = '-v2';
    nameEnd = '';
    if ops2.processed
        nameEnd = [nameEnd '_proc'];
    end
    if ops2.newFile
        nameEnd = [nameEnd '_new'];
    end
    nameEndRed = [nameEnd '_new.mat'];
    nameEnd = [nameEnd '.mat'];
    for iPlane = planesToProcess
        list = dir(fullfile(ops2.ResultsSavePath, ...
            sprintf(['F_%s_%s_plane%d' nameEnd], ...
            ops2.mouse_name, ops2.date, iPlane)));
        if length(list) ~= 1
            display('WARNING: files found:')
            list.names
            continue
        end
        filename = list(1).name;
        allInfos = preproc.dat2info(fullfile(ops2.ResultsSavePath, filename), ops2);
        allInfosRed = preproc.dat2info(fullfile(ops2.ResultsSavePath, 'red', filename), ops2);
        nInfos = length(allInfos);
        
        for iInfo = 1:nInfos
            meta = allInfos{iInfo};
            m = allInfosRed{iInfo};
            meta.R = m.F;
            targetFolder = fullfile(ops.InfoStorage, db(iSet).mouse_name, ...
                db(iSet).date, num2str(iInfo));
            targetFile = fullfile(targetFolder, ...
                sprintf('%s_%d_%s_2P_plane%03d_ROI.mat', db(iSet).date, ...
                iInfo, db(iSet).mouse_name, iPlane));
%             targetFile = fullfile(meta.folderProcessed, [meta.basenamePlane, '_ROI.mat']);
            if ~exist(targetFolder, 'dir')
                mkdir(targetFolder);
            end
            save(targetFile, 'meta');
        end
    end
end

% merge red and green meta structure if they exist already
% for iSet = 1:length(db)
%     for iExp = 1:length(db(iSet).expts)
%         for iPlane = 1:length(db(iSet).planesToProcess)
%             data = load(fullfile(ops.InfoStorage, db(iSet).mouse_name, ...
%                 db(iSet).date, num2str(db(iSet).expts(iExp)), ...
%                 sprintf('%s_%d_%s_2P_plane%03d_ROI.mat', db(iSet).date, ...
%                 db(iSet).expts(iExp), db(iSet).mouse_name, ...
%                 db(iSet).planesToProcess(iPlane))));
%             meta = data.meta;
%             data = load(fullfile(ops.InfoStorage, 'red', db(iSet).mouse_name, ...
%                 db(iSet).date, num2str(db(iSet).expts(iExp)), ...
%                 sprintf('%s_%d_%s_2P_plane%03d_ROI.mat', db(iSet).date, ...
%                 db(iSet).expts(iExp), db(iSet).mouse_name, ...
%                 db(iSet).planesToProcess(iPlane))));
%             meta.R = data.meta.R;
%         end
%     end
% end
   
%% Preprocess traces (each experiment separately)
% % parameters chosen on 19.12.2016
% timeWindow = 180; % in sec (for neuron imaging)
% prctileNeuron = 8;
% prctileNpil = 40;
% windowRed = 10; % in sec
% for iSet = 1:length(db)
%     fprintf('Dataset %d of %d\n', iSet, length(db))
%     for iExp = 1:length(db(iSet).expts)
%         fprintf('  experiment %d of %d\n', iExp, length(db(iSet).expts))
%         for iPlane = 1:length(db(iSet).planesToProcess)
%             fprintf('    plane %d of %d\n', iPlane, length(db(iSet).planesToProcess))
%             % load meta (with green and red traces)
%             data = load(fullfile(ops.InfoStorage, db(iSet).mouse_name, ...
%                 db(iSet).date, num2str(iExp), sprintf( ...
%                 '%s_%d_%s_2P_plane%03d_ROI.mat', db(iSet).date, ...
%                 iExp, db(iSet).mouse_name, ...
%                 db(iSet).planesToProcess(iPlane))));
%             meta = data.meta;
%             % get smapling rate
%             frameTimes = ppbox.getFrameTimes(meta);
%             samplingRate = 1/median(diff(frameTimes));
%             
%             % Correct for neuropil contamination
%             [F_filt, F0] = preproc.removeSlowDrift(meta.F, samplingRate, ...
%                 timeWindow, prctileNeuron);
%             [Npil_filt, Npil0] = preproc.removeSlowDrift(meta.Npil, ...
%                 samplingRate, timeWindow, prctileNpil);
%             [~,t] = min((F0-Npil0) ./ Npil0, [], 1);
%             meta.opsNpil.subtractDrift = 1;
%             meta.opsNpil.driftWindow = timeWindow;
%             meta.opsNpil.prctileNeuron = prctileNeuron;
%             meta.opsNpil.prctileNpil = prctileNpil;
%             opt.minNp = 5;
%             opt.maxNp = 95;
%             opt.constrainedFit = 1;
%             opt.verbose = 0;
%             meta.opsNpil.minNp = opt.minNp;
%             meta.opsNpil.maxNp = opt.maxNp;
%             meta.opsNpil.constrainedFit = opt.constrainedFit;
%             ind = sub2ind(size(F0), t, 1:size(F0,2));
%             [~, parameters] = estimateNeuropil(bsxfun(@plus, F_filt, F0(ind))', ...
%                 bsxfun(@plus, Npil_filt, Npil0(ind))', opt);
%             F_woNpil = meta.F - bsxfun(@times, parameters.corrFactor(:,2)', ...
%                 meta.Npil);          
%             meta.NpilSlopes = permute(parameters.corrFactor(:,2,:),[3 1 2]);
%             meta.NpilIntercepts = permute(parameters.corrFactor(:,1,:),[3 1 2]);
%             meta.F_woNpil = F_woNpil;
%             % remove slow drift from neuropil corrected traces
%             [deltaF,meta.F0] = preproc.removeSlowDrift(F_woNpil, samplingRate, ...
%                 timeWindow, prctileNeuron);
%             
%             % calculate delta-F-over-F
%             meta.F_delta = bsxfun(@rdivide, deltaF, max(1, mean(meta.F0, 1)));
% 
%             % for boutons only -------------------------------------------
% %             % remove exponential decay
% %             R = NaN(size(meta.R));
% %             ind = round(150 * samplingRate);
% %             for iCell = 1:size(meta.R,2)
% %                 y = meta.R(1:ind,iCell);
% %                 f = fit(frameTimes(1:ind)', y, ...
% %                     @(a,b,c,d,e,x) a+b.*exp(-x./c)+d.*exp(-x./e), ...
% %                     'Lower', [0 0 0 0 0], ...
% %                     'Upper', [max(y) max(y) 500 max(y) 500], ...
% %                     'StartPoint', [min(y) mean(y) 50 mean(y) 5]);
% %                 R(:,iCell) = meta.R(:,iCell) - f(frameTimes) + f.a;
% %             end
% %             % remove slow drift and filter red traces
% %             [R_noDrift, R0] = preproc.removeSlowDrift(R, ...
% %                 samplingRate, timeWindow, prctileNeuron);
%             % for neurons ------------------------------------------------
%             % remove slow drift and filter red traces
%             [R_noDrift, R0] = preproc.removeSlowDrift(meta.R, ...
%                 samplingRate, timeWindow, prctileNeuron);
%             % end (for neurons or boutons) -------------------------------
%             meta.redFiltWindow = windowRed;
%             nRF = round(windowRed * samplingRate);
%             R_final = medfilt1(R_noDrift, nRF, [], 1, 'includenan', 'truncate');
%             slopes = zeros(1, size(meta.F,2));
%             F_final = NaN(size(meta.F));
%             for c = 1:size(meta.F,2)
%                 if all(any(isnan([R_final(:,c) meta.F_delta(:,c)]),2))
%                     continue
%                 end
%                 mdl = fitlm(R_final(:,c), meta.F_delta(:,c), 'RobustOpts', 'on');
%                 slopes(c) = mdl.Coefficients.Estimate(2);
%                 F_final(:,c) = meta.F_delta(:,c) - slopes(c).*R_final(:,c);
%             end
%             meta.redFiltSlopes = slopes;
%             meta.R0 = R0;
%             meta.R_final = R_final;
%             meta.F_final = F_final;
%             
%             if isfield(meta, 'badFrames')
%                 meta.F_final(min(meta.badFrames, size(meta.F,1)),:) = NaN;
%             end
%             
%             if isfield(meta, 'Fcorr')
%                 meta = rmfield(meta, 'Fcorr');
%             end
%             
%             % save results
%             targetFolder = fullfile(ops.InfoStorage, db(iSet).mouse_name, ...
%                 db(iSet).date, num2str(iExp));
%             targetFile = fullfile(targetFolder, ...
%                 sprintf('%s_%d_%s_2P_plane%03d_ROI.mat', db(iSet).date, ...
%                 iExp, db(iSet).mouse_name, db(iSet).planesToProcess(iPlane)));
% %             targetFolder = meta.folderProcessed;
% %             targetFile = fullfile(meta.folderProcessed, [meta.basenamePlane, '_ROI.mat']);
%             if ~exist(targetFolder, 'dir')
%                 mkdir(targetFolder);
%             end
%             save(targetFile, 'meta');
%         end
%     end
% end

%% Preprocess traces across experiments!
% parameters chosen on 19.12.2016
timeWindow = 180; % in sec (for neuron imaging)
prctileNeuron = 8;
prctileNpil = 40;
windowRed = 10; % in sec
for iSet = 1:length(db)
    fprintf('Dataset %d of %d\n', iSet, length(db))
    for iPlane = 1:length(db(iSet).planesToProcess)
        fprintf('  plane %d of %d\n', iPlane, length(db(iSet).planesToProcess))
        clear meta
        for iExp = 1:length(db(iSet).expts)
            fprintf('    experiment %d of %d\n', iExp, length(db(iSet).expts))
            % load meta (with green and red traces)
            data = load(fullfile(ops.InfoStorage, db(iSet).mouse_name, ...
                db(iSet).date, num2str(iExp), sprintf( ...
                '%s_%d_%s_2P_plane%03d_ROI.mat', db(iSet).date, ...
                iExp, db(iSet).mouse_name, ...
                db(iSet).planesToProcess(iPlane))));
            meta(iExp) = data.meta;
        end
        % get smapling rate
        frameTimes = ppbox.getFrameTimes(meta(1));
        samplingRate = 1/median(diff(frameTimes));
            
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
        for iExp = 1:length(meta)
            [F_filt, F0] = preproc.removeSlowDrift(meta(iExp).F, samplingRate, ...
                timeWindow, prctileNeuron);
            F_filt_all = [F_filt_all; F_filt];
            F0_all = [F0_all; F0];
            [Npil_filt, Npil0] = preproc.removeSlowDrift(meta(iExp).Npil, ...
                samplingRate, timeWindow, prctileNpil);
            Npil_filt_all = [Npil_filt_all; Npil_filt];
            Npil0_all = [Npil0_all; Npil0];
            
            meta(iExp).opsNpil.subtractDrift = 1;
            meta(iExp).opsNpil.driftWindow = timeWindow;
            meta(iExp).opsNpil.prctileNeuron = prctileNeuron;
            meta(iExp).opsNpil.prctileNpil = prctileNpil;
            meta(iExp).opsNpil.minNp = opt.minNp;
            meta(iExp).opsNpil.maxNp = opt.maxNp;
            meta(iExp).opsNpil.constrainedFit = opt.constrainedFit;
        end
        % (2) find common neuropil correction parameters across experiments
        % and perform correction on each experiment
        [~,t] = min((F0_all-Npil0_all) ./ Npil0_all, [], 1);
        ind = sub2ind(size(F0_all), t, 1:size(F0_all,2));
        [~, parameters] = preproc.estimateNeuropil(bsxfun(@plus, F_filt_all, F0_all(ind))', ...
            bsxfun(@plus, Npil_filt_all, Npil0_all(ind))', opt);
        for iExp = 1:length(meta)
            F_woNpil = meta(iExp).F - bsxfun(@times, parameters.corrFactor(:,2)', ...
                meta(iExp).Npil);          
            meta(iExp).NpilSlopes = permute(parameters.corrFactor(:,2,:),[3 1 2]);
            meta(iExp).NpilIntercepts = permute(parameters.corrFactor(:,1,:),[3 1 2]);
            meta(iExp).F_woNpil = F_woNpil;
        end
        % (3) remove slow drift from neuropil corrected traces in each
        % experiment
        F0_all = [];
        deltaF = cell(1, length(meta));
        for iExp = 1:length(meta)
            [deltaF{iExp}, meta(iExp).F0] = preproc.removeSlowDrift( ...
                meta(iExp).F_woNpil, samplingRate, timeWindow, prctileNeuron);
            F0_all = [F0_all; meta(iExp).F0];
        end
        % (4) calculate delta-F-over-F with same denominator across
        % experiments
        for iExp = 1:length(meta)
            meta(iExp).F_delta = bsxfun(@rdivide, deltaF{iExp}, ...
                max(1, mean(F0_all, 1)));
        end

        % (5) remove slow drift of red traces for each experiment
        R_all = [];
        F_all = [];
        for iExp = 1:length(meta)
            % for boutons only (not adjusted for multiple experiments)-----
            % remove exponential decay
            R = NaN(size(meta(iExp).R));
            ind = round(150 * samplingRate);
            for iCell = 1:size(meta(iExp).R,2)
                y = meta(iExp).R(1:ind,iCell);
                f = fit((1:ind)', y, ...
                    @(a,b,c,d,e,x) a+b.*exp(-x./c)+d.*exp(-x./e), ...
                    'Lower', [0 0 0 0 0], ...
                    'Upper', [max(y) max(y) 500 max(y) 500], ...
                    'StartPoint', [min(y) mean(y) 50 mean(y) 5]);
                R(:,iCell) = meta(iExp).R(:,iCell) - f(1:size(R,1)) + f.a;
            end
            % remove slow drift and filter red traces
            [R_noDrift, meta(iExp).R0] = preproc.removeSlowDrift(R, ...
                samplingRate, timeWindow, prctileNeuron);
            % for neurons ------------------------------------------------
%             [R_noDrift, meta(iExp).R0] = preproc.removeSlowDrift(meta(iExp).R, ...
%                 samplingRate, timeWindow, prctileNeuron);
            % end for neurons or boutons -------------------------------
            meta(iExp).redFiltWindow = windowRed;
            nRF = round(windowRed * samplingRate);
            meta(iExp).R_final = medfilt1(R_noDrift, nRF, [], 1, 'includenan', 'truncate');
            R_all = [R_all; meta(iExp).R_final];
            F_all = [F_all; meta(iExp).F_delta];
        end
        % (6) fit red traces to green traces across all experiments at once
        slopes = zeros(1, size(meta(1).F,2));
        for c = 1:size(meta(1).F,2)
            if all(any(isnan([R_all(:,c) F_all(:,c)]),2))
                continue
            end
            mdl = fitlm(R_all(:,c), F_all(:,c), 'RobustOpts', 'on');
            slopes(c) = mdl.Coefficients.Estimate(2);
        end
        % (7) subtract weighted red traces from green traces
        for iExp = 1:length(meta)
            meta(iExp).F_final = meta(iExp).F_delta - slopes .* meta(iExp).R_final;
            meta(iExp).redFiltSlopes = slopes;
        end
        
        meta_all = meta;
        for iExp = 1:length(meta_all)
            meta = meta_all(iExp);
            if isfield(meta, 'badFrames')
                meta.F_final(min(meta.badFrames, size(meta.F,1)),:) = NaN;
            end
            
            % save results
            targetFolder = fullfile(ops.InfoStorage, db(iSet).mouse_name, ...
                db(iSet).date, num2str(iExp));
            targetFile = fullfile(targetFolder, ...
                sprintf('%s_%d_%s_2P_plane%03d_ROI.mat', db(iSet).date, ...
                iExp, db(iSet).mouse_name, db(iSet).planesToProcess(iPlane)));
%             targetFolder = meta.folderProcessed;
%             targetFile = fullfile(meta.folderProcessed, [meta.basenamePlane, '_ROI.mat']);
            if ~exist(targetFolder, 'dir')
                mkdir(targetFolder);
            end
            save(targetFile, 'meta');
        end
    end
end

%% Plot results of preprocessing
folder = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\Preprocessing\plots';
for iSet = 1:length(db)
    for iExp = 1 %:length(db(iSet).expts)
        for iPlane = 1:length(db(iSet).planesToProcess)
            % load meta (with green traces)
            data = load(fullfile(ops.InfoStorage, db(iSet).mouse_name, ...
                db(iSet).date, num2str(iExp), sprintf( ...
                '%s_%d_%s_2P_plane%03d_ROI.mat', db(iSet).date, ...
                iExp, db(iSet).mouse_name, ...
                db(iSet).planesToProcess(iPlane))));
            meta = data.meta;
            % get sampling rate
            frameTimes = ppbox.getFrameTimes(meta);
            % FOR BOUTONS -------------------------------------------------
            % load ops structure
%             data = load(fullfile(ops.ResultsSavePath, db(iSet).mouse_name, ...
%                 db(iSet).date, str, sprintf( ...
%                 'F_%s_%s_plane%d.mat', db(iSet).mouse_name, ...
%                 db(iSet).date, db(iSet).planesToProcess(iPlane))));
%             procOps = data.ops;
%             numPl = 5;
            % -------------------------------------------------------------
            
            numPl = 4;
            fPlots = fullfile(folder, [db(iSet).mouse_name '_' ...
                db(iSet).date '_' num2str(db(iSet).expts(iExp))], ...
                num2str(db(iSet).planesToProcess(iPlane)));
            if ~exist(fPlots, 'dir')
                mkdir(fPlots);
            end
            for c = 1:size(meta.F,2) %min(size(meta.F,2), 50)
                ax = [];
                
                figure('Position', [1 41 1920 1083])
                subplot(numPl,1,1)
                hold on
                plot(frameTimes, meta.F(:,c), 'k')
                plot(frameTimes, meta.F_woNpil(:,c), 'Color', [0 .5 .5])
                plot(frameTimes, meta.F0(:,c), 'Color', [0 1 1])
                legend('raw','npil. corr.','F0')
                axis tight
                ax(end+1) = gca;
                
                subplot(numPl,1,2)
                hold on
                plot(frameTimes, meta.F_delta(:,c), 'Color', [0 0 .5])
                plot(frameTimes, meta.F_final(:,c), 'Color', [0 .5 0])
                legend('\DeltaF/F','F final')
                axis tight
                ax(end+1) = gca;
                
                subplot(numPl,1,3)
                hold on
                plot(frameTimes, meta.R_final(:,c), 'Color', [.5 0 0])
                legend('R final')
                axis tight
                ax(end+1) = gca;
                
                subplot(numPl,1,4)
                hold on
                plot(frameTimes, meta.R(:,c), 'Color', [.5 0 .5])
                plot(frameTimes, meta.R0(:,c), 'Color', [1 0 1])
                legend('raw R','R0')
                axis tight
                ax(end+1) = gca;
                
                % FOR BOUTONS ---------------------------------------------
%                 subplot(numPl,1,5)
%                 imagesc(frameTimes([1 end]), [1 size(procOps.usedPlanes,2)], ...
%                     procOps.usedPlanes(sum(procOps ...
%                     .Nframes(1:iExp-1))+(1:procOps.Nframes(iExp)),:)')
%                 colormap gray
%                 ylabel('Used planes')
%                 axis tight
%                 ax(end+1) = gca;
                % ---------------------------------------------------------
                
                linkaxes(ax, 'x')
                
                fig = gcf;
                fig.PaperPositionMode = 'auto';
                print(fullfile(fPlots, sprintf('preprocessing_%03d.jpg', c)), ...
                    '-djpeg','-r0')
                
                close gcf
            end
        end
    end
end

%% Find duplicate cells
clear opt
opt.minCorr = 0.4; % neurons: 0.4 (was 0.35 before 01.10.2017), boutons: 0.5;
opt.maxDistXY = 5; % neurons: 5, boutons: 2;
opt.maxDistZ = 20; % neurons: 20, boutons: 3;
opt.filtWindow = 5; % in samples
for iSet = 1:length(db)
    folder = fullfile(ops.InfoStorage, db(iSet).mouse_name, db(iSet).date);
    files = [db(iSet).date '_%d_' db(iSet).mouse_name '_2P_plane%03d_ROI.mat'];
    planes = db(iSet).planesToProcess;
    expts = db(iSet).expts;
    ROIpos = cell(1, max(planes));
    F = cell(max(expts), max(planes));
    for iExp = 1:length(expts)
        for iPlane = planes
            data = load(fullfile(folder, num2str(iExp), ...
                sprintf(files, iExp, iPlane)));
            meta = data.meta;
            if iExp == 1
                ROIpos{iPlane} = cat(1, meta.ROI.CellXYZMicrons{:});
            end
            if iPlane == planes(1)
                frameTimes = ppbox.getFrameTimes(meta);
                if length(frameTimes) ~= size(meta.F_final,1)
                    fprintf('WARNING: number of timestamps is less than number of frames! Last frames ignored.\n');
                end
                F{iExp,iPlane} = meta.F_final(1:length(frameTimes),:);
            else
                tr = NaN(length(frameTimes), size(meta.F_final,2));
                ind = ~all(isnan(meta.F_final),1);
                t = ppbox.getFrameTimes(meta);
                if length(t) ~= size(meta.F_final,1)
                    fprintf('WARNING: number of timestamps is less than number of frames! Last frames ignored.\n');
                end
                tr(:,ind) = interp1(t, meta.F_final(1:length(t),ind), ...
                    frameTimes, 'pchip');
                F{iExp,iPlane} = tr;
            end
        end
    end
    minis = min(cellfun(@size, F(1:length(expts),planes), ...
        num2cell(ones(size(F(1:length(expts),planes))))), [], 2);
    F(1:length(expts),planes) = cellfun(@(x,y) (x(1:y,:)), F(1:length(expts),planes), ...
        num2cell(repmat(minis,1,length(planes))), 'UniformOutput', false);
%     F_ = cell(1, max(planes));
%     for iPlane = planes
%         F_{iPlane} = cat(1, F{:,iPlane});
%     end
    [isDuplicate, duplicates] = preproc.findDuplicateCells(ROIpos, F, opt);
    
    for iPlane = planes
        for iExp = 1:length(expts)
            file = fullfile(folder, num2str(iExp), sprintf(files, iExp, iPlane));
            data = load(file);
            meta = data.meta;
            meta.ROI.isDuplicate = isDuplicate{iPlane};
            meta.ROI.duplicates = duplicates(iPlane).ROI;
            save(file, 'meta');
        end
    end
end

%% Plot traces of duplicate cells
% for iSet = 1:length(db)
%     folder = fullfile(ops.InfoStorage, db(iSet).mouse_name, db(iSet).date);
%     files = [db(iSet).date '_%d_' db(iSet).mouse_name '_2P_plane%03d_ROI.mat'];
%     planes = db(iSet).planesToProcess;
%     expts = 1:length(db(iSet).expts);
%     for iExp = expts(1)
%         clear meta
%         for iPlane = planes
%             data = load(fullfile(folder, num2str(iExp), ...
%                 sprintf(files, iExp, iPlane)));
%             meta(iPlane) = data.meta;
%         end
%         frameTimes = ppbox.getFrameTimes(meta(planes(end)));
%         figure('position',[3 420 1915 630])
%         hold on
%         for iPlane = planes
%             for iCell = 1:size(meta(iPlane).F_final,2)
%                 if meta(iPlane).ROI.isDuplicate(iCell)==1 || ...
%                         isempty(meta(iPlane).ROI.duplicates(iCell).plane)
%                     continue
%                 end
%                 tr = NaN(length(frameTimes), max(planes));
%                 tr(:,iPlane) = meta(iPlane).F_final(1:length(frameTimes),iCell);
%                 for d = 1:length(meta(iPlane).ROI.duplicates(iCell).plane)
%                     plane = meta(iPlane).ROI.duplicates(iCell).plane(d);
%                     id = meta(iPlane).ROI.duplicates(iCell).ID(d);
%                     tr(:,plane) = meta(plane).F_final(1:length(frameTimes),id);
%                 end
% %                 tr = bsxfun(@rdivide, bsxfun(@minus, tr, nanmean(tr,1)), ...
% %                     nanstd(tr,0,1));
%                 offset = 0;
%                 pl = [];
%                 depths = [];
%                 for k = 1:size(tr,2)
%                     if all(isnan(tr(:,k)))
%                         continue
%                     end
%                     col = zeros(1,3);
%                     if k ~= iPlane
%                         col = ones(1,3).*0.5;
%                     end
%                     maxi = max(tr(:,k));
%                     mini = min(tr(:,k));
%                     plot(frameTimes,tr(:,k)+offset-maxi, 'Color', col)
%                     depths(end+1) = offset - .5*(maxi-mini);
%                     pl(end+1) = k;
%                     offset = offset - (maxi-mini);
%                 end
%                 depths = flip(depths);
%                 pl = flip(pl);
%                 set(gca, 'YTick', depths, 'YTickLabel', pl,'box','off')
%                 title(sprintf('%s %s %d, plane %d, cell %d', ...
%                     db(iSet).mouse_name, db(iSet).date, iExp, iPlane, ...
%                     iCell), 'interpreter', 'none')
%                 ylabel('F_{final}')
%                 axis tight
%                 
%                 pause
%                 cla
%             end
%         end
%         close gcf
%     end
% end

%% Find cells with sudden and long lasting high calcium values
% CHANGE: have same cells as switch-on across all experiments
clear opt
opt.switchOnTime = 25; % in sec (should be shorter than half of timeWindow given in "Preprocess traces"
opt.threshold = .6;
opt.maxRatio = -0.1;
for iSet = 1:length(db)
    folder = fullfile(ops.InfoStorage, db(iSet).mouse_name, db(iSet).date);
    files = [db(iSet).date '_%d_' db(iSet).mouse_name '_2P_plane%03d_ROI.mat'];
    planes = db(iSet).planesToProcess;
    expts = 1:length(db(iSet).expts);
    for iPlane = planes
        for iExp = expts
            data = load(fullfile(folder, num2str(iExp), ...
                sprintf(files, iExp, iPlane)));
            meta = data.meta;
            frameTimes = ppbox.getFrameTimes(meta);
            sr = 1 / median(diff(frameTimes));
            switchOnSamples = round(opt.switchOnTime * sr);
            traces = medfilt1(meta.F_woNpil(2:end,:), 3, [], 1) - meta.F0(2:end,:);
            minis = min(traces, [], 1);
            maxis = max(traces, [], 1);
            threshs = minis + opt.threshold * (maxis - minis);
            isSwitchOn = zeros(size(traces,2), 1);
            for iCell = 1:size(traces,2)
                if isnan(threshs(iCell))
                    continue
                end
                ind = find(traces(:,iCell) < threshs(iCell));
                maxPeriod = max(diff(ind));
                medRatio = (median(traces(:,iCell)) - ...
                    (minis(iCell)+ 0.5*(maxis(iCell) - minis(iCell)))) / ...
                    (maxis(iCell) - minis(iCell));
                if maxPeriod > switchOnSamples && medRatio < opt.maxRatio
                    isSwitchOn(iCell) = 1;
                end
            end
            meta.ROI.isSwitchOn = isSwitchOn;
            save(fullfile(folder, num2str(iExp), ...
                sprintf(files, iExp, iPlane)), 'meta')
        end
    end
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

%% Plot raw and processed Calcium traces to check results (OLD)
% figure('Position', [3 260 1915 835]);
% for iSet = 1:length(db)
%     folder = fullfile(ops.InfoStorage, db(iSet).mouse_name, db(iSet).date);
%     files = [db(iSet).date '_%d_' db(iSet).mouse_name '_2P_plane%03d_ROI.mat'];
%     planes = db(iSet).planesToProcess;
%     expts = db(iSet).expts;
%     for iPlane = planes
%         F = [];
%         Npil = [];
%         Fcorr = [];
%         deltaF =[];
%         time = [];
%         t = 0;
%         for iExp = expts
%             file = fullfile(folder, num2str(iExp), sprintf(files, iExp, iPlane));
%             data = load(file);
%             meta = data.meta;
%             F = [F; meta.F];
%             Npil = [Npil; meta.Npil];
%             Fcorr = [Fcorr; meta.Fcorr];
%             deltaF = [deltaF; bsxfun(@rdivide, (meta.Fcorr - meta.F0), ...
%                 max(1, mean(meta.F0,1)))];
%             frameTimes = ppbox.getFrameTimes(meta);
%             time = [time; frameTimes(:)+t];
%             t = t + frameTimes(end);
%         end
%         fprintf('%s %s plane %d\n', ...
%             db(iSet).mouse_name, db(iSet).date, iPlane)
%         for iCell = 1:size(F,2)
%             if meta.ROI.isDuplicate(iCell) == 1
%                 continue
%             end
%             subplot(2,1,1)
%             hold off
%             plot(time, F(:,iCell))
%             hold on
%             d = prctile(F(:,iCell),3);
%             plot(time, Npil(:,iCell)+d-prctile(Npil(:,iCell),97))
%             d = d-diff(prctile(Npil(:,iCell),[3 97]));
%             plot(time, Fcorr(:,iCell)+d-prctile(Fcorr(:,iCell),97));
%             ylim([d-(prctile(Fcorr(:,iCell),97)-min(Fcorr(:,iCell))) max(F(:,iCell))])
%             xlim(time([1 end]))
%             legend('raw F','Npil','corr. F')
%             title(sprintf('%s %s plane %d neuron %d', ...
%                 db(iSet).mouse_name, db(iSet).date, iPlane, iCell), ...
%                 'Interpreter', 'none')
%             ax1 = gca;
%             subplot(2,1,2)
%             plot(time, deltaF(:,iCell))
%             xlabel('Time (in s)')
%             xlim(time([1 end]))
%             legend('\DeltaF/F_0')
%             linkaxes([ax1 gca], 'x')
%             pause
%         end
%     end
% end