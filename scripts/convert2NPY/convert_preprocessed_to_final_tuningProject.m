%% Folders
folderSource = 'Z:\ProcessedData';
folderSave = 'D:\Data\Running2pPaper';
folderTarget = 'C:\Users\Sylvia\OneDrive - University of Sussex\Projects\2025_RunningChangesTuning_SC_2P';
folderTools = 'C:\dev\toolboxes';
foldersRepo = 'C:\dev\workspaces\CortexLab';

%% Parameters
setNames = {'boutons', 'neurons'};
zooms = [1, 1.5, 2, 4, 8, 16];
fovSize = [730, 490, 360, 180, 95, 50];
f_zoom = fit(zooms', fovSize', 'linear');

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(genpath(fullfile(foldersRepo)))

%% Setup Python packages
np = py.importlib.import_module('numpy');

%% Copy data
% logFile = fopen(fullfile(folderSave, 'Data', 'log_copyData.txt'), 'w');
% for set = 1:length(setNames)
%     recs = readtable(fullfile(folderTarget, 'DataDescriptions', ...
%         sprintf("fitting_%s.csv", setNames{set})));
%     for rec = 1:size(recs,1)
%         fprintf(logFile, "%s %s\n", recs.Subject{rec}, recs.Date(rec));
%         fs = fullfile(folderSource, recs.Subject{rec}, ...
%             string(recs.Date(rec)));
%         ft = fullfile(folderSave, "Data", setNames{set}, ...
%             recs.Subject{rec}, string(recs.Date(rec)));
%         if ~isfolder(ft)
%             mkdir(ft);
%         end
% 
%         %% Neural data
%         % copy ROI information
%         pl = readNPY(fullfile(fs, "rois.planes.npy"));
%         if isfile(fullfile(fs, "rois.xyz.npy"))
%             xyz = readNPY(fullfile(fs, "rois.xyz.npy"));
%             % original position: (1) vertical pos (row of pixels 
%             % transformed to um), (2) horizontal pos (column of pixels 
%             % transformed to um), (3) depth (in um)
%             xyz = xyz(:, [2 1 3]);
%         else
%             xyz = NaN(length(pl),1);
%             fprintf(logFile, "  NOTE! File rois.xyz.npy missing.\n");
%         end
%         if isfile(fullfile(fs, "rois.id.npy"))
%             id = readNPY(fullfile(fs, "rois.id.npy"));
%         else
%             id = NaN(length(pl),1);
%             fprintf(logFile, "  NOTE! File rois.id.npy missing.\n");
%         end
%         % sz = [size(id,1), size(pl,1)];
%         sz = [size(id,1), size(pl,1), size(xyz,1)];
%         if set == 2
%             prof = readNPY(fullfile(fs, "rois.zprofiles.npy"));
%             if size(prof,2) == size(pl,1)
%                 if size(prof,1) == size(pl,1)
%                     fprintf(logFile, "  Cannot resolve correct format of zprofiles (rows == columns).\n");
%                 else
%                     fprintf(logFile, "  rois.zprofiles: dimensions swapped!\n");
%                     prof = prof';
%                 end
%             end
%             sz = [sz, size(prof,1)];
%             if isfile(fullfile(fs, "rois.isInh.npy"))
%                 isInh = readNPY(fullfile(fs, "rois.isInh.npy"));
%                 sz = [sz, size(isInh,1)];
%             elseif isfile(fullfile(fs, "rois.isInhibitory.npy"))
%                 isInh = readNPY(fullfile(fs, "rois.isInhibitory.npy"));
%                 sz = [sz, size(isInh,1)];
%             end
%         end
%         if ~all(sz(1) == sz(2:end))
%             fprintf(logFile, "  Unequal number of rows: roi files\n");
%             return
%         end
%         writeNPY(id, fullfile(ft, "2pRois.ids.npy"));
%         writeNPY(pl, fullfile(ft, "2pRois.2pPlanes.npy"));
%         if set == 2
%             writeNPY(prof, fullfile(ft, "2pRois.zProfiles.npy"));
%             if isfile(fullfile(fs, "rois.isInh.npy")) || ...
%                     isfile(fullfile(fs, "rois.isInhibitory.npy"))
%                 writeNPY(isInh, fullfile(ft, "2pRois.isInhibitory.npy"));
%             end
%         end
%         clear prof isInh
% 
%         % check ROI planes with suite2p plane folders
%         planes = unique(pl)';
%         dirsSuite2p = dir(fullfile(fs, "suite2p", "plane*"));
%         delay = readNPY(fullfile(fs, "planes.delay.npy"));
%         if set == 2
%             planesFolders = NaN(1, length(dirsSuite2p));
%             for p = 1:length(dirsSuite2p)
%                 planesFolders(p) = str2double(dirsSuite2p(p).name(6:end));
%             end
%             if ~all(ismember(planes, planesFolders))
%                 fprintf(logFile, '  NOTE! ROI planes do not match suite2p plane folders!\n');
%                 if length(planes) == length(planesFolders)
%                     fprintf(logFile, ['  --> Replace ROI planes (' num2str(planes) ...
%                         ') by folder planes (' num2str(planesFolders) ')\n']);
%                     pl_new = NaN(size(pl));
%                     for p = 1:length(planes)
%                         pl_new(pl == planes(p)) = planesFolders(p);
%                     end
%                     writeNPY(pl_new, fullfile(ft, "2pRois.2pPlanes.npy"))
%                     pl = pl_new;
%                     planes = unique(pl)';
%                 else
%                     fprintf(logFile, '  --> Correct ROI planes cannot be inferred. --> IGNORE DATASET!\n');
%                     continue
%                 end
%             end
%             zCorr = readNPY(fullfile(fs, "planes.zcorrelation.npy"));
%             zTrace = readNPY(fullfile(fs, "planes.zTrace.npy"));
%             if size(delay,1) > size(zCorr,1)
%                 indPl = ismember(planesFolders, planes);
%                 if ~all(indPl) && size(zCorr,1) == length(planes)
%                     zCorr_new = NaN([length(delay), size(zCorr,2:3)]);
%                     zCorr_new(indPl,:,:) = zCorr;
%                     zCorr = zCorr_new;
%                     zTrace_new = NaN([length(delay), size(zTrace,2)]);
%                     zTrace_new(indPl,:) = zTrace;
%                     zTrace = zTrace_new;
%                 elseif ~ismember(0, planesFolders) || ...
%                         length(delay) == size(zCorr,1)+1
%                     zCorr = padarray(zCorr, 1, NaN, "pre");
%                     zTrace = padarray(zTrace, 1, NaN, "pre");
%                 else
%                     fprintf(logFile, '  Mismatch between #planes in delays and zcorrelation cannot be corrected. --> IGNORE DATASET!\n');
%                     continue
%                 end
%             end
%             sz = [size(delay,1), size(zCorr,1), size(zTrace,1)];
%             if ~all(sz(1) == sz(2:end))
%                 fprintf(logFile, "  Unequal number of rows: plane files\n");
%                 return
%             end
%             writeNPY(delay, fullfile(ft, "2pPlanes.delay.npy"));
%             writeNPY(zCorr, fullfile(ft, "2pPlanes.zCorrelations.npy"));
%             writeNPY(zTrace, fullfile(ft, "2pPlanes.zTraces.npy"));
%         end
%         clear zCorr zTrace
% 
%         % copy calcium traces
%         dff = readNPY(fullfile(fs, "calcium.dff.npy"));
%         t = readNPY(fullfile(fs, "calcium.timestamps.npy"));
%         sz = [size(dff,1), size(t,1)];
%         if ~all(sz(1) == sz(2:end))
%             if sz(2) == sz(1) + 1
%                 t(end) = [];
%             else
%                 fprintf(logFile, "  Unequal number of rows: calcium files --> ignore dataset\n");
%                 rmdir(ft)
%                 continue
%             end
%         end
%         if set == 1 % boutons
%             t = t + delay(planes + 1);
%         end
%         writeNPY(dff, fullfile(ft, "2pCalcium.dff.npy"));
%         writeNPY(t, fullfile(ft, "2pCalcium.timestamps.npy"));
%         clear dff t
% 
%         %% Suite2p data
%         xPix = cell(length(id),1);
%         yPix = cell(length(id),1);
%         zPos = NaN(length(id),1);
%         for p = planes
%             if set == 1 % boutons
%                 ind = strcmp({dirsSuite2p.name}, "plane") | ...
%                     strcmp({dirsSuite2p.name}, "plane8");
%                 index_plane = 1;
%             else
%                 ind = strcmp({dirsSuite2p.name}, sprintf("plane%d", p));
%                 index_plane = p;
%             end
%             dSuite2p = dirsSuite2p(ind);
%             stat_py = np.load(fullfile(dSuite2p.folder, dSuite2p.name, ...
%                 "stat.npy"), pyargs("allow_pickle", true));
%             stat_py = pyrun("stats = c.tolist()", "stats", c = stat_py);
%             ops_py = np.load(fullfile(dSuite2p.folder, dSuite2p.name, ...
%                 "ops.npy"), pyargs("allow_pickle", true));
%             ops_py = pyrun("ops = c.item()", "ops", c = ops_py);
% 
%             mi = double(ops_py{"meanImg"});
%             if p == planes(1)
%                 % check whether red channel was used
%                 ch = double(ops_py{"nchannels"});
%                 if ch == 1
%                     hasRedChannel = 0;
%                 else
%                     hasRedChannel = 1;
%                 end
%                 % use zoom to calculate size of FOV in um
%                 tiffPath = ops_py{"filelist"};
%                 tiffPath = char(tiffPath{end});
%                 tiffPath(1) = 'Z';
%                 tiffObj = Tiff(tiffPath);
%                 settings = getTag(tiffObj, 'Artist');
%                 zoom = split(regexp(settings, 'scanZoomFactor": [\d]', 'match'));
%                 zoom = str2double(zoom{end});
%                 fovSizeMicrons = [1 1] .* f_zoom(zoom);
%                 % initialize mean image of planes
%                 if set == 1 % boutons
%                     meanImg = NaN([1, size(mi)]);
%                 else
%                     meanImg = NaN([length(delay), size(mi)]);
%                 end
%                 % get imaging depths across planes
%                 [piezoDepth, piezoPlane] = spatial.getPiezoMovement(ops_py);
%             end
%             % mean image for each plane
%             meanImg(index_plane,:,:) = permute(mi, [3 1 2]);
% 
%             ind_rois = find(pl == p)';
%             % get number of lines per frame
%             numLines = double(ops_py{'Ly'});
%             for k = ind_rois
%                 roi = stat_py{id(k) + 1};
%                 % ROI masks
%                 xPix{k} = double(roi{"xpix"});
%                 yPix{k} = double(roi{"ypix"});
%                 % ROI depth
%                 yPos = cell(roi{'med'});
%                 yPos = double(yPos{1});
%                 relativeVolPos = p + yPos / numLines;
%                 zPos(k) = interp1(piezoPlane, piezoDepth, relativeVolPos);
%             end
%         end
%         maxLen = max(cellfun(@length, xPix));
%         maskX = NaN(length(id), maxLen);
%         maskY = NaN(length(id), maxLen);
%         for k = 1:length(id)
%             maskX(k,1:length(xPix{k})) = xPix{k};
%             maskY(k,1:length(yPix{k})) = yPix{k};
%         end
%         zPos = zPos + recs.Depth(rec);
%         xyz(:,3) = zPos;
% 
%         writeNPY(xyz, fullfile(ft, "2pRois.xyz.npy"));
%         writeNPY(maskX, fullfile(ft, "2pRois.maskX.npy"));
%         writeNPY(maskY, fullfile(ft, "2pRois.maskY.npy"));
%         writeNPY(meanImg, fullfile(ft, "2pPlanes.meanFrame.npy"));
%         writeNPY(hasRedChannel, fullfile(ft, "recordings.hasRedChannel.npy"));
%         writeNPY(fovSizeMicrons, fullfile(ft, "recordings.fovSizeMicrons.npy"));
% 
%         %% Stimulus data
%         % copy grating information
%         if isfile(fullfile(fs, "gratings.startTime.npy"))
%             t0 = readNPY(fullfile(fs, "gratings.startTime.npy"));
%             t1 = readNPY(fullfile(fs, "gratings.endTime.npy"));
%             dirs = readNPY(fullfile(fs, "gratings.direction.npy"));
%             sf = readNPY(fullfile(fs, "gratings.spatialF.npy"));
%             tf = readNPY(fullfile(fs, "gratings.temporalF.npy"));
%             contrast = readNPY(fullfile(fs, "gratings.contrast.npy"));
%             reward = readNPY(fullfile(fs, "gratings.reward.npy"));
%             if ~all(isnan(reward))
%                 fprintf(logFile, "  Reward was given during gratings!\n");
%                 return
%             end
%             if isfile(fullfile(fs, "gratings.st.updated.npy"))
%                 t0 = readNPY(fullfile(fs, "gratings.st.updated.npy"));
%                 t1 = readNPY(fullfile(fs, "gratings.et.updated.npy"));
%                 if isfile(fullfile(fs, "gratings.direction.updated.npy"))
%                     dirs = readNPY(fullfile(fs, "gratings.direction.updated.npy"));
%                     sf = readNPY(fullfile(fs, "gratings.spatialF.updated.npy"));
%                     tf = readNPY(fullfile(fs, "gratings.temporalF.updated.npy"));
%                     contrast = readNPY(fullfile(fs, "gratings.contrast.updated.npy"));
%                 end
%             end
%             sz = [size(t0,1), size(t1,1), size(dirs,1), size(sf,1), ...
%                 size(tf,1), size(contrast,1)];
%             if ~all(sz(1) == sz(2:end))
%                 fprintf(logFile, "  Unequal number of rows: grating files\n");
%                 return
%             end
%             intv = readNPY(fullfile(fs, "gratingsExp.intervals.npy"));
%             if size(intv,2) == 1 % was incorrectly saved into single vector
%                 intv = reshape(intv, 2, [])';
%             end
%             writeNPY([t0 t1], fullfile(ft, "grating.intervals.npy"));
%             writeNPY(dirs, fullfile(ft, "grating.directions.npy"));
%             writeNPY(sf, fullfile(ft, "grating.spatialFrequencies.npy"));
%             writeNPY(tf, fullfile(ft, "grating.temporalFrequencies.npy"));
%             writeNPY(contrast, fullfile(ft, "grating.contrasts.npy"));
%             writeNPY(intv, fullfile(ft, "gratingExp.intervals.npy"));
%             clear t0 t1 dirs sf tf contrast reward intv
%         end
% 
%         % copy sparse noise information
%         if isfile(fullfile(fs, "sparse.startTime.npy"))
%             t = readNPY(fullfile(fs, "sparse.startTime.npy"));
%             map = readNPY(fullfile(fs, "sparse.map.npy"));
%             edges = readNPY(fullfile(fs, "sparseExp.edges.npy"))';
%             intv = readNPY(fullfile(fs, "sparseExp.intervals.npy"));
%             if size(intv, 2) ~= 2
%                 fprintf(logFile, "  Incorrect interval format: sparse noise files\n");
%                 return
%             end
%             if size(edges,2) < 4 % only number of rows and columns were saved
%                 fprintf(logFile, "  Sparse noise: position of stimulus edges missing. -> Save NaNs.\n");
%                 edges = NaN(size(edges,1), 4);
%             else
%                 if size(intv,1) > 1 && size(edges,1) == 1
%                     edges = reshape(edges', [], size(intv,1))';
%                 end
%                 if size(t,1) ~= size(map,1) || size(edges,1) ~= size(intv,1)
%                     fprintf(logFile, "  Unequal number of rows: sparse noise files\n");
%                     return
%                 end
%                 if size(edges,2) > 6
%                     fprintf(logFile, "  Sparse noise: stimulus edge data incorrect.\n");
%                     return
%                 end
%                 % original edges: [rows, columns, bottom, top, left, right]
%                 % or [columns, bottom, top, left, right]
%                 % where bottom and top are positive if above horizon;
%                 % transform to: [left, right, top, bottom] where bottom and
%                 % top are negative if above horizon
%                 edges = edges(:, end-3 : end);
%                 edges = [edges(:, [3 4]) -edges(:, [2 1])];
%             end
%             map(map == 0) = -1;
%             map(map == 0.5) = 0;
%             writeNPY(t, fullfile(ft, "sparseNoise.times.npy"));
%             writeNPY(map, fullfile(ft, "sparseNoise.map.npy"));
%             writeNPY(edges, fullfile(ft, "sparseNoiseExp.edges.npy"));
%             writeNPY(intv, fullfile(ft, "sparseNoiseExp.intervals.npy"));
%             clear t map edges intv
%         end
% 
%         % copy full-field stimulus information (chirps)
%         if isfile(fullfile(fs, "fullField.startTime.npy"))
%             t = readNPY(fullfile(fs, "fullField.startTime.npy"));
%             descr_py = np.load(fullfile(fs, "fullField.stim.npy"), ...
%                 pyargs("allow_pickle", true));
%             numStim = double(descr_py.size);
%             stimDescr = cell(numStim, 1);
%             for st = 1:numStim
%                 descr = pyrun("b = c[ind,0]", "b", c = descr_py, ind = int32(st-1));
%                 stimDescr{st} = char(descr);
%             end
%             if size(t,1) ~= numStim
%                 fprintf(logFile, "  Unequal number of rows: fullField files\n");
%                 return
%             end
%             intv = readNPY(fullfile(fs, "fullFieldExp.intervals.npy"));
%             if size(intv, 2) ~= 2
%                 fprintf(logFile, "  Incorrect interval format: fullField files\n");
%                 return
%             end
%             writeNPY(t, fullfile(ft, "fullField.times.npy"));
%             writecell(stimDescr, fullfile(ft, "fullField.stimNames.csv"))
%             writeNPY(intv, fullfile(ft, "fullFieldExp.intervals.npy"));
%             clear t descr_py numStim stimDescr descr st intv
%         end
% 
%         % copy circle information
%         if isfile(fullfile(fs, "circles.startTime.npy"))
%             t = readNPY(fullfile(fs, "circles.startTime.npy"));
%             if isfile(fullfile(fs, "circles.startTime.updated.npy"))
%                 t = readNPY(fullfile(fs, "circles.startTime.updated.npy"));
%             end
%             diam = readNPY(fullfile(fs, "circles.diameter.npy"));
%             x = readNPY(fullfile(fs, "circles.x.npy"));
%             y = readNPY(fullfile(fs, "circles.y.npy"));
%             isWhite = readNPY(fullfile(fs, "circles.isWhite.npy"));
%             if isfile(fullfile(fs, "circles.diameter.updated.npy"))
%                 diam = readNPY(fullfile(fs, "circles.diameter.updated.npy"));
%                 x = readNPY(fullfile(fs, "circles.x.updated.npy"));
%                 y = readNPY(fullfile(fs, "circles.y.updated.npy"));
%                 isWhite = readNPY(fullfile(fs, "circles.isWhite.updated.npy"));
%             end
%             sz = [size(t,1), size(diam,1), size(x,1), size(y,1), ...
%                 size(isWhite,1)];
%             if ~all(sz(1) == sz(2:end))
%                 fprintf(logFile, "  Unequal number of rows: circle files\n");
%                 fprintf(logFile, "  %d ", sz);
%                 fprintf(logFile, "\n");
%             else
%                 intv = readNPY(fullfile(fs, "circlesExp.intervals.npy"));
%                 if size(intv, 2) ~= 2
%                     fprintf(logFile, "  Incorrect interval format: circles files\n");
%                     return
%                 end
%                 writeNPY(t, fullfile(ft, "circles.times.npy"));
%                 writeNPY(diam, fullfile(ft, "circles.diameters.npy"));
%                 writeNPY(x, fullfile(ft, "circles.xPos.npy"));
%                 writeNPY(y, fullfile(ft, "circles.yPos.npy"));
%                 writeNPY(isWhite, fullfile(ft, "circles.isWhite.npy"));
%                 writeNPY(intv, fullfile(ft, "circlesExp.intervals.npy"));
%             end
%             clear t diam x y isWhite intv
%         end
% 
%         %% Behaviour data
%         % copy running data (also interpolate data to regular time points)
%         speed = readNPY(fullfile(fs, "wheel.velocity.npy"));
%         t = readNPY(fullfile(fs, "wheel.timestamps.npy"));
%         l = [length(t), length(speed)];
%         if l(1) ~= l(2)
%             fprintf(logFile, '  NOTE! Unequal length of entries: running data\n');
%         end
%         dt = diff(t);
%         t_new = (min(t):median(dt, "omitnan"):max(t))';
%         % % find gaps in data, remove them from new timestamps
%         % gaps = find(dt > 3 * median(dt));
%         % for g = 1:length(gaps)
%         %     ind1 = find(t_new > t(gaps(g)), 1);
%         %     ind2 = find(t_new < t(gaps(g)+1), 1, 'last');
%         %     t_new(ind1:ind2) = [];
%         % end
%         mini = min(length(t), length(speed));
%         speed = general.interpolateWithNaNs(t(1:mini), speed(1:mini), ...
%             t_new, 'linear');
%         writeNPY(t_new, fullfile(ft, "running.timestamps.npy"))
%         writeNPY(speed, fullfile(ft, "running.speed.npy"))
% 
%         % copy eye data
%         if isfile(fullfile(fs, "eye.timestamps.npy"))
%             if ~isfile(fullfile(fs, "eye.diameter.npy"))
%                 fprintf(logFile, "  NOTE! Eye data exists (timestamps) but wasn't processed! Ignored for now.\n");
%             else
%                 t = readNPY(fullfile(fs, "eye.timestamps.npy"));
%                 d = readNPY(fullfile(fs, "eye.diameter.npy"));
%                 pos = readNPY(fullfile(fs, "eye.xyPos.npy"));
%                 dt = diff(t);
%                 t_new = (min(t):median(dt, "omitnan"):max(t))';
%                 % % find gaps in data, remove them from new timestamps
%                 % gaps = find(dt > 3 * median(dt));
%                 % for g = 1:length(gaps)
%                 %     ind1 = find(t_new > t(gaps(g)), 1);
%                 %     ind2 = find(t_new < t(gaps(g)+1), 1, 'last');
%                 %     t_new(ind1:ind2) = [];
%                 % end
%                 l = [length(t), length(d), length(pos)];
%                 if ~all(l(1) == l(2:3))
%                     fprintf(logFile, '  NOTE! Unequal length of entries: eye data\n');
%                 end
%                 mini = min(l);
%                 d = general.interpolateWithNaNs(t(1:mini), ...
%                     d(1:mini), t_new, 'linear');
%                 pos = general.interpolateWithNaNs(t(1:mini), ...
%                     pos(1:mini,:), t_new, 'linear');
%                 writeNPY(d, fullfile(ft, "eye.diameter.npy"))
%                 writeNPY(pos, fullfile(ft, "eye.xyPos.npy"))
%                 writeNPY(t_new, fullfile(ft, "eye.timestamps.npy"))
% 
%                 dirs = dir(fullfile(fs, "pupil", "xyPos_diameter"));
%                 mmPerPx = NaN(length(dirs), 1);
%                 for d = 1:length(dirs)
%                     f = fullfile(dirs(d).folder, dirs(d).name, ...
%                         "eyeVideo.mmPerPixel.npy");
%                     if startsWith(dirs(d).name, '.') || ~isfile(f)
%                         continue
%                     end
%                     mmPerPx(d) = readNPY(f);
%                 end
%                 writeNPY(median(mmPerPx, 1, "omitnan"), ...
%                     fullfile(ft, "eyeVideo.mmPerPixel.npy"))
%             end
%         else
%             fprintf(logFile, "  NOTE! Eye data does not exist!\n");
%         end
%     end
% end
% fclose(logFile);

%% Test validity of data
logFile = fopen(fullfile(folderSave, 'Data', 'log_testData.txt'), 'w');
binSizeGratings = 0.25;
fprintf(logFile, '\nTesting data validity\n');

fprintf(logFile, '\nGratings\n');
for set = 1:length(setNames)
    dirs_animals = dir(fullfile(folderSave, "Data", setNames{set}));
    for animal = 1:length(dirs_animals)
        if startsWith(dirs_animals(animal).name, '.')
            continue
        end
        dirs_dates = dir(fullfile(dirs_animals(animal).folder, ...
            dirs_animals(animal).name));
        for date = 1:length(dirs_dates)
            if startsWith(dirs_dates(date).name, '.')
                continue
            end
            ft = fullfile(folderSave, "Data", setNames{set}, ...
                dirs_animals(animal).name, dirs_dates(date).name);

            if isfile(fullfile(ft, "grating.intervals.npy"))
                t = readNPY(fullfile(ft, "grating.intervals.npy"));
                durs = diff(t, 1, 2);
                c = histcounts(durs, -binSizeGratings/2 : binSizeGratings : 5);
                ind = c < 0.1 * sum(c);
                if sum(c(ind)) > 0
                    fprintf(logFile, "  %s %s\n", dirs_animals(animal).name, ...
                        dirs_dates(date).name);
                    fprintf(logFile, '    %d of %d stimuli outside tolerated time range.\n', ...
                        sum(c(ind)), sum(c));
                end
                % mini = min(durs);
                % maxi = max(durs);
                % if mini < 1.8 || maxi > 3.2
                %     fprintf('    NOTE!!! Range: [%.2f %.2f]\n', mini, maxi)
                % else
                %     fprintf('    ok\n')
                % end
            end
        end
    end
end

fprintf(logFile, '\nSparse noise:\n');
for set = 1:length(setNames)
    dirs_animals = dir(fullfile(folderSave, "Data", setNames{set}));
    for animal = 1:length(dirs_animals)
        if startsWith(dirs_animals(animal).name, '.')
            continue
        end
        dirs_dates = dir(fullfile(dirs_animals(animal).folder, ...
            dirs_animals(animal).name));
        for date = 1:length(dirs_dates)
            if startsWith(dirs_dates(date).name, '.')
                continue
            end
            ft = fullfile(folderSave, "Data", setNames{set}, ...
                dirs_animals(animal).name, dirs_dates(date).name);

            % check consistency of sparse noise frame durations
            if isfile(fullfile(ft, "sparseNoise.times.npy"))
                t = readNPY(fullfile(ft, "sparseNoise.times.npy"));
                intv = readNPY(fullfile(ft, "sparseNoiseExp.intervals.npy"));
                durs = [];
                for k = 1:size(intv,1)
                    ind = t >= intv(k,1) & t <= intv(k,2);
                    durs = [durs; diff(t(ind))];
                end
                durMed = median(durs);
                c = histcounts(durs, [0 0.5*durMed 1.5*durMed max(durs)+1]);
                out = sum(c([1 3]));
                if out > 0
                    fprintf(logFile, "  %s %s\n", dirs_animals(animal).name, ...
                        dirs_dates(date).name);
                    fprintf(logFile, '    %d of %d stimuli outside tolerated time range.\n', ...
                        out, sum(c));
                end
                % if any(durs < durMed/2 | durs > 2*durMed)
                %     fprintf('    NOTE!!! Range: [%.3f %.3f]. Median: %.3f\n', ...
                %         min(durs), max(durs), durMed)
                % else
                %     fprintf('    ok\n')
                % end
            end
        end
    end
end

fprintf(logFile, '\nFull field chirps:\n');
for set = 1:length(setNames)
    dirs_animals = dir(fullfile(folderSave, "Data", setNames{set}));
    for animal = 1:length(dirs_animals)
        if startsWith(dirs_animals(animal).name, '.')
            continue
        end
        dirs_dates = dir(fullfile(dirs_animals(animal).folder, ...
            dirs_animals(animal).name));
        for date = 1:length(dirs_dates)
            if startsWith(dirs_dates(date).name, '.')
                continue
            end
            ft = fullfile(folderSave, "Data", setNames{set}, ...
                dirs_animals(animal).name, dirs_dates(date).name);

            % check length of full field times and consistency across
            % repetitions
            if isfile(fullfile(ft, "fullField.times.npy"))
                t = readNPY(fullfile(ft, "fullField.times.npy"));
                intv = readNPY(fullfile(ft, "fullFieldExp.intervals.npy"));
                % there are 13 stimulus switches per repetition
                onsets = [];
                for k = 1:size(intv,1)
                    ind = t >= intv(k,1) & t <= intv(k,2);
                    t_k = t(ind);
                    t_k = padarray(t_k, mod(13 - mod(length(t_k),13), 13), ...
                        NaN, "post");
                    t_k = reshape(t_k, 13, []);
                    onsets = [onsets, t_k];
                end
                durs = diff(onsets, 1, 1);
                dursMed = median(durs, 2, "omitnan");
                dev = (durs - dursMed) ./ dursMed;
                if any(dev(:) < -0.5 | dev(:) > 1)
                    fprintf(logFile, "  %s %s\n", dirs_animals(animal).name, ...
                        dirs_dates(date).name);
                    fprintf(logFile, '    NOTE!!! Range: [%.3f %.3f] (deviation from median)\n', ...
                        min(dev(:)), max(dev(:)));
                end
            end
        end
    end
end

fprintf(logFile, '\nCircles\n');
for set = 1:length(setNames)
    dirs_animals = dir(fullfile(folderSave, "Data", setNames{set}));
    for animal = 1:length(dirs_animals)
        if startsWith(dirs_animals(animal).name, '.')
            continue
        end
        dirs_dates = dir(fullfile(dirs_animals(animal).folder, ...
            dirs_animals(animal).name));
        for date = 1:length(dirs_dates)
            if startsWith(dirs_dates(date).name, '.')
                continue
            end
            ft = fullfile(folderSave, "Data", setNames{set}, ...
                dirs_animals(animal).name, dirs_dates(date).name);

            % check consistency of circle frame durations
            if isfile(fullfile(ft, "circles.times.npy"))
                t = readNPY(fullfile(ft, "circles.times.npy"));
                intv = readNPY(fullfile(ft, "circlesExp.intervals.npy"));
                durs = [];
                for k = 1:size(intv,1)
                    ind = t >= intv(k,1) & t <= intv(k,2);
                    durs = [durs; diff(t(ind))];
                end
                if isempty(durs)
                    fprintf(logFile, "  %s %s\n", dirs_animals(animal).name, ...
                        dirs_dates(date).name);
                    fprintf(logFile, '    Times outside protocol intervals!\n');
                    continue
                end
                durMed = median(durs);
                c = histcounts(durs, [0 0.5*durMed 1.5*durMed max(durs)+1]);
                out = sum(c([1 3]));
                if out > 0
                    fprintf(logFile, "  %s %s\n", dirs_animals(animal).name, ...
                        dirs_dates(date).name);
                    fprintf(logFile, '    %d of %d stimuli outside tolerated time range.\n', ...
                        out, sum(c));
                end
                % if any(durs < durMed/2 | durs > 2*durMed)
                %     fprintf('    NOTE!!! Range: [%.3f %.3f]. Median: %.3f\n', ...
                %         min(durs), max(durs), durMed)
                % else
                %     fprintf('    ok\n')
                % end
            end
        end
    end
end
fclose(logFile);
