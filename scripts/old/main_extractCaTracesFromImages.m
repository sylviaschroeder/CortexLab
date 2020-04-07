%% Register the most important/longest/most stimulating experiment
subject='M150410_SS044';
expDate = '2015-04-28';
exp=3;
planes=1:5;

% initialize
info=ppbox.infoPopulate(subject, expDate, exp);

% load the datasets and divide into sub-datasets with different planes
% Go to local folder first!!! (to save temporary files)
gcp;
options.nFramesPerChunk = 4096;
options.nFrames4TargetSelection = 300;
options.resolution = 4;
options.autoRect = false;
options.nY = 1;
options.nX = 1;
options.rectOverlap = 50;
% options.satMin = -700;
% options.satMax = 4500;

for iPlane = planes
    options.iPlane = iPlane;
    
    % extract single plane from the tiff files
    infoRaw = ppbox.extractSinglePlane(info, options);
    % if the unregistered .bin file already exist, then instead of
    % extracting single plane just get information about that file
%     [~, ~, infoRaw] = loadArrInfo(fullfile(info.folderProcessed, ...
%         sprintf('%s_plane%03d_raw', info.basename2p, iPlane)));

    % register one plane
    infoReg = ppbox.registerSinglePlanePiecewise(infoRaw, options);
end
delete(gcp);

%% Generate tiffs and check how good the registration was
subject='M150410_SS044';
expDate = '2015-04-28';
exp=3;
planes=3;
channels = 2;

% specify rectangles to be used
rects = cell(length(planes), 1);
% rects{1}={2,[1 2]};
% rects{2}={2,1};
% rects{3}={2,[1 2]};
% rects{4}={2,1:3};
% rects{5}={2,1:3};
for iPlane=1:length(planes)
    rects{iPlane} = {1,1};
%     rects{iPlane} = {2,[1 2]};
end
% rects = [];

% initialize
info=ppbox.infoPopulate(subject, expDate, exp);

for iPlane = 1:length(planes)
    % load meta, arrSize and arrPrecision
    if isempty(rects)
        load(fullfile(info.folderProcessed, ...
            sprintf('%s_plane%03d_registered.mat', info.basename2p, planes(iPlane))));
%           load(fullfile(info.folderProcessed, ...
%               sprintf('%s_plane%03d_raw.mat', info.basename2p, planes(iPlane))));
        for iCh = 1:channels
            ppbox.bin2tiff(fullfile(meta.folderProcessed, ...
                [meta.chData(iCh).basename '_registered.bin']), ...
                arrSize, arrPrecision);
%             ppbox.bin2tiff(fullfile(meta.folderProcessed, ...
%                 [meta.chData(iCh).basename '_raw.bin']), ...
%                 arrSize, arrPrecision);
        end
    else
        for iPass = 1:size(rects{iPlane},1)
            for iRect = 1:length(rects{iPlane}{iPass, 2})
                load(fullfile(info.folderProcessed, ...
                    sprintf('%s_plane%03d_rect%02d_%02d_registered.mat',...
                    info.basename2p, planes(iPlane), rects{iPlane}{iPass,1}, ...
                    rects{iPlane}{iPass,2}(iRect))));
                for iCh = channels
                    ppbox.bin2tiff(fullfile(info.folderProcessed, ...
                        [meta.chData(iCh).basename '_registered.bin']), ...
                        arrSize, arrPrecision);
                end
            end
        end
    end
end


%% If registration was satisfactory, delete all raw .bin files
for iPlane = planes
    [~, ~, infoRaw] = loadArrInfo(fullfile(info.folderProcessed, ...
        sprintf('%s_plane%03d_raw', info.basename2p, iPlane)));
    for iCh = 1:infoRaw.nChannels
        delete(fullfile(infoRaw.folderProcessed, ...
            [infoRaw.chData(iCh).basename, '_raw.bin']));
    end
end

%% Register other experiments in same session
subject='M150323_SS042';
expDate = '2015-04-14';
refExp=2;
experiments = 1:4;
planes=1:5;

% specify rectangles to be used
rects = cell(length(planes), 1);
rects{1}={2,[1 2]};
rects{2}={2,1};
rects{3}={2,[1 2]};
rects{4}={2,1:3};
rects{5}={2,1:3};
% for iPlane=1:length(planes)
%     rects{iPlane} = {1,1};
% end
% rects{1} = {1, [1 2]; 2, [3]}; % 1st col: pass, 2nd col: rectangles

% get target frame and crop positions from reference experiment
targetFrames = cell(length(planes),1);
cropPositions = cell(length(planes),1);
regChan = [];
infoRef = ppbox.infoPopulate(subject, expDate, refExp);
for iPlane=1:length(planes)
    cropPos = [];
    for iPass=1:size(rects{iPlane},1)
        for iRect=1:length(rects{iPlane}{iPass,2})
            [~,~,infoReg]=loadArrInfo(fullfile(infoRef.folderProcessed, ...
                sprintf('%s_plane%03d_rect%02d_%02d_registered', ...
                infoRef.basename2p, planes(iPlane), rects{iPlane}{iPass,1}, ...
                rects{iPlane}{iPass,2}(iRect))));
            if isempty(regChan)
                regChan=infoReg.registrationChannel;
            end
            if isempty(targetFrames{iPlane})
                targetFrames{iPlane} = infoReg.targetFrame;
            end
            if isempty(cropPos)
                cropPos = infoReg.cropPosSource;
            else
                cropPos(end+1)=infoReg.cropPosSource;
            end
        end
    end
    cropPositions{iPlane}=cropPos;
end

experiments = setdiff(experiments, refExp);
options.nFramesPerChunk = 4096;
options.nFrames4TargetSelection = 300;
options.resolution = 2;
options.autoRect = false;
options.referenceExperiment = refExp;

gcp;
    
for iExp = 1:length(experiments)
    exp=experiments(iExp);
    info = ppbox.infoPopulate(subject, expDate, exp);
    
    for iPlane = 1:length(planes)
        options.targetFrameForAlignment = targetFrames{iPlane};
        options.iPlane = planes(iPlane);
        options.cropPosSource = cropPositions{iPlane};
        options.registrationChannel = regChan;
        options.rectIDs = rects{iPlane};
        
        % extract single plane from the tiff files
        infoRaw = ppbox.extractSinglePlane(info, options);
        % if the unregistered .bin file already exist, then instead of
        % extracting single plane just get information about that file
%         [~, ~, infoRaw] = loadArrInfo(fullfile(info.folderProcessed, ...
%             sprintf('%s_plane%03d_raw', info.basename2p, planes(iPlane))));
        
        % register one plane
        infoReg = ppbox.registerSinglePlanePiecewise(infoRaw, options);
        % delete all raw .bin files
        for iCh = 1:infoRaw.nChannels
            delete(fullfile(infoRaw.folderProcessed, ...
                [infoRaw.chData(iCh).basename, '_raw.bin']));
        end
    end
end
delete(gcp);

%% now generate tiffs and check how good the registration was
subject='M150323_SS042';
expDate = '2015-04-14';
refExp=2;
experiments = 1:4;
planes=1:5;
channels = 1;

% specify rectangles to be used
rects = cell(length(planes), 1);
rects{1}={2,[1 2]};
rects{2}={2,1};
rects{3}={2,[1 2]};
rects{4}={2,1:3};
rects{5}={2,1:2};
% for iPlane=1:length(planes)
%     rects{iPlane} = {1,1};
% end
% rects{1} = {1,[1 2]};
% rects{1} = {1, [1 2]; 2, [3]}; % 1st col: pass, 2nd col: rectangles

experiments = setdiff(experiments, refExp);
for exp = experiments
    % initialize
    info=ppbox.infoPopulate(subject, expDate, exp);
    
    for iPlane = 1:length(planes)
        % load meta, arrSize and arrPrecision
        if isempty(rects)
            load(fullfile(info.folderProcessed, ...
                sprintf('%s_plane%03d_registered.mat', info.basename2p, ...
                planes(iPlane))));
            for iCh = 1:channels
                ppbox.bin2tiff(fullfile(meta.folderProcessed, ...
                    [meta.chData(iCh).basename '_registered.bin']), ...
                    arrSize, arrPrecision);
            end
        else
            for iPass = 1:size(rects{iPlane},1)
                for iRect = 1:length(rects{iPlane}{iPass, 2})
                    load(fullfile(info.folderProcessed, ...
                        sprintf('%s_plane%03d_rect%02d_%02d_registered.mat',...
                        info.basename2p, planes(iPlane), rects{iPlane}{iPass,1}, ...
                        rects{iPlane}{iPass,2}(iRect))));
                    for iCh = channels
                        ppbox.bin2tiff(fullfile(info.folderProcessed, ...
                            [meta.chData(iCh).basename '_registered.bin']), ...
                            arrSize, arrPrecision);
                    end
                end
            end
        end
    end
end

% NOW delete all tiffs

%% Compare registered average frames (of registration channel) between planes
% just plotting
subject='M150323_SS042';
expDate = '2015-04-14';
experiments = 1:4;
planes=1:5;

% specify rectangles to be used
rects = cell(length(planes), 1);
rects{1}={2,1};
rects{2}={2,1};
rects{3}={2,[1 2]};
rects{4}={2,1:2};
rects{5}={2,1:2};
% for iPlane=1:length(planes)
%     rects{iPlane} = {1,1};
% end
% rects{1} = {1,[1 2]};
% rects{1} = {1, [1 2]; 2, [3]}; % 1st col: pass, 2nd col: rectangles

figureHandles = cell(1, length(planes));
screenSize = get(0, 'ScreenSize');
rows = floor(sqrt(length(experiments)));
cols = ceil(length(experiments)/rows);
for exp = 1:length(experiments)
    % initialize
    info=ppbox.infoPopulate(subject, expDate, experiments(exp));
    
    for iPlane = 1:length(planes)
        % load meta, arrSize and arrPrecision
        k = 1;
        for iPass = 1:size(rects{iPlane},1)
            for iRect = 1:length(rects{iPlane}{iPass, 2})
                load(fullfile(info.folderProcessed, ...
                    sprintf('%s_plane%03d_rect%02d_%02d_registered.mat',...
                    info.basename2p, planes(iPlane), rects{iPlane}{iPass,1}, ...
                    rects{iPlane}{iPass,2}(iRect))));
                regChan = meta.registrationChannel;
                
                if exp == 1
                    figureHandles{iPlane}(k) = figure('Position',[10 50 ...
                        screenSize(3)-20 screenSize(4)-140]);
                end
                figure(figureHandles{iPlane}(k))
                subplot(rows,cols,exp)
                imagesc(meta.chData(regChan).registeredAVG)
                colormap gray
                axis image
                set(gca,'XTick',[],'YTick',[])
                if exp == 1
                    title(sprintf('Plane %d, rect. %02d-%02d',planes(iPlane), ...
                        rects{iPlane}{iPass,1}, rects{iPlane}{iPass,2}(iRect)))
                end
                ylabel(sprintf('Exp. %d',experiments(exp)))
                k=k+1;
            end
        end
    end
end

%% Align rectangles across experiments
subject='M150410_SS044';
expDate = '2015-04-28';
refExp=3;
experiments = 4;
planes=3;

options.resolution = 4;

% specify rectangles to be used
rects = cell(length(planes), 1);
% rects{1}={2,1};
% rects{2}={2,1};
% rects{3}={2,[1 2]};
% rects{4}={2,1:2};
% rects{5}={2,1:2};
for iPlane=1:length(planes)
    rects{iPlane} = {1,1};
end
% rects{1} = {1,[1 2]};
% rects{1} = {1, [1 2]; 2, [3]}; % 1st col: pass, 2nd col: rectangles

% initialize (ref. exp.)
infoRef=ppbox.infoPopulate(subject, expDate, refExp);

experiments = setdiff(experiments, refExp);
for iPlane = 1:length(planes)
    for iPass = 1:size(rects{iPlane},1)
        for iRect = 1:length(rects{iPlane}{iPass, 2})
            % load info structures of registered rectangles (ref. exp.)
            data = load(fullfile(infoRef.folderProcessed, ...
                sprintf('%s_plane%03d_rect%02d_%02d_registered.mat',...
                infoRef.basename2p, planes(iPlane), rects{iPlane}{iPass,1}, ...
                rects{iPlane}{iPass,2}(iRect))));
            infoRefReg = data.meta;
            for exp = experiments
                % initialize (other exp.)
                infoOth=ppbox.infoPopulate(subject, expDate, exp);
                % load info structures of registered rectangles (other exp.)
                filePath = fullfile(infoOth.folderProcessed, ...
                    sprintf('%s_plane%03d_rect%02d_%02d_registered.mat',...
                    infoOth.basename2p, planes(iPlane), rects{iPlane}{iPass,1}, ...
                    rects{iPlane}{iPass,2}(iRect)));
                data = load(filePath);
                infoOthReg = data.meta;
                
                % align rectangle to that of reference experiment
                [dx, dy] = ppbox.alignRectangles( ...
                    infoRefReg.chData(infoRefReg.registrationChannel).registeredAVG, ...
                    infoRefReg.cropPosDest, infoRefReg.validX, infoRefReg.validY, ...
                    infoOthReg.chData(infoOthReg.registrationChannel).registeredAVG, ...
                    infoOthReg.cropPosDest, infoOthReg.validX, infoOthReg.validY, ...
                    options);
                
                if abs(dx) <= 0.5
                    dx = 0;
                end
                if abs(dy) <= 0.5
                    dy = 0;
                end
                infoOthReg.dxAlign = round(dx);
                infoOthReg.dyAlign = round(dy);
                % Added on 02.11.2015 -------------------------------------
                infoOthReg.cropPosDest.ymin = infoOthReg.cropPosDest.ymin - round(dy);
                infoOthReg.cropPosDest.ymax = infoOthReg.cropPosDest.ymax - round(dy);
                infoOthReg.cropPosDest.xmin = infoOthReg.cropPosDest.xmin - round(dx);
                infoOthReg.cropPosDest.xmax = infoOthReg.cropPosDest.xmax - round(dx);
                % ---------------------------------------------------------
                data.meta=infoOthReg;
                save(filePath, '-struct', 'data');
            end
        end
    end
end

%% Detect movement (e.g. in z-direction) in movie
% save indices of detected frames in info.movingFrames
subject='M150323_SS042';
expDate = '2015-04-14';
experiments = 1:4;
planes=1:5;

% specify rectangles to be used
rects = cell(length(planes), 1);
rects{1}={2,1};
rects{2}={2,1};
rects{3}={2,[1 2]};
rects{4}={2,1:2};
rects{5}={2,1:2};
% for iPlane=1:length(planes)
%     rects{iPlane} = {1,1};
% end
% rects{1} = {1,[1 2]};
% rects{1} = {1, [1 2]; 2, [3]}; % 1st col: pass, 2nd col: rectangles

options.upperThreshold = 1.5;
options.lowerThreshold = 0.8;
for iExp = 1:length(experiments)
    exp=experiments(iExp);
    info = ppbox.infoPopulate(subject, expDate, exp);
    
    for iPlane = 1:length(planes)
        for iPass=1:size(rects{iPlane},1)
            for iRect=1:length(rects{iPlane}{iPass,2})
                filePath = fullfile(info.folderProcessed, ...
                    sprintf('%s_plane%03d_rect%02d_%02d_registered', ...
                    info.basename2p, planes(iPlane), rects{iPlane}{iPass,1}, ...
                    rects{iPlane}{iPass,2}(iRect)));
                [~,~, infoReg]=loadArrInfo(filePath);
                movie = loadArr(fullfile(infoReg.folderProcessed, ...
                    [infoReg.chData(2).basename '_registered']));
                times = ppbox.getFrameTimes(infoReg);
                sampleRate = 1 / median(diff(times));
                detectedFrames = ppbox.detectMovingFrames(movie, sampleRate, options);
                infoReg.movingFrames = detectedFrames;
                meta = infoReg;
                save(filePath, 'meta', '-append')
            end
        end
    end
end

%% Detect ROIs and extract Ca traces in each rectangle of reference experiment
subject='M150410_SS045';
expDate = '2015-05-04';
exp=2;
iPlane = 1;
iRect='_rect01_01';
% iRect='';

info=ppbox.infoPopulate(subject, expDate, exp);
CellRadius = 5; % 4 for zoom 2 experiments, 6 for zoom 3
DontPause = 0;
% frames2load = 16000; % this number should be tweaked on every system (as many frames as possible, but while staying within the memory limits
posMain = [-280 410 2400 700];
posTrace = [-250 10 2400 380];

% for iPlane = planes
basenameRegistered = sprintf('%s_plane%03d%s_registered', ...
    info.basename2p, iPlane, iRect);
filePath = fullfile(info.folderProcessed, basenameRegistered);    
    
[sz, datatype, infoReg] = loadArrInfo(filePath);
nFrames = sz(3);
filePathMovie = fullfile(info.folderProcessed, [infoReg.chData(1).basename ...
    '_registered']);
% [data, infoReg] = loadMovieFrames(filePathMovie, 1, min(nFrames, frames2load), sz, datatype);

[ROI.CellMaps, ROI.CellClasses] = AutoROIDev(filePathMovie, [], ...
    CellRadius, [], [], DontPause, 0, posMain, posTrace, 0.6);
% [ROI.CellMaps, ROI.CellClasses] = AutoROIDev(data, [], CellRadius, [], [], DontPause, 0, posMain, posTrace);

% clear data;
infoROI = infoReg;
infoROI.ROI = ROI;
% the following rectangles define the top-left and bottom-right corners
% of the FOVs
rectSource = [infoReg.validY(1), infoReg.validX(1), infoReg.validY(end), infoReg.validX(end)];
rectDest = [infoReg.cropPosDest.ymin infoReg.cropPosDest.xmin ...
    infoReg.cropPosDest.ymax infoReg.cropPosDest.xmax];
infoROI.targetFrameROI = ppbox.remapROI(ROI.CellMaps, rectSource, rectDest);

% now apply ROI chunk-by-chunk
framesPerChunk = 8192;
nChunks = ceil(nFrames/framesPerChunk);
    
for iCh = 1:infoReg.nChannels
    chanFilePath = fullfile(info.folderProcessed, [infoReg.chData(iCh).basename ...
        '_registered']);
    for iChunk = 1:nChunks
        frameFrom = framesPerChunk*(iChunk-1)+1;
        frameTo = min(nFrames, framesPerChunk*iChunk);
        [data, ~] = loadMovieFrames(chanFilePath, frameFrom, frameTo, sz, datatype);
        if iChunk == 1;
            infoROI.chData(iCh).F = ppbox.applyROI(data, ROI.CellMaps);
        else
            infoROI.chData(iCh).F = [infoROI.chData(iCh).F; ppbox.applyROI(data, ROI.CellMaps)];
        end
        if iCh == 1 % green channel
            infoROI.F = infoROI.chData(iCh).F;
        end
        clear data;
    end
end

s.arrPrecision = datatype;
s.arrSize = sz;
s.meta = infoROI;
fileSavePath = strrep(filePath, 'registered', 'ROI');
save(fileSavePath, '-struct', 's');

% save tiff-file with cell masks
spatial.plotCellMasks(infoROI, 1, 1);
% end

%% apply the same ROIs on a different dataset 
% NOTE: these datasets MUST have the same targetFrame
% get existing ROIs
subject='M150410_SS044';
expDate = '2015-04-28';
refExp=3;
experiments = 4;
planes=3;
options.resolution = 2;
channels = 1;

% specify rectangles to be used
rects = cell(length(planes), 1);
for iPlane=1:length(planes)
    rects{iPlane} = {1,1};
end
% rects{1} = {1,[1 2]};
% rects{1} = {1, [1 2]; 2, [3]}; % 1st col: pass, 2nd col: rectangles

experiments = setdiff(experiments, refExp);

oldInfo = ppbox.infoPopulate(subject, expDate, refExp);


for iPlane = 1:length(planes)
    for iPass = 1:size(rects{iPlane},1)
        for iRect = 1:length(rects{iPlane}{iPass, 2})
            % calculate ROI masks for other data sets (with the same targetFrame)
%             oldBasename = sprintf('%s_plane%03d_rect%02d_%02d_ROI', oldInfo.basename2p, ...
%                 planes(iPlane), rects{iPlane}{iPass,1}, rects{iPlane}{iPass,2}(iRect));
            oldBasename = sprintf('%s_plane%03d_ROI', oldInfo.basename2p, ...
                planes(iPlane));
            data = load(fullfile(oldInfo.folderProcessed, oldBasename));
            oldPlaneInfo = data.meta;
            
            for exp = experiments
                newInfo = ppbox.infoPopulate(subject, expDate, exp);
                basename = sprintf('%s_plane%03d_rect%02d_%02d_registered', ...
                    newInfo.basename2p, planes(iPlane), rects{iPlane}{iPass,1}, ...
                    rects{iPlane}{iPass,2}(iRect));
                filePath = fullfile(newInfo.folderProcessed, basename);
                [sz, datatype, infoReg] = loadArrInfo(filePath);
            
                infoROI = infoReg;
                infoROI.ROI.CellClasses = oldPlaneInfo.ROI.CellClasses;
%                 % the following rectangles define the top-left and bottom-right corners
%                 % of the FOVs
%                 shiftY = 1 - min(1,oldPlaneInfo.cropPosDest.ymin);
%                 shiftX = 1 - min(1,oldPlaneInfo.cropPosDest.xmin);
%                 if isfield(oldPlaneInfo, 'dyAlign') % if the reference was not also
%                     % used as reference for registration
%                     shiftY = shiftY + oldPlaneInfo.dyAlign;
%                     shiftX = shiftX + oldPlaneInfo.dxAlign;
%                 end
%                 rectSource = [oldPlaneInfo.validY(1) + shiftY, ...
%                     oldPlaneInfo.validX(1) + shiftX, ...
%                     oldPlaneInfo.validY(end) + shiftY, ...
%                     oldPlaneInfo.validX(end) + shiftX];
%                 shiftY = 1 - min(1,infoReg.cropPosDest.ymin);
%                 shiftX = 1 - min(1,infoReg.cropPosDest.xmin);
%                 if isfield(infoReg, 'dyAlign')
%                     shiftY = shiftY + infoReg.dyAlign;
%                     shiftX = shiftX + infoReg.dxAlign;
%                 end
%                 rectDest = [infoReg.validY(1)+shiftY, infoReg.validX(1)+shiftX, ...
%                     infoReg.validY(end)+shiftY, infoReg.validX(end)+shiftX];
%                 infoROI.ROI.CellMaps = ppbox.remapROI(oldPlaneInfo.ROI.CellMaps, ...
%                     rectSource, rectDest);

                infoROI.targetFrameROI = oldPlaneInfo.targetFrameROI;
                
                % ADDED 02.11.2015 instead of commented code above --------
                % remap targetFrameROI indices, which are given relative to
                % posDest, to indices relative to the cropped movie (validX
                % and validY)
                rectSource = [infoROI.cropPosDest.ymin ...
                    infoROI.cropPosDest.xmin infoROI.cropPosDest.ymax ...
                    infoROI.cropPosDest.xmax];
                rectDest = [infoROI.validY(1) infoROI.validX(1) ...
                    infoROI.validY(end) infoROI.validX(end)];
                infoROI.ROI.CellMaps = ppbox.remapROI(infoROI.targetFrameROI, ...
                    rectSource, rectDest);
                % ---------------------------------------------------------
                
                % and then apply these ROI masks chunk-by-chunk
                nFrames = sz(3);
                framesPerChunk = 8192;
                nChunks = ceil(nFrames/framesPerChunk);
                
                avgFrame = [];
                for iCh = channels
                    
                    chanFilePath = fullfile(infoROI.folderProcessed, ...
                        [infoROI.chData(iCh).basename '_registered']);
                    for iChunk = 1:nChunks
                        frameFrom = framesPerChunk*(iChunk-1)+1;
                        frameTo = min(nFrames, framesPerChunk*iChunk);
                        [data, ~] = loadMovieFrames(chanFilePath, frameFrom, ...
                            frameTo, sz, datatype);
                        if iChunk == 1;
                            infoROI.chData(iCh).F = ppbox.applyROI(data, infoROI.ROI.CellMaps);
                            if iCh == 1
                                avgFrame = mean(data,3);
                                avgFrame = avgFrame - min(avgFrame(:));
                                avgFrame = avgFrame ./ max(avgFrame(:));
                            end
                        else
                            infoROI.chData(iCh).F = [infoROI.chData(iCh).F; ...
                                ppbox.applyROI(data, infoROI.ROI.CellMaps)];
                        end
                        clear data;
                    end
                end
                infoROI.F = infoROI.chData(1).F;
                
                s.arrPrecision = datatype;
                s.arrSize = sz;
                s.meta = infoROI;
                % if only one crop for this plane was used
                if size(rects{iPlane},1) == 1 && ...
                        length(rects{iPlane}{iPass, 2}) == 1
                    basename = sprintf('%s_plane%03d_ROI', ...
                        newInfo.basename2p, planes(iPlane));
                    filePath = fullfile(newInfo.folderProcessed, basename);
                else
                    filePath = strrep(filePath, 'registered', 'ROI');
                end
                save(filePath, '-struct', 's');
                
                % save tiff-file with cell masks
                [maskImage, colors] = spatial.plotCellMasks(infoROI, 1, 1);
                
                if ~isempty(avgFrame)
                    mergedIm = spatial.mergeImageWithMask(avgFrame, maskImage, colors);
                    figure
                    imshow(mergedIm)
                    
                    basename = sprintf('%s_plane%03d_rect%02d_%02d_avgFrame.tiff', ...
                        newInfo.basename2p, planes(iPlane), rects{iPlane}{iPass,1}, ...
                        rects{iPlane}{iPass,2}(iRect));
                    filePath = fullfile(newInfo.folderProcessed, basename);
                    imwrite(avgFrame, filePath)
                end
                    
            end
        end
    end
end

%% Regress traces from red channel out of traces from green channel
% to correct for movement artifacts in z-dimension
subject='M150410_SS044';
expDate = '2015-04-28';
experiments = 3;
planes=3;

% specify rectangles to be used
% rects = cell(length(planes), 1);
% for iPlane=1:length(planes)
%     rects{iPlane} = {1,1};
% end
rects = [];

folder = 'C:\Temp2p\Results\preprocessing\';

for iExp = 1:length(experiments)
    exp=experiments(iExp);
    info = ppbox.infoPopulate(subject, expDate, exp);
    
    for iPlane = 1:length(planes)
        if isempty(rects)
            filePath = fullfile(info.folderProcessed, ...
                sprintf('%s_plane%03d_ROI', ...
                info.basename2p, planes(iPlane)));
            [sz, datatype, infoROI]=loadArrInfo(filePath);
            times = ppbox.getFrameTimes(infoROI);
%             sampleRate = 1 / median(diff(times));
%             movieFilePath = fullfile(infoROI.folderProcessed, ...
%                 [infoROI.chData(2).basename '_registered']);
%             movie = loadMovieFrames(movieFilePath, 1, sz(3), sz, datatype);
%             [~, errorFilt, error] = ppbox.detectMovingFrames(movie, ...
%                 sampleRate);
%             errorFilt = [];
%             error = [];
            infoCorrectedF = infoROI;
%             cleanF = ssLocal.regressRedFromGreenTraces(infoROI.F, ...
%                 infoROI.chData(2).F, times);
            cleanF = ssLocal.regressRedFromGreenTraces(infoROI.F, ...
                infoROI.chData(2).F, times, 1,  fullfile(folder, subject, ...
                expDate, num2str(exp)));
            infoCorrectedF.F = cleanF;
            
            s.arrPrecision = datatype;
            s.arrSize = sz;
            s.meta = infoCorrectedF;
            filePath = strrep(filePath, 'ROI', 'correctedF');
            save(filePath, '-struct', 's');
        else
            for iPass=1:size(rects{iPlane},1)
                for iRect=1:length(rects{iPlane}{iPass,2})
                    filePath = fullfile(info.folderProcessed, ...
                        sprintf('%s_plane%03d_rect%02d_%02d_ROI', ...
                        info.basename2p, planes(iPlane), rects{iPlane}{iPass,1}, ...
                        rects{iPlane}{iPass,2}(iRect)));
                    [sz, datatype, infoROI]=loadArrInfo(filePath);
                    times = ppbox.getFrameTimes(infoROI);
                    infoCorrectedF = infoROI;
                    cleanF = ssLocal.regressRedFromGreenTraces(infoROI.F, ...
                        infoROI.chData(2).F, times);
                    infoCorrectedF.F = cleanF;
                    
                    s.arrPrecision = datatype;
                    s.arrSize = sz;
                    s.meta = infoCorrectedF;
                    filePath = strrep(filePath, 'ROI', 'correctedF');
                    save(filePath, '-struct', 's');
                end
            end
        end
    end
end

%% Merge ROIs from all rectangles in reference experiment
