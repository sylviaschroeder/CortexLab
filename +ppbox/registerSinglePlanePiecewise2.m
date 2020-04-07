function infoReg = registerSinglePlanePiecewise2(info, options)

% this function is doing the registration of the multiple planes of a 2p
% dataset

% 26.03.15 -- LFR added the option to provide Crops and targetFrameCrops to register
% rectangles according to a target experiment acquired in the same
% conditions, by specifying
% options.cropPosition;
% options.targetFrameCrop;
% options.nCrops;
% options.collage;

if nargin<2 || ~isfield(options, 'targetFrame')
    options.targetFrame='auto';
end
if nargin<2 || ~isfield(options, 'nFrames4TargetSelection')
    options.nFrames4TargetSelection=100;
end
if nargin<2 || ~isfield(options, 'nFramesPerChunk')
    options.nFramesPerChunk=512;
end
if nargin<2 || ~isfield(options, 'doClipping')
    options.doClipping=true;
end
if nargin<2 || ~isfield(options, 'nParThreads')
    % define this >1 if you want to use parrallel processing
    options.nParThreads=4;
end
if nargin<2 || ~isfield(options, 'fastSave')
    % true if using fast binary saving/loading
    options.fastSave=true;
end
if nargin<2 || ~isfield(options, 'noTiff')
    % set to true if do not want to save tiffs
    options.noTiff=false;
end
if nargin<2 || ~isfield(options, 'iPlane')
    % register all planes by default
    options.iPlane=info.iPlane;
end
if nargin<2 || ~isfield(options, 'resolution')
    options.resolution = 1; % [1/pix],
    % e.g. if options.resolution = 10 the registration resolution will be
    % down to 1/10 pixels
end
if nargin<2 || ~isfield(options, 'autoRect')
    options.autoRect = false;
end
if nargin<2 || ~isfield(options, 'nY')
    options.nY = 3;
end
if nargin<2 || ~isfield(options, 'nX')
    options.nX = 3;
end
if nargin<2 || ~isfield(options, 'rectOverlap')
    options.rectOverlap = 20; % pixels
end
if nargin<2 || ~isfield(options, 'satMin')
    options.satMin = -Inf;
end
if nargin<2 || ~isfield(options, 'satMax')
    options.satMax = Inf;
end
satMin = int16(options.satMin);
satMax = int16(options.satMax);

if nargin<2 || ~isfield(options, 'registrationChannel')
    % register red channel by default (if it exists), otherwise the first
    % analysed channel
    regCh = find(strcmp({info.chData.color}, 'red'));
    tmp = 0;
    while isempty(regCh)
        tmp = tmp + 1;
        if ~isempty(info.chData(tmp).tiffFrames)
            regCh = tmp;
        end
    end
    options.registrationChannel = regCh;
end

regCh = options.registrationChannel;
if info.nChannels>1
    chString = sprintf('_channel%03d', regCh);
else
    chString = '';
end

try
    filePath = fullfile(info.folderProcessed, [info.chData(regCh).basename '_raw']);
catch
    info.chData(regCh).basename = sprintf('%s_plane%03d%s', ...
        info.basename2p, options.iPlane, chString);
    filePath = fullfile(info.folderProcessed, [info.chData(regCh).basename '_raw']);
    
    % for backward compatibility
    greenCh = strcmp({info.chData.color}, 'green');
    if info.nChannels>1
        chString = sprintf('_channel%03d', greenCh);
    else
        chString = '';
    end
    info.basenameRaw = sprintf('%s_plane%03d%s', ...
        info.basename2p, options.iPlane, chString);
end

% this is the Gaussian filter for registration frames
hGauss = fspecial('gaussian', [7 7], 1);

regOptions.nParThreads = options.nParThreads;
regOptions.resolution = options.resolution;
regOptions.spFilter = hGauss;
regOptions.nFramesPerChunk = options.nFramesPerChunk;

nFrames2Skip=[1, 0]; % number of frames to skip in the beginning/end of the stack
% the first 1-2 piezo cycles may be unstable, and the last cycle can be corrupt

fprintf('registering file %s_raw...\n', info.chData(regCh).basename);
[sz, prec, info] = loadArrInfo(fullfile(info.folderProcessed, [info.basenamePlane '_raw']));
nFrames = sz(3);
nFrames2Use = min(nFrames, options.nFrames4TargetSelection);
% select random frames
% frames2Use = randperm(nFrames-1, nFrames2Use)+1;
% or evenly spaced frames
frames2Use = round(linspace(0.5*nFrames/nFrames2Use, ...
    (nFrames2Use-0.5)*nFrames/nFrames2Use, nFrames2Use));

m = memmapfile([filePath, '.bin']);
m.Format =  {'int16', sz, 'frames'};

% cutting the required frames out
indStart=nFrames2Skip(1)+1;
indEnd=nFrames-nFrames2Skip(2);
nFrames=nFrames-sum(nFrames2Skip);
% updating the planeFrame indices after removing the undesired frames
info.planeFrames=info.planeFrames(indStart:indEnd);

% for backward compatibility
info.meanIntensity = info.meanIntensity(indStart:indEnd);

for iCh = 1:info.nChannels
    if isempty(info.chData(iCh).tiffFrames)
        continue
    end
    info.chData(iCh).meanIntensity = info.chData(iCh).meanIntensity(indStart:indEnd);
    info.chData(iCh).tiffFrames = info.chData(iCh).tiffFrames(indStart:indEnd);
end

% the following code is adopted (and then adapted) from regTranslations()
if isequal(options.targetFrame, 'auto')
    
    %     data = nan([sz(1:2), nFrames2Use]);
    %     for iFrame = 1:nFrames2Use
    %         data(:,:,iFrame) =  loadMovieFrames(filePath, frames2Use(iFrame), frames2Use(iFrame));
    %     end
    
    data = m.Data.frames(:,:,frames2Use);
    clear m;
    
    options.targetFrame = selectTargetDev(1, data, nFrames2Use, floor(nFrames2Use/10), 1, hGauss);
    %     clear data;
    
elseif isnumeric(options.targetFrame) && length(options.targetFrame) == 1
    
    data = m.Data.frames(:,:,options.targetFrame);
    clear m
    %Gaussian filter the target image
    options.targetFrame = single(imfilter(data, hGauss, 'same', 'replicate'));
    
else
    
    data = m.Data.frames(:,:,frames2Use);
    clear m;

    % we already have a target frame provided, do nothing
    
end


%%%%%%%%%%%%%% LFR added if clause 26/3/15

if isfield(options, 'cropPosition') && isfield(options, 'targetFrameCrop')
    posSource = options.cropPosition;
    targetBank = options.targetFrameCrop;
    nRects = options.nCrops;
    collage =  options.collage;
elseif isfield(options, 'cropPosSource') && ~isempty(options.cropPosSource)
    posSource = options.cropPosSource;
    nRects = length(posSource);
else
    figure
    imagesc(options.targetFrame);
    axis equal tight;
    colormap gray
    hold on;
    
    if options.autoRect
        nY = options.nY;
        nX = options.nX;
        halfOverlap = ceil(options.rectOverlap/2);
        yy = round(linspace(0, sz(1), nY+1));
        xx = round(linspace(0, sz(2), nX+1));
        iRect = 0;
        posSource = struct([]);
        for iY = 1:nY
            for iX = 1:nX
                iRect = iRect+1;
                posSource(iRect).ymin = max(1, yy(iY)+1-halfOverlap);
                posSource(iRect).ymax = min(sz(1), yy(iY+1)+halfOverlap);
                posSource(iRect).xmin = max(1, xx(iX)+1-halfOverlap);
                posSource(iRect).xmax = min(sz(2), xx(iX+1)+halfOverlap);
            end
        end
        
        
    else
        fprintf('\nDraw the subregions that will be registered separately!\n\n')
        posSource = struct([]);
        for iRect = 1:100 % 100 looks like a reasonable maximum number of
            h = imrect(gca);
            pos = wait(h);
            set(h, 'Visible', 'off');
            hold on;
            if (prod(pos(3:4))==0)
                %             if the size of the rectangle is 0 - exit the loop
                drawnow;
                break;
            else
                posSource(iRect).ymin = max(1, floor(pos(2)));
                posSource(iRect).ymax = min(sz(1), ceil(pos(2)+pos(4)));
                posSource(iRect).xmin = max(1, floor(pos(1)));
                posSource(iRect).xmax = min(sz(2), ceil(pos(1)+pos(3)));
                
                %     plot(cumsum([pos(1) pos(3)]), cumsum([pos(2) pos(4)]), '.r')
                plot([pos(1) pos(1) pos(1)+pos(3) pos(1)+pos(3) pos(1)], [pos(2) pos(2)+pos(4) pos(2)+pos(4) pos(2) pos(2)], 'r:');
                
                %             yInd = round(pos(2)):round(pos(2)+pos(4));   % 26.03.15 -- LFR commented out. Seems unused
                %             xInd = round(pos(1)):round(pos(1)+pos(3));
                drawnow;
            end
        end
    end
    
    % if the user didn't define any rectangle - use the whole frame
    if isempty(posSource)
        posSource(1).ymin = 1;
        posSource(1).ymax = sz(1);
        posSource(1).xmin = 1;
        posSource(1).xmax = sz(2);
    end
    
    nRects = length(posSource);
    
end

for iRect = 1:nRects
    p = posSource(iRect);
    plot([p.xmin, p.xmax, p.xmax, p.xmin, p.xmin], [p.ymin, p.ymin, p.ymax, p.ymax, p.ymin], 'g:')
end
drawnow;

collage = zeros(sz(1), sz(2));
targetBank = cell(nRects, 1);
for iRect = 1:nRects
    p = posSource(iRect);
    targetBank{iRect} = ...
        selectTargetDev(1, data(p.ymin:p.ymax, p.xmin:p.xmax, :), nFrames2Use, floor(nFrames2Use/10), 1, hGauss);
    regOptions.spFilter = 1;
    [dx, dy] = img.regTranslationsNew(targetBank{iRect}, options.targetFrame(p.ymin:p.ymax, p.xmin:p.xmax), regOptions);
    dx = round(dx);
    dy = round(dy);
    regOptions.spFilter = hGauss;
    posDest(iRect).ymin = p.ymin + dy;
    %         if posDest(iRect).ymin < 1
    %             toTrim = 1-posDest(iRect).ymin;
    %             posSource(iRect).ymin = posSource(iRect).ymin + toTrim;
    %             targetBank{iRect} = targetBank{iRect}(toTrim+1:end, :);
    %             posDest(iRect).ymin = 1;
    %         end
    posDest(iRect).ymax = p.ymax + dy;
    %         if posDest(iRect).ymax > sz(1)
    %             toTrim = posDest(iRect).ymax-sz(1);
    %             posSource(iRect).ymax = posSource(iRect).ymax - toTrim;
    %             targetBank{iRect} = targetBank{iRect}(1:end-toTrim, :);
    %             posDest(iRect).ymax = sz(1);
    %         end
    posDest(iRect).xmin = p.xmin + dx;
    %         if posDest(iRect).xmin < 1
    %             toTrim = 1-posDest(iRect).xmin;
    %             posSource(iRect).xmin = posSource(iRect).xmin + toTrim;
    %             targetBank{iRect} = targetBank{iRect}(:, toTrim+1:end);
    %             posDest(iRect).xmin = 1;
    %         end
    posDest(iRect).xmax = p.xmax + dx;
    %         if posDest(iRect).xmax > sz(2)
    %             toTrim = posDest(iRect).xmax-sz(2);
    %             posSource(iRect).xmax = posSource(iRect).xmax - toTrim;
    %             targetBank{iRect} = targetBank{iRect}(:, 1:end-toTrim);
    %             posDest(iRect).xmax = sz(2);
    %         end
    yy = posDest(iRect).ymin:posDest(iRect).ymax;
    yyValidIdx = find(yy>=1 & yy<=sz(1));
    xx = posDest(iRect).xmin:posDest(iRect).xmax;
    xxValidIdx = find(xx>=1 & xx<=sz(2));
    collage(yy(yyValidIdx), xx(xxValidIdx)) = ...
        targetBank{iRect}(yyValidIdx, xxValidIdx);
    
    figure
    imshowpair(options.targetFrame, collage);
    axis equal tight;
    %     colormap gray;
    drawnow;
    
end %%%%% LFR 26.3.15
% pause;

if options.nParThreads>1
    try
        tmp = gcp('nocreate');
        if isempty(tmp)
            ppl=parpool(options.nParThreads);
        end
    catch
        % for older vesions of matlab (prior to 2013b)
        tmp = matlabpool('size');
        if ~tmp
            ppl = [];
            matlabpool;
        end
    end
end


nChunks = ceil(nFrames/options.nFramesPerChunk);

fprintf('Calculating registration parameters...\n');

allDx = nan(nFrames, nRects);
allDy = nan(nFrames, nRects);
for iChunk = 1:nChunks
    
    nChars = fprintf('chunk %d/%d\n', iChunk, nChunks);
    idx = (iChunk-1)*options.nFramesPerChunk+1:...
        min(iChunk*options.nFramesPerChunk, nFrames);
    
    frameStart = indStart + (iChunk-1)*options.nFramesPerChunk;
    frameEnd = min(frameStart + options.nFramesPerChunk - 1, indEnd);
    planeData = loadMovieFrames(filePath, frameStart, frameEnd, sz, prec);
    %     planeData = m.Data.frames(:,:,frameStart:frameEnd);
    
    for iRect = 1:nRects
        yInd = posSource(iRect).ymin:posSource(iRect).ymax;
        xInd = posSource(iRect).xmin:posSource(iRect).xmax;
        planeDataCropped = planeData(yInd, xInd, :);
        options.targetFrameCropped = targetBank{iRect};
        
        options.targetFrameCropped(options.targetFrameCropped<satMin) = satMin;
        options.targetFrameCropped(options.targetFrameCropped>satMax) = satMax;
        planeDataCropped(planeDataCropped<satMin) = satMin;
        planeDataCropped(planeDataCropped>satMax) = satMax;
        
        [dx, dy] = img.regTranslationsNew(single(planeDataCropped), options.targetFrameCropped, regOptions);
        
        allDx(idx, iRect) = dx(:);
        allDy(idx, iRect) = dy(:);
    end
    fprintf(repmat('\b', 1, nChars));
end

% If requested, clip the frames to the maximum fully valid region
dxMax = max(0, ceil(max(allDx)));
dxMin = min(0, floor(min(allDx)));
dyMax = max(0, ceil(max(allDy)));
dyMin = min(0, floor(min(allDy)));
% [h, w, ~] = size(planeData);
if options.doClipping
    for iRect = 1:nRects
        w = posSource(iRect).xmax - posSource(iRect).xmin + 1;
        h = posSource(iRect).ymax - posSource(iRect).ymin + 1;
        validX{iRect} = (1 + dxMax(iRect)):(w + dxMin(iRect));
        validY{iRect} = (1 + dyMax(iRect)):(h + dyMin(iRect));
    end
else
    for iRect = 1:nRects
        w = posSource(iRect).xmax - posSource(iRect).xmin + 1;
        h = posSource(iRect).ymax - posSource(iRect).ymin + 1;
        validX{iRect} = 1:w;
        validY{iRect} = 1:h;
    end
end

% now do the translations required
fprintf('Applying registration and saving...\n');

% check how often this function was run on the same file and create
% filenames accordingly
pass = 1;
existingFiles = dir(fullfile(info.folderProcessed, ...
    [info.basenamePlane '_rect*_registered.mat']));
if ~isempty(existingFiles)
    rectNames = regexp({existingFiles.name}, 'rect\d*', 'match');
    maxPass = 1;
    for iPass = 1:length(rectNames)
        num = sscanf(rectNames{iPass}{1}, 'rect%d');
        if num > maxPass
            maxPass = num;
        end
    end
    pass = maxPass + 1;
end

fids = cell(info.nChannels, nRects);
for iCh = 1:info.nChannels
    if isempty(info.chData(iCh).tiffFrames)
        continue
    end
    for iRect = 1:nRects
        if isfield(options, 'basenameRect') && ~isempty(options.basenameRect)
            fids{iCh, iRect} = fopen(fullfile(info.folderProcessed, ...
                sprintf('%s_registered.bin', options.basenameRect)), 'w');
        else
            fids{iCh, iRect} = fopen(fullfile(info.folderProcessed, ...
                sprintf('%s_rect%02d_%02d_registered.bin', ...
                info.chData(iCh).basename, pass, iRect)), 'w');
        end
    end
end

dataPrecision = 'int16';
for iCh = 1:info.nChannels
    if isempty(info.chData(iCh).tiffFrames)
        continue
    end
    filePath = fullfile(info.folderProcessed, [info.chData(iCh).basename '_raw']);
    regMIP = cell(nRects, nChunks, info.nChannels);
    for iChunk=1:nChunks
        nChars = fprintf('chunk %d/%d\n', iChunk, nChunks);
        idx = (iChunk-1)*options.nFramesPerChunk+1:...
            min(iChunk*options.nFramesPerChunk, nFrames);
        frameStart = indStart + (iChunk-1)*options.nFramesPerChunk;
        frameEnd = min(frameStart + options.nFramesPerChunk - 1, indEnd);
        planeData = loadMovieFrames(filePath, frameStart, frameEnd, sz, prec);
        %     planeData = m.Data.frames(:,:,frameStart:frameEnd);
        
        for iRect = 1:nRects
            yInd = posSource(iRect).ymin:posSource(iRect).ymax;
            xInd = posSource(iRect).xmin:posSource(iRect).xmax;
            planeDataRect = planeData(yInd, xInd, :);
            
            [mov, ~, ~]=img.translate(single(planeDataRect), allDx(idx, iRect), allDy(idx, iRect));
            
            try
                fwrite(fids{iCh, iRect}, int16(mov(validY{iRect}, validX{iRect}, :)), dataPrecision);
            catch ex
                fclose(fids{iCh, iRect});
                rethrow(ex);
            end
            
            regMIP{iRect, iChunk, iCh} = length(idx)*double(mean(mov(validY{iRect}, validX{iRect}, :), 3));
            
        end
        fprintf(repmat('\b', 1, nChars));
    end
end

% closing the data files and also saving the additional information in the .mat files
basenamePlane = info.basenamePlane;
channelBasenames = {info.chData.basename};
for iRect = 1:nRects
    for iCh = 1:info.nChannels
        fclose(fids{iCh, iRect});
        if isfield(options, 'basenameRect') && ~isempty(options.basenameRect)
            info.chData(iCh).basename = options.basenameRect;
        else
            info.chData(iCh).basename = sprintf('%s_rect%02d_%02d', ...
                channelBasenames{iCh}, pass, iRect);
        end
    end
    info.chData(regCh).targetFrame = options.targetFrame;
    info.registrationChannel = regCh;
    if isfield(options, 'basenameRect') && ~isempty(options.basenameRect)
        info.basenameRect = options.basenameRect;
    else
        info.basenameRect = sprintf('%s_rect%02d_%02d', basenamePlane, ...
            pass, iRect);
    end
    
    % for backward compatibility
    if isfield(options, 'basenameRegistered') && ~isempty(options.basenameRegistered)
        info.basenameRegistered = options.basenameRegistered;
    else
        info.basenameRegistered = sprintf('%s_registered', info.basenameRect);
    end
    info.targetFrame = options.targetFrame;
    info.targetFrameCrop = targetBank{iRect};
    info.cropPosition = posSource(iRect);
    info.cropPosSource = posSource(iRect);
    info.cropPosDest = posDest(iRect);
    info.nCrops = nRects;
    info.iCrop = iRect;
    info.collage = collage;
    info.validX = validX{iRect} + posDest(iRect).xmin - 1;
    info.validY = validY{iRect} + posDest(iRect).ymin - 1;
    info.dx = allDx(:, iRect);
    info.dy = allDy(:, iRect);
    info.regOptions = options;
    
    for iCh = 1:info.nChannels
        info.chData(iCh).registeredAVG = ...
            int16(sum(cell2mat(permute(regMIP(iRect, :, iCh), [1, 3, 2])), 3)/nFrames);
    end
    
    
    s.arrSize = [length(info.validY), length(info.validX), length(info.dx)];
    s.arrPrecision = dataPrecision;
    s.meta = info;
    %     save(fullfile(info.folderProcessed, sprintf('%s_rect%02d_%02d_registered', ...
    %         info.basenamePlane, pass, iRect)), '-struct', 's');
    save(fullfile(info.folderProcessed, info.basenameRegistered), '-struct', 's');
    
    if iRect == 1
        infoReg = info;
    else
        infoReg(iRect) = info;
    end
end


fprintf('Finished with plane %d\n\n', options.iPlane);


% close parallel pool, if it was opened here
if options.nParThreads>1 && exist('ppl', 'var')
    if isempty(ppl)
        % for older versions of matlab
        matlabpool close;
    else
        delete(ppl);
    end
    
end


