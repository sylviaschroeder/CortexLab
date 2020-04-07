function infoReg = registerSinglePlanePiecewise(info, options)

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
if nargin<2 || ~isfield(options, 'targetFrameForAlignment')
    options.targetFrameForAlignment=[];
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
if nargin<2 || ~isfield(options, 'cropPosSource')
    options.cropPosSource = []; % pixels
end
if nargin<2 || ~isfield(options, 'rectIDs')
    options.rectIDs = [];
end
if nargin<2 || ~isfield(options, 'referenceExperiment')
    options.referenceExperiment = 0;
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
    data = m.Data.frames(:,:,frames2Use);
    clear m;
    options.targetFrame = selectTargetDev(1, data, nFrames2Use, floor(nFrames2Use/10), 1, hGauss);
    
elseif isnumeric(options.targetFrame) && length(options.targetFrame) == 1
    data = m.Data.frames(:,:,options.targetFrame);
    clear m
    %Gaussian filter the target image
    options.targetFrame = single(imfilter(data, hGauss, 'same', 'replicate'));
    
else
    data = m.Data.frames(:,:,frames2Use);
    clear m;
    
end

if isempty(options.cropPosSource) % define new positions for rectangles that will be registered
    figure;
    imagesc(options.targetFrame);
    axis equal tight;
    colormap gray
    hold on;
    
    if options.autoRect % finding positions of rectangles is predefined (number and overlap)
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
                % if the size of the rectangle is 0 - exit the loop
                drawnow;
                break;
            else
                posSource(iRect).ymin = max(1, floor(pos(2)));
                posSource(iRect).ymax = min(sz(1), ceil(pos(2)+pos(4)));
                posSource(iRect).xmin = max(1, floor(pos(1)));
                posSource(iRect).xmax = min(sz(2), ceil(pos(1)+pos(3)));
                
                plot([pos(1) pos(1) pos(1)+pos(3) pos(1)+pos(3) pos(1)], ...
                    [pos(2) pos(2)+pos(4) pos(2)+pos(4) pos(2) pos(2)], 'r:');
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
    for iRect = 1:nRects
        p = posSource(iRect);
        plot([p.xmin, p.xmax, p.xmax, p.xmin, p.xmin], ...
            [p.ymin, p.ymin, p.ymax, p.ymax, p.ymin], 'g:')
    end
else % positions of rectangles were provided

    % if a target frame from another experiment is given
    % (targetFrameForAlignment), shift the current target frame to match the
    % given one
    posSource = options.cropPosSource;
    nRects = length(posSource);
    if ~isempty(options.targetFrameForAlignment)
        regOptions.spFilter = 1; % don't filter the targetFrame
        [dx_targetAlign, dy_targetAlign] = img.regTranslationsNew(options.targetFrame, ...
            options.targetFrameForAlignment, regOptions);
        % if dx or dy are x.5, add a constant to make sure that rounding
        % will result in next larger number (not the case for -x.5)
        if mod(dx_targetAlign, 1) == 0.5
            dx_targetAlign = dx_targetAlign + 0.1;
        end
        if mod(dy_targetAlign, 1) == 0.5
            dy_targetAlign = dy_targetAlign + 0.1;
        end
        regOptions.spFilter = hGauss;
        
        figure('Position', [260 678 1640 420])
        subplot(1,3,1)
        imshowpair(options.targetFrameForAlignment( ...
            max(1, round(1+dy_targetAlign)) : min(sz(1), round(sz(1)+dy_targetAlign)), ...
            max(1, round(1+dx_targetAlign)) : min(sz(2), round(sz(2)+dx_targetAlign))), ...
            options.targetFrame( ...
            max(1, round(1-dy_targetAlign)) : min(sz(1), round(sz(1)-dy_targetAlign)), ...
            max(1, round(1-dx_targetAlign)) : min(sz(2), round(sz(2)-dx_targetAlign))))
        title(sprintf('Target frames of plane %d: reference (green) vs. exp. %d (red)', ...
            options.iPlane, info.exp))
        subplot(1,3,2)
        imagesc(options.targetFrameForAlignment)
        colormap gray
        axis image
        hold on
        for iRect = 1:nRects
            p = posSource(iRect);
            plot([p.xmin, p.xmax, p.xmax, p.xmin, p.xmin], ...
                [p.ymin, p.ymin, p.ymax, p.ymax, p.ymin], 'g:')
        end
        title('Reference')
        subplot(1,3,3)
        imagesc(options.targetFrame)
        colormap gray
        axis image
        title(sprintf('Exp. %d', info.exp))
        drawnow
    else
        display(['WARNING! Crop positions for rectangles are given, ' ...
            'but target frame (options.targetFrameForAlignment) was not provided!']);
        dx_targetAlign = 0;
        dy_targetAlign = 0;
        
        figure
        imagesc(options.targetFrame)
        axis equal tight;
        colormap gray
        hold on;
        for iRect = 1:nRects
            p = posSource(iRect);
            plot([p.xmin, p.xmax, p.xmax, p.xmin, p.xmin], ...
                [p.ymin, p.ymin, p.ymax, p.ymax, p.ymin], 'g:')
        end
    end
    
    for iRect = 1:length(posSource)
        posSource(iRect).ymin = round(posSource(iRect).ymin - dy_targetAlign);
        posSource(iRect).ymax = round(posSource(iRect).ymax - dy_targetAlign);
        posSource(iRect).xmin = round(posSource(iRect).xmin - dx_targetAlign);
        posSource(iRect).xmax = round(posSource(iRect).xmax - dx_targetAlign);
    end
    
end
    
collage = zeros(sz(1), sz(2));
targetBank = cell(nRects, 1);
targetBankTrimmed = cell(nRects, 1);
posDest = struct('ymin', cell(nRects,1), 'ymax', cell(nRects,1), ...
    'xmin', cell(nRects,1), 'xmax', cell(nRects,1));
pDest = struct('ymin', cell(nRects,1), 'ymax', cell(nRects,1), ...
    'xmin', cell(nRects,1), 'xmax', cell(nRects,1));
tmp = mat2cell(zeros(nRects,1), ones(1,nRects));
toTrim = struct('ymin', tmp, 'ymax', tmp, 'xmin', tmp, 'xmax', tmp);
for iRect = 1:nRects
    p = posSource(iRect);
    p.ymin = max(1, p.ymin);
    p.ymax = min(sz(1), p.ymax);
    p.xmin = max(1, p.xmin);
    p.xmax = min(sz(2), p.xmax);
    pDest(iRect) = p;
    targetBank{iRect} = selectTargetDev(1, data(p.ymin:p.ymax, p.xmin:p.xmax, :), ...
        nFrames2Use, floor(nFrames2Use/10), 1, hGauss);
    targetBankTrimmed{iRect} = targetBank{iRect};
    regOptions.spFilter = 1;
    [dx, dy] = img.regTranslationsNew(targetBank{iRect}, ...
        options.targetFrame(p.ymin:p.ymax, p.xmin:p.xmax), regOptions);
    dx = round(dx);
    dy = round(dy);
    pDest(iRect).ymin = max(1, pDest(iRect).ymin+dy);
    pDest(iRect).ymax = min(sz(1), pDest(iRect).ymax+dy);
    pDest(iRect).xmin = max(1, pDest(iRect).xmin+dx);
    pDest(iRect).xmax = min(sz(2), pDest(iRect).xmax+dx);
    regOptions.spFilter = hGauss;
    posDest(iRect).ymin = posSource(iRect).ymin + dy;
    if posDest(iRect).ymin < 1
        toTrim(iRect).ymin = 1-posDest(iRect).ymin;
%         posSource(iRect).ymin = posSource(iRect).ymin + toTrim;
%         targetBank{iRect} = targetBank{iRect}(toTrim(iRect).ymin+1:end, :);
%         posDest(iRect).ymin = 1;
        if dy < 0
            targetBankTrimmed{iRect} = targetBankTrimmed{iRect}(1-dy:end,:);  %### was missing, LFR 30.06.15
        end
    end
    posDest(iRect).ymax = posSource(iRect).ymax + dy;
    if posDest(iRect).ymax > sz(1)
        toTrim(iRect).ymax = posDest(iRect).ymax-sz(1);
%         posSource(iRect).ymax = posSource(iRect).ymax - toTrim;
        if dy > 0
            targetBankTrimmed{iRect} = targetBankTrimmed{iRect}(1:end-dy, :);
        end
%         posDest(iRect).ymax = sz(1);
    end
    posDest(iRect).xmin = posSource(iRect).xmin + dx;
    if posDest(iRect).xmin < 1
        toTrim(iRect).xmin = 1-posDest(iRect).xmin;
%         posSource(iRect).xmin = posSource(iRect).xmin + toTrim;
        if dx < 0
            targetBankTrimmed{iRect} = targetBankTrimmed{iRect}(:, 1-dx:end);
        end
%         posDest(iRect).xmin = 1;
    end
    posDest(iRect).xmax = posSource(iRect).xmax + dx;
    if posDest(iRect).xmax > sz(2)
        toTrim(iRect).xmax = posDest(iRect).xmax-sz(2);
%         posSource(iRect).xmax = posSource(iRect).xmax - toTrim;
        if dx > 0
            targetBankTrimmed{iRect} = targetBankTrimmed{iRect}(:, 1:end-dx);
        end
%         posDest(iRect).xmax = sz(2);
    end
    collage(pDest(iRect).ymin : pDest(iRect).ymax, ...
        pDest(iRect).xmin : pDest(iRect).xmax) = ...
        targetBankTrimmed{iRect};
end

figure
imshowpair(options.targetFrame, collage);
axis equal tight;
title('Collage of rectangle target frames')
drawnow;
 
if options.nParThreads>1
    try
        tmp = gcp('nocreate');
        if isempty(tmp)
            ppl=parpool(options.nParThreads);
        end
    catch
        % for older versions of matlab (prior to 2013b)
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
    
    for iRect = 1:nRects
        yInd = pDest(iRect).ymin : pDest(iRect).ymax;
        xInd = pDest(iRect).xmin : pDest(iRect).xmax;
%         yInd = posSource(iRect).ymin+toTrim(iRect).ymin : ...
%             posSource(iRect).ymax-toTrim(iRect).ymax;
%         xInd = posSource(iRect).xmin+toTrim(iRect).xmin : ...
%             posSource(iRect).xmax-toTrim(iRect).xmax;
        planeDataCropped = planeData(yInd, xInd, :);
        targetFrameCropped = targetBankTrimmed{iRect};
        targetFrameCropped(targetFrameCropped < satMin) = satMin;
        targetFrameCropped(targetFrameCropped > satMax) = satMax;
        planeDataCropped(planeDataCropped < satMin) = satMin;
        planeDataCropped(planeDataCropped > satMax) = satMax;
        
        [dx, dy] = img.regTranslationsNew(single(planeDataCropped), ...
            targetFrameCropped, regOptions);
        
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
validX = cell(nRects,1);
validY = cell(nRects,1);
if options.doClipping
    for iRect = 1:nRects
%         [h,w] = size(targetBank{iRect});
%         w = posSource(iRect).xmax - posSource(iRect).xmin + 1;
%         h = posSource(iRect).ymax - posSource(iRect).ymin + 1;
        validX{iRect} = (pDest(iRect).xmin + dxMax(iRect)): ...
            (pDest(iRect).xmax + dxMin(iRect));
        validY{iRect} = (pDest(iRect).ymin + dyMax(iRect)): ...
            (pDest(iRect).ymax + dyMin(iRect));
%         validX{iRect} = (1 + dxMax(iRect)):(w + dxMin(iRect));
%         validY{iRect} = (1 + dyMax(iRect)):(h + dyMin(iRect));
    end
else
    for iRect = 1:nRects
%         [h,w] = size(targetBank{iRect});
%         w = posSource(iRect).xmax - posSource(iRect).xmin + 1;
%         h = posSource(iRect).ymax - posSource(iRect).ymin + 1;
        validX{iRect} = pDest(iRect).xmin:pDest(iRect).xmax;
        validY{iRect} = pDest(iRect).ymin:pDest(iRect).ymax;
%         validX{iRect} = 1:w;
%         validY{iRect} = 1:h;
    end
end

% now do the translations required
fprintf('Applying registration and saving...\n');

if isempty(options.rectIDs)
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
else
    pass = options.rectIDs{1};
end

fids = cell(info.nChannels, nRects);
rectIDs = 1:nRects;
if ~isempty(options.rectIDs) && length(options.rectIDs{2}) == nRects
    rectIDs = options.rectIDs{2};
end
for iCh = 1:info.nChannels
    if isempty(info.chData(iCh).tiffFrames)
        continue
    end
    for iRect = 1:length(rectIDs)
        fids{iCh, iRect} = fopen(fullfile(info.folderProcessed, ...
            sprintf('%s_rect%02d_%02d_registered.bin', ...
            info.chData(iCh).basename, pass, rectIDs(iRect))), 'w');
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
        
        for iRect = 1:nRects
            yInd = pDest(iRect).ymin : pDest(iRect).ymax;
            xInd = pDest(iRect).xmin : pDest(iRect).xmax;
%             yInd = posDest(iRect).ymin+toTrim(iRect).ymin : ...
%                 posDest(iRect).ymax-toTrim(iRect).ymax;
%             xInd = posDest(iRect).xmin+toTrim(iRect).xmin : ...
%                 posDest(iRect).xmax-toTrim(iRect).xmax;
%             yInd = posSource(iRect).ymin+toTrim(iRect).ymin : ...
%                 posSource(iRect).ymax-toTrim(iRect).ymax;
%             xInd = posSource(iRect).xmin+toTrim(iRect).xmin : ...
%                 posSource(iRect).xmax-toTrim(iRect).xmax;
            planeDataRect = planeData(yInd, xInd, :);
            
            [mov, ~, ~]=img.translate(single(planeDataRect), allDx(idx, iRect), allDy(idx, iRect));

            try
                fwrite(fids{iCh, iRect}, int16(mov(1+dyMax(iRect):end+dyMin(iRect), ...
                    1+dxMax(iRect):end+dxMin(iRect), :)), dataPrecision);
%                 fwrite(fids{iCh, iRect}, int16(mov(validY{iRect}, validX{iRect}, :)), dataPrecision);
            catch ex
                fclose(fids{iCh, iRect});
                rethrow(ex);
            end
            
            regMIP{iRect, iChunk, iCh} = length(idx)* ...
                double(mean(mov(1+dyMax(iRect):end+dyMin(iRect), ...
                1+dxMax(iRect):end+dxMin(iRect), :), 3));
            
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
        info.chData(iCh).basename = sprintf('%s_rect%02d_%02d', ...
            channelBasenames{iCh}, pass, rectIDs(iRect));
    end
    info.chData(regCh).targetFrame = options.targetFrame;
    info.registrationChannel = regCh;
    info.basenameRect = sprintf('%s_rect%02d_%02d', basenamePlane, ...
        pass, rectIDs(iRect));
    info.referenceExperiment = options.referenceExperiment;
    
    % for backward compatibility
    info.basenameRegistered = sprintf('%s_registered', info.basenameRect);
    info.targetFrame = options.targetFrame;
    
    info.targetFrameCrop = targetBank{iRect};
    info.cropPosSource = posSource(iRect);
    info.cropPosDest = posDest(iRect);
    info.nCrops = nRects;
    info.iCrop = iRect;
    info.collage = collage;
    info.validX = validX{iRect} + posDest(iRect).xmin - 1;
    info.validY = validY{iRect} + posDest(iRect).ymin - 1;
%     info.validX = validX{iRect} + posDest(iRect).xmin - 1;
%     info.validY = validY{iRect} + posDest(iRect).ymin - 1;
    info.dx = allDx(:, iRect);
    info.dy = allDy(:, iRect);
    
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