function [dx, dy, target, error] = regTranslationsNew(movie, target, regOptions)
% registration Movie frame registration to a target frame using
% parallelisation and array vectorisation for efficiency
%
%   [DX, DY, TARGET] = rapidReg(MOVIE, TARGET,...) find registration of
%   MOVIE (an (Y,X,T) array) to TARGET (either the target image frame (Y,X)
%   or 'auto', to find one automatically). Returns the translations required
%   to register each frame, DX and DY, and TAGRET, the target frame.
%   Optionally takes, 'noparallel', meaning use single-threaded codepath
%   instead of parallel.

% 2013-07 CB created (heavily plagiarised from Mario Dipoppa's code)
% 2014-02 MK 'FramesPerChunk' parameter supported (to do the analysis
%         chunk-by-chunk). Supply a string, which contains the value
%         immediately afterwards:
%         [DX, DY, TARGET] = regTranslations(MOVIE, TARGET,..., 'FramesPerChunk', '512')
% 2014-03 MK target can now be also defined by defining the frame number from the dataset


[h, w, nFrames] = size(movie);
temporaryPool = false;

% fprintf('Registering..');

%% Setup
%create a Gaussian filter for filtering registration frames
if nargin<2 || ~isfield(regOptions, 'spFilter')
    regOptions.spFilter = fspecial('gaussian', [7 7], 1);
end

if nargin<2 || ~isfield(regOptions, 'nParThreads')
    regOptions.nParThreads = 4;
    parallel = true;
else
    if regOptions.nParThreads == 1
        parallel = false;
    else
        parallel = true;
    end
end

if nargin<2 || ~isfield(regOptions, 'nFramesPerChunk')
    regOptions.nFramesPerChunk = nFrames;
end

if nargin<2 || ~isfield(regOptions, 'resolution')
    regOptions.resolution = 1;
end

%% If requested, compute best target frame
if strcmpi(target, 'auto')
    
    fprintf('finding target..');
    %first compute a smoothed mean of each frame
    meanF = smooth(mean(reshape(movie, h*w, nFrames)));
    %now look in the middle third of the image frames for the minimum
    fromFrame = round(nFrames*1/3);
    toFrame = round(nFrames*2/3);
    [~, idx] = min(meanF(fromFrame:toFrame));
    minFrame = fromFrame + idx;
    %Gaussian filter the target image
    target = imfilter(movie(:,:,minFrame), regOptions.spFilter, 'same', 'replicate');
elseif isnumeric(target) && length(target)==1
    % if a specific frame is defined to be a target frame
    fprintf('frame %d will be the target frame...\n', target);
    target = imfilter(movie(:,:,target), regOptions.spFilter, 'same', 'replicate');
end

%% Fourier transform the filtered movie frames for registration
% fprintf('filtering..\n');
% if ~isnan(nFramesPerChunk)
%     nChunks=ceil(nFrames/nFramesPerChunk);
%     for iChunk=1:nChunks
%         fprintf('Chunk %d/%d\n', iChunk, nChunks);
%         idx=nFramesPerChunk*(iChunk-1)+1:min(nFramesPerChunk*iChunk, nFrames);
%         movie(:,:,idx) = fft2(imfilter(movie(:,:,idx), hGauss, 'same', 'replicate'));
%     end
% else
%     movie = fft2(imfilter(movie, hGauss, 'same', 'replicate'));
% end

%% Compute required displacement and register each frame
dx = zeros(1, nFrames);
dy = zeros(1, nFrames);
error = zeros(1, nFrames);
% fprintf('registering..\n');

if parallel
    %% Register in parallel
    poolobj = gcp('nocreate');
    if isempty(poolobj) || poolobj.NumWorkers == 1
        poolobj = parpool(regOptions.nParThreads);
        temporaryPool = true;
    end
    try
        %do parallel loops in chunks of data to prevent matlab choking
        chunkSize = regOptions.nFramesPerChunk; %frames
        nChunks = ceil(nFrames/chunkSize);
%         progressStepSize = 100;
        %     pctRunOnAll('mattoolsJava');
        %     ppm = ParforProgMon('Registration: ', nFrames, progressStepSize, 300, 80);
        fftTarget = fft2(target);
        spFilter = regOptions.spFilter;
        subRes = regOptions.resolution;
        for i = 0:(nChunks - 1)
            if nChunks>1
            fprintf('Chunk %d/%d\n', i+1, nChunks);
            end
            sidx = i*chunkSize + 1;
            eidx = min((i + 1)*chunkSize, nFrames);
            n = eidx - sidx + 1;
            parfor t = sidx:eidx
                %find the best registration translation
                fftFrame = fft2(imfilter(movie(:,:,t), spFilter, 'same', 'replicate'));
                output = dftregistrationDev(fftTarget, fftFrame, subRes);
                error(t) = output(1);
                if length(output) >= 4
                    dx(t) = output(4);
                    dy(t) = output(3);
                end
                
                %         if mod(t, progressStepSize) == 0
                %           ppm.increment();
                %         end
            end
        end
        %     ppm.delete();
        if temporaryPool
            delete(poolobj); %close worker pool
        end
    catch ex
        if temporaryPool
            %in case of error, ensure temporary worker pool is closed
            delete(poolobj);
        end
        rethrow(ex)
    end
else
    %% Register sequentially
    nChars = 0;
    targetFFT = fft2(target);
    for t = 1:nFrames
        %find the best registration translation
        fftFrame = fft2(imfilter(movie(:,:,t), regOptions.spFilter, 'same', 'replicate'));
        output = dftregistration(targetFFT, fftFrame, regOptions.resolution);
        fprintf(repmat('\b', 1, nChars));
        nChars = fprintf('frame %d/%d', t, nFrames);
        error(t) = output(1);
        if length(output) >= 4
            dx(t) = output(4);
            dy(t) = output(3);
        end
    end
end

% fprintf('.done\n');

end

