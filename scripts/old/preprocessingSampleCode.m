% this script is showing how to use functions included in +ppbox package to
% do basic pre-processing of the 2p data
% it was written having single-channel multiple-plane B-Scope data in mind
% adaptation to 2-channel tiffs is in the plans

% run this script section-by-section to figure out how the functions should
% be used

% by Michael Krumin, started early 2014

% check ppbox.getTmpFolder function and configure it for you computer to make sure
% you have enough space for all the processed files there (By default, you will
% need twice the amount of space as the raw tiff files. You can cut this in
% half if you choose not to generate registered tiffs automatically).
% Do not use network drives, this will low down the process. If you have
% SSD - use it, it will make things faster. The ppbox code tries to be
% memory-efficient (this is required for large datasets), which means it is
% writing to disk more than minimally needed.

Experiments=4;
subject=repmat({'M150410_SS045'}, length(Experiments), 1);
expDate = repmat({'2015-05-04'}, length(Experiments), 1);
planes = 3;
refExp = 3; % use this exp for target frame selection for registration. 
% try to choose the best experiment, which is supposed to drive as many cells as possible


for exp = 1:length(Experiments)
    % initializing
%     refExpInd = find(Experiments == refExp);
    info=ppbox.infoPopulate(subject{exp}, expDate{exp}, Experiments(exp));
    
    % load the datasets and divide into sub-datasets with different planes
    
    % let's start parallel pool
    gcp;
    clear options;
    options.nFramesPerChunk = 4096;
    options.nFrames4TargetSelection=300;
    options.resolution = 2;
    options.autoRect = false;
    options.nY = 1;
    options.nX = 1;
    options.rectOverlap = 80; % pixels
    
    for iPlane = planes;
        options.iPlane = iPlane;
        
        % extract single plane from the tiff files
        infoRaw = ppbox.extractSinglePlane(info, options);
        % if the unregistered .bin file already exist, then instead of
        % extracting single plane just get information about that file
        % [~, ~, infoRaw] = loadArrInfo(fullfile(info.folderProcessed, sprintf('%s_plane%03d_raw', info.basename2p, iPlane)));
        
        % register one plane
        infoReg = ppbox.registerSinglePlanePiecewise(infoRaw, options);
        
        % we don't need the raw data .bin file anymore. If needed, can be
        % rebuilt again from the original tiffs
        
        %     delete(fullfile(infoRaw.folderProcessed, [infoRaw.basenameRaw, '.bin']));
    end
end
delete(gcp);
return;

%% now generate tiffs (if not done yet) and check how good the registration was
% then, delete the raw data, probably also plane 1 and continue
% use ppbox.bin2tiff to generate tiffs

%% now register other experiments using the same targetFrames
% all the comments from the previous section apply here as well

% getting the target frames from already analysed dataset

info=ppbox.infoPopulate(subject{refExpInd}, expDate{refExpInd}, refExp)

clear infoReg;

for iPlane = Planes
    [~, ~, infoReg{iPlane}] = loadArrInfo(fullfile(info.folderProcessed, sprintf('%s_plane%03d_registered', info.basename2p, iPlane)));
end

% gcp

for exp = Experiments(Experiments~=refExp)
    
    % and this is the experiment dataset to analyse now
    expInd = find(Experiments == exp);
    info=ppbox.infoPopulate(subject{expInd}, expDate{expInd}, exp);
    
%     tic
    clear options;
    options.nFramesPerChunk = 1024;
    options.quickReg = true;
    for iPlane = Planes
        options.iPlane = iPlane;
        infoRaw = ppbox.extractSinglePlane(info, options);
        options.targetFrame = infoReg{iPlane}.targetFrame;
        infoRegistered = ppbox.registerSinglePlane(infoRaw, options);
        delete(fullfile(infoRaw.folderProcessed, [infoRaw.basenameRaw, '.bin']));
%         ppbox.bin2tiff(fullfile(infoRegistered.folderProcessed, [infoRegistered.basenameRegistered, '.bin']));
    end
%     toc
    
end

delete(gcp);

toc
return;

%% now generate tiffs (if not done yet) and check how good the registration was
% then, delete the raw data, probably also plane 1 and continue
% use ppbox.bin2tiff to generate tiffs

