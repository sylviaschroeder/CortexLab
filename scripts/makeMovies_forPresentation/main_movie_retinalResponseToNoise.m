%% Parameters
% outputFrameRrate = 30;
% spatialAvg = 8; % in pixels
spatialAvg = 1; % in pixels
tempAvg = 1;
movieStart = -1; % relative to stimulus start (in sec)
movieDur = 40; % in sec
hasTwoROIs = false;
% cutOffs = [-8000 -3000];
cutOffs = [-16000 -14000];

% if reading binary movie
isBinary = false;
noRows = 1024;
noFrames = 15 * 46;

%% Paths
% folders
% folderData = 'Z:\RawData\SS113\2022-06-30\1';
folderData = 'C:\Data\SS113_2022-06-30_1';
% folderData = 'C:\Data\Hedes_2022-03-01_1';
folderTools = 'C:\dev\toolboxes';
folderMyRepos = 'C:\dev\workspaces';

% add paths
addpath(genpath(fullfile(folderMyRepos, 'CortexLab')))
addpath(genpath(fullfile(folderTools, 'npy-matlab')))

%% Load data
% load 2P movie
if isBinary
    fileName = fullfile(folderData, 'data.bin');
    fid = fopen(fileName);
    neuroFrames = fread(fid,  noRows * noRows * noFrames, 'int16');
    fclose(fid);
    neuroFrames = reshape(neuroFrames, noRows, noRows, []);
else
    list = dir(fullfile(folderData, 'file_*.tif'));
    neuroFrames = [];
    for k = 1 :length(list)
        neuroFrames = cat(3, neuroFrames, tiffreadVolume(fullfile(folderData, ...
            list(k).name))); % [rows x columns x time]
    end
end

% load frame times of tiff movie
% neuroTime = readNPY(fullfile(folderData, '2pCalcium.timestamps.npy')); % in msec
neuroTime = readNPY(fullfile(folderData, 'frameTimes.npy')); % in sec

% load stimulus frames
% stimFrames = readNPY(fullfile(folderData, 'sparseNoise.mapArray.npy')); % [time x rows x columns]
stimFrames = readNPY(fullfile(folderData, 'sparse.npy')); % [time x rows x columns]

% load frame times of stimulus
% stimTime = readNPY(fullfile(folderData, 'sparseNoise.times.npy')); % in msec
stimTime = readNPY(fullfile(folderData, 'stimTimes.npy')); % in sec

%% Process data
if size(neuroFrames,3) < length(neuroTime)
    neuroTime = neuroTime(1:size(neuroFrames,3));
elseif size(neuroFrames,3) > length(neuroTime)
    neuroFrames(:,:,length(neuroTime)+1 : end) = [];
end
ind = neuroTime >= stimTime(1) + movieStart & ...
    neuroTime <= stimTime(1) + movieStart + movieDur;
neuroTime(~ind) = [];
neuroFrames(:,:,~ind) = [];
if isBinary
    neuroFrames = permute(neuroFrames, [2 1 3]);
end

if hasTwoROIs
    neuroFrames = cat(2, neuroFrames(1 : size(neuroFrames,1)/2,:,:), ...
        neuroFrames(size(neuroFrames,1)/2+1 : end,:,:));
end
% average tiff frames locally and temporally
if spatialAvg > 1 || tempAvg > 1
    movieNeuro = imresize3(neuroFrames, round([ ...
        size(neuroFrames,[1 2])./spatialAvg, size(neuroFrames,3)/tempAvg]));
else
    movieNeuro = neuroFrames;
end

% resample tiff frames to match output rate
if tempAvg > 1
    newRate = size(neuroFrames,3) / size(movieNeuro,3);
    movieTime = interp1(double(neuroTime(1:size(neuroFrames,3))), ...
        (1+newRate)/2 + (0:newRate:size(movieNeuro,3)))';
    ind = movieTime >= stimTime(1) + movieStart & ...
        movieTime <= stimTime(1) + movieStart + movieDur;
    movieTime = movieTime(ind);
    movieNeuro = movieNeuro(:,:,ind);
else
    movieTime = neuroTime;
end

% resample stimulus to match output rate
stimFrames = permute(stimFrames, [2 3 1]); %[rows x columns x time]
stimFrames(stimFrames(:) == 0) = 1;
stimFrames(stimFrames(:) == -128) = 0;
[~,~,stimBins] = histcounts(movieTime, stimTime);
movieStim = zeros([size(stimFrames, [1 2]), length(movieTime)], 'int8');
ind1 = find(stimBins > 0, 1);
movieStim(:,:,ind1:end) = stimFrames(:,:,stimBins(ind1:end));

figure
histogram(movieNeuro(:))
title('Choose cutOff accordingly')

%% Make movie
v = VideoWriter(fullfile(folderData, 'movie.avi'));
v.FrameRate = round(1/median(diff(movieTime)));
open(v);

figure('Position', [640 110 1148 712])
% figure('Position', [640 110 1070 712])
tiledlayout(3,4, 'TileSpacing','tight','Padding','tight')
nexttile([3 3]);
axis image off ij
colormap gray
set(gca, 'nextplot', 'replacechildren')
nexttile(4);
axis image off
colormap gray
set(gca, 'nextplot', 'replacechildren')
for f = 1:length(movieTime)
    nexttile(1)
    imagesc(movieNeuro(:,:,f), cutOffs)
%     if f == 1
%         set(gca, 'nextplot', 'replacechildren')
%     end
    nexttile(4)
    imagesc(movieStim(:,:,f), [-1 1])
%     if f == 1
%         set(gca, 'nextplot', 'replacechildren')
%     end
    
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v)
close gcf