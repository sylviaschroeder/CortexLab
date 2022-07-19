%% Parameters
% outputFrameRrate = 30;
% spatialAvg = 8; % in pixels
spatialAvg = 1; % in pixels
tempAvg = 1;
movieStart = -1; % relative to stimulus start (in sec)
movieDur = 40; % in sec
hasTwoROIs = true;
% cutOffs = [-8000 0];
cutOffs = [-16400 -14800];

%% Paths
% folders
% folderData = 'Z:\RawData\SS109\2022-03-24\1';
% folderData = 'C:\Data\SS109_2022-03-24_1';
folderData = 'C:\Data\SS109_2022-05-09_1';
folderTools = 'C:\dev\toolboxes';
folderMyRepos = 'C:\dev\workspaces';

% add paths
addpath(genpath(fullfile(folderMyRepos, 'CortexLab')))
addpath(genpath(fullfile(folderTools, 'npy-matlab')))

%% Load data
% load tiff of 2P movie
list = dir(fullfile(folderData, 'file_*.tif'));
neuroFrames = [];
for k = 1:length(list)
    neuroFrames = cat(3, neuroFrames, tiffreadVolume(fullfile(folderData, ...
        list(k).name))); % [rows x columns x time]
%     neuroFrames = cat(3, neuroFrames, tiffreadVolume(fullfile(folderData, list(1).name), ...
%         'PixelRegion', {[1 Inf], [1 Inf], [1 Inf]})); % [rows x columns x time]
end

% load frame times of tiff movie
neuroTime = readNPY(fullfile(folderData, '2pCalcium.timestamps.npy')); % in msec

% load stimulus frames
stimFrames = readNPY(fullfile(folderData, 'sparseNoise.mapArray.npy')); % [time x rows x columns]

% load frame times of stimulus
stimTime = readNPY(fullfile(folderData, 'sparseNoise.times.npy')); % in msec

%% Process data
ind = neuroTime >= stimTime(1) + movieStart & ...
    neuroTime <= stimTime(1) + movieStart + movieDur;
neuroTime(~ind) = [];
neuroFrames(:,:,~ind) = [];

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

% figure('Position', [680 560 1070 420])
% figure('Position', [640 110 1070 712])
figure('Position', [8 110 1702 658])
tiledlayout(2,4, 'TileSpacing','tight','Padding','tight')
nexttile([2 3]);
axis image off
colormap gray
set(gca, 'nextplot', 'replacechildren')
nexttile(4);
axis image off
colormap gray
set(gca, 'nextplot', 'replacechildren')
for f = 1:length(movieTime)
    nexttile(1)
    imagesc(movieNeuro(:,:,f), cutOffs)
    nexttile(4)
    imagesc(movieStim(:,:,f), [-1 1])
    
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v)
close gcf