% filename = '\\ZSERVER.cortexlab.net\Data\Subjects\SS048\2015-09-03\1\2015-09-03_1_SS048_2P_001_003.tif';
% missing = 81:82; % start and end frame that are missing
% copy = 75:76;
% filename = '\\ZSERVER.cortexlab.net\Data\Subjects\SS048\2015-09-26\1\2015-09-26_1_SS048_2P_001_013.tif';
% missing = 3619:3620; % start and end frame that are missing
% copy = 3613:3614;
% filename = '\\ZSERVER.cortexlab.net\Data\Subjects\SS052\2015-12-01\1\2015-12-01_1_SS052_2P_001_015.tif';
% missing = [17 18]; % start and end frame that are missing
% copy = [9 10];
filename = '\\ZSERVER.cortexlab.net\Data\Subjects\M150803_SS052\2015-12-03\1\2015-12-03_1_M150803_SS052_2P_001_005.tif';
missing = [3827 3832];
copy = [3819:3824];
% filename = 'J:\test\2015-09-03_1_SS048_2P_001_003.tif';
% missing = [81 82];
% copy = 75:76;

data = loadFramesBuff(filename);

% BETTER: copy manually to make sure the resulting file is EXACTLY the same
% [path, name, ext] = fileparts(filename);
% mkdir(fullfile(path, 'original_missedFrames'))
% TiffWriter(int16(data), fullfile(path, 'original_missedFrames', ...
%     [name, ext]), 16);

% assumes that missing frames are successive frames!!!
data = cat(3, data(:,:,1:missing(1)-1), data(:,:,copy), ...
    data(:,:,missing(1):end));
TiffWriter(data, filename, 16);
clear