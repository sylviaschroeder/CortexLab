subject='M141007_SS029';
expDate = '2014-11-12';
exp=4;
info=ppbox.infoPopulate(subject, expDate, exp);
planes = 2:3;

for iPlane = planes

    basenameRegistered = sprintf('%s_plane%03d_registered', info.basename2p, iPlane);
    filePath = fullfile(info.folderProcessed, basenameRegistered);
    data = loadArr(filePath);
    
    averageFrame = mean(data, 3);
    averageFrame = averageFrame - min(averageFrame(:));
    averageFrame = averageFrame ./ max(averageFrame(:));
    
    clear data
    
    save([filePath '_average.mat'], 'averageFrame');
    imwrite(averageFrame, [filePath '_average.tif'])
end