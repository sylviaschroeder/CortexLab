subject = 'SS047';
date = '2015-11-12';
exp = 1;
planes = 2:4;
expStr = sprintf('%d_', exp);
expStr(end) = [];
% expStr = '1_2_2';

for iPlane = planes
    fid = fopen(fullfile('J:\DATA\tempreg', sprintf('%s_%s_%s_plane%d_red.bin', ...
        subject, date, expStr, iPlane)),'r');
    load(fullfile('C:\STORAGE\OneDrive - University College London\Lab\DATA\F\', ...
        subject, date, expStr, sprintf('regops_%s_%s.mat', subject, date)));
    mov = fread(fid,ops1{iPlane}.Ly * ops1{iPlane}.Lx * sum(ops1{iPlane}.Nframes),'*int16');
    m = reshape(mov, ops1{iPlane}.Ly, ops1{iPlane}.Lx, sum(ops1{iPlane}.Nframes));
    img = mean(m, 3);
    figure,imagesc(img)
    axis image
    for iExp = exp
        load(fullfile('C:\STORAGE\OneDrive - University College London\Lab\DATA\InfoStructs', ...
            subject, date, num2str(iExp), sprintf('%s_%d_%s_2P_plane%03d_ROI.mat', ...
            date, iExp, subject, iPlane)))
        figure('Position', [1254 687 560 420]), imagesc(meta.targetFrame)
        axis image
        meta.targetFrame=img;
        save(fullfile('C:\STORAGE\OneDrive - University College London\Lab\DATA\InfoStructs', ...
            subject, date, num2str(iExp), sprintf('%s_%d_%s_2P_plane%03d_ROI.mat', ...
            date, iExp, subject, iPlane)), 'meta')
    end
end