folderEye = '\\ZSERVER.cortexlab.net\Data\EyeCamera';
folderSpikes = 'C:\DATA\electrophys\';
folderResults = 'J:\Eye\';

db_ephys_driftingGratings

for k = 2:length(db)
    fprintf('Dataset %d of %d\n', k, length(db))
    vr = VideoReader(fullfile(folderEye, db(k).subject, db(k).date, ...
        num2str(db(k).eyeFolder), 'eye.mj2'));
    % load eye movie alignment with Timeline (tVid)
    load(fullfile(folderEye, db(k).subject, db(k).date, ...
        num2str(db(k).eyeFolder), 'eye_timeStamps.mat'))
    % load time vector of spike data (time)
    load(fullfile(folderSpikes, db(k).subject, db(k).date, ...
        sprintf('%02d_data.mat', db(k).exp)), 'time', 'timelineToEphys')
    
    eyeTime = tVid .* timelineToEphys(1)+timelineToEphys(2);
    start = find(eyeTime>time(1),1) - ceil(vr.FrameRate);
    stop = find(eyeTime>time(end),1) + ceil(vr.FrameRate);
    
    movie = read(vr, [start stop]);
    
    fldr = fullfile(folderResults, db(k).subject, db(k).date);
    if exist(fldr, 'dir') == 0
        mkdir(fldr);
    end
    vw = VideoWriter(fullfile(folderResults, db(k).subject, db(k).date, ...
        sprintf('%02d_eye.mj2', db(k).exp)),'Motion JPEG 2000');
    vw.MJ2BitDepth = 8;
    vw.LosslessCompression = true;
    vw.FrameRate = vr.FrameRate;
    
    open(vw);
    for f=1:size(movie,4)
        writeVideo(vw,movie(:,:,:,f));
    end
    close(vw)
    
    eyeTime = eyeTime(start:stop);
    save(fullfile(folderResults, db(k).subject, db(k).date, ...
        sprintf('%02d_eyeTime.mat', db(k).exp)), 'eyeTime')
    
    clear movie vr vw time tVid eyeTime
end

% SS061: add frames 1-8631