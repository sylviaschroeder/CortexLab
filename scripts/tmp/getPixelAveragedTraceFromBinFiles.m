binFolder = 'J:\DATA\tempreg';
Ffolder = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\F';
metaFolder = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\InfoStructs';

batchSize = 5000;
timeWindow = 180; % in sec (for neuron imaging)
prctileNeuron = 8;

make_db_wheelTask_visualNoise;

for iSet = 2:length(db)
    fprintf('Dataset %d of %d\n', iSet, length(db))

    % load regops struct
    expStr = sprintf('%d_', db(iSet).expts);
    expStr(end) = [];
    data = load(fullfile(Ffolder, db(iSet).mouse_name, db(iSet).date, ...
        expStr, sprintf('regops_%s_%s.mat', db(iSet).mouse_name, ...
        db(iSet).date)));
    ops = data.ops1;
    
    planes = db(iSet).planesToProcess;
    
    for k = planes
        Lx = ops{k}.Lx;
        Ly = ops{k}.Ly;
        yrange = ops{k}.yrange;
        xrange = ops{k}.xrange;
        
        % get pixel averaged trace
        fid = fopen(fullfile(binFolder, sprintf('%s_%s_%s_plane%d.bin', ...
            db(iSet).mouse_name, db(iSet).date, expStr, k)), 'r');
        exp = db(iSet).noiseExp;
        trace = [];
        fr = 0;
        frames = sum(ops{k}.Nframes(1:exp-1)) + (1:ops{k}.Nframes(exp));
        while true
            mov = fread(fid, Ly*Lx*batchSize, '*int16');
            if isempty(mov)
                break
            end
            mov = reshape(mov, Ly, Lx, []);
            ind = fr + (1:size(mov,3));
            if frames(1) > ind(end)
                continue
            elseif frames(end) < ind(1)
                break
            end
            j = ismember(ind, frames);
            mov = mov(yrange, xrange, j);
            mov = reshape(mov, length(yrange)*length(xrange), []);
            mov = mean(mov, 1);
            trace = [trace; mov'];
            fr = ind(end);
        end
        fclose(fid);
        
        % translate to meta structure
        filename = fullfile(Ffolder, db(iSet).mouse_name, db(iSet).date, ...
            expStr, sprintf('F_%s_%s_plane%d.mat', db(iSet).mouse_name, ...
            db(iSet).date, k));
        ops2.makeCompatibleWith = '-v2';
        ops2.microID = db(iSet).microID;
        allExpInfos = preproc.dat2info(filename, ops2);
        meta = allExpInfos{db(iSet).noiseExp};
        
        % remove slow drift from trace
        frameTimes = ppbox.getFrameTimes(meta);
        samplingRate = 1/median(diff(frameTimes));
        F = preproc.removeSlowDrift(trace, samplingRate, timeWindow, prctileNeuron);
        meta.F_final = F;
        
        targetFolder = fullfile(metaFolder, db(iSet).mouse_name, db(iSet).date, ...
            num2str(db(iSet).expts(db(iSet).noiseExp)));
        targetFile = fullfile(targetFolder, ...
            sprintf('%s_%d_%s_2P_plane%03d_ROI.mat', db(iSet).date, ...
            db(iSet).expts(db(iSet).noiseExp), db(iSet).mouse_name, k));
        if ~isdir(targetFolder)
            mkdir(targetFolder)
        end
        save(targetFile, 'meta');
        
    end
end