db_boutons;
alignPl = 9;
framesForMean = 500;
minCorr = 0.9;

for iExp = 1
    exps = sprintf('%d_',db(iExp).exps);
    exps(end)= [];
    j = find(db(iExp).exp==db(iExp).exps);
    d = load(fullfile('C:\DATA\F\', db(iExp).subject, db(iExp).date, exps, ...
        sprintf('F_%s_%s_plane%d_Nk200.mat', db(iExp).subject, ...
        db(iExp).date, alignPl)));
    ops = d.ops;
    nFrames = ops.Nframes(j);
    fid = fopen(fullfile('J:\DATA\registered', db(iExp).subject, ...
        db(iExp).date, exps, ...
        sprintf('plane%d.bin',alignPl)));
    fread(fid, ops.Ly*ops.Lx*sum(ops.Nframes(1:(j-1))), '*int16');
    movie = fread(fid, ops.Ly*ops.Lx*nFrames, '*int16');
    movie = reshape(movie, ops.Ly*ops.Lx, nFrames);
    fclose(fid);
    
    meanInd = ceil((nFrames-framesForMean)/2)+(1:framesForMean);
    cutOff = median(ops.CorrFrame(meanInd))-1.253*3*mad(ops.CorrFrame(meanInd));
    meanInd(ops.CorrFrame(meanInd) < cutOff) = [];
    mimg = mean(movie(:,meanInd),2);
    
    corrAcrossPlanes = NaN(nFrames,length(db(iExp).planes));
    corrAcrossPlanes(:,alignPl) = corr(mimg, double(movie));
    for p = setdiff(1:length(db(iExp).planes),alignPl)
        fprintf('Plane %d\n',p)
        d = load(fullfile('C:\DATA\F\', db(iExp).subject, db(iExp).date, exps, ...
            sprintf('F_%s_%s_plane%d_Nk200.mat', db(iExp).subject, ...
            db(iExp).date, db(iExp).planes(p))));
        ops = d.ops;
        nFr = ops.Nframes(j);
        fid = fopen(fullfile('J:\DATA\registered', db(iExp).subject, ...
            db(iExp).date, exps, ...
            sprintf('plane%d.bin',db(iExp).planes(p))));
        fread(fid, ops.Ly*ops.Lx*sum(ops.Nframes(1:(j-1))), '*int16');
        movie = fread(fid, ops.Ly*ops.Lx*nFr, '*int16');
        movie = double(reshape(movie,ops.Ly*ops.Lx,nFr));
        fclose(fid);
        
        last = min(nFrames,nFr);
        corrAcrossPlanes(1:last,p) = corr(mimg, movie(:,1:last));
    end
    
    corrNorm = bsxfun(@rdivide, corrAcrossPlanes, max(corrAcrossPlanes, [], 2));
    corrNorm(corrNorm < minCorr) = NaN;
    corrNorm = bsxfun(@minus, corrNorm, min(corrNorm,[],2));
    corrNorm = bsxfun(@rdivide, corrNorm, max(corrNorm,[],2));
end