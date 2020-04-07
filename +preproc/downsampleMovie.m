function ops1 = downsampleMovie(ops)

Ly = ops.Ly;
Lx = ops.Lx;
ntotframes  = ceil(sum(ops.Nframes));
nimgbatch = 2000;

fid = fopen(ops.RegFile, 'r');
fid2 = fopen([ops.RegFile(1:end-4) '_small.bin'], 'w');

Ly2 = floor(Ly/2);
Lx2 = floor(Ly/2);

for cyc = 1:ceil(ntotframes/nimgbatch)
    data = fread(fid,  Ly*Lx*nimgbatch, '*int16');
    dataClass = class(data);
    if isempty(data)
       break; 
    end
    data = single(data);
    data = reshape(data, Ly, Lx, []);
    data2 = zeros(Ly2, Lx2, size(data,3), dataClass);
    for f = 1:size(data,3)
        tmp = imresize(data(:,:,f), 0.5);
        data2(:,:,f) = tmp;
    end
    
    fwrite(fid2, data2, dataClass);
end
fclose(fid);
fclose(fid2);

ops1 = ops;
ops1.yrange = ceil(ops.yrange(1)/2):floor(ops.yrange(end)/2);
ops1.xrange = ceil(ops.xrange(1)/2):floor(ops.xrange(end)/2);
ops1.mimg = imresize(ops.mimg, 0.5);
ops1.mimg1 = imresize(ops.mimg1, 0.5);