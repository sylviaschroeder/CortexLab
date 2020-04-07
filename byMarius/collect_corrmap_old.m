%% parameters and user-modifiable input

% takes every skp-th frame
skp = 1; 

% replace this with a top level directory
root_dir = 'C:\Temp2p\M141002_SS026\2014-10-29\';

% replace these with the full paths to the subfolders with the data in the
% root_dir. For example, there should be tif files in the directory 
% sprintf('%s%s', root_dir, folders{1}{1})
folders = {};
folders{1} = {'1\plane002\', '2\plane002\', '3\plane002\', '4\plane002\'};
folders{2} = {'1\plane003\', '2\plane003\', '3\plane003\', '4\plane003\'};

% how many planes per each of the experiments, should be same size as "folders"
nplanes = [5 5 5 5 5 5 5 5 5 5 5 5 5]; 

% replace this with an output location for variable "mI" that includes the means and
% correlation maps for all the datasets
savepath = 'F:\CalciumImaging\temp';

% codepath = 'C:\Dropbox\Code\autoDonut';
% addpath(genpath(codepath))

% donut_diameter = 7; % 7 for cortex, 10 for SC, donut diameter at zoom 2

% zoom = [2 2 3 2];
% max_cells = ceil(1e5./(donut_diameter^2*zoom.^2)); 
% max_cells = [400 400 200 400];
%% collect mean images
lbox = [0 9];
lsig = [0 4];
box = cell(1, length(lbox));
for i = 2:length(lbox)
    xs = repmat(1:lbox(i), lbox(i), 1);
    xs = xs - mean(xs(:));
    ys = xs';
    rs = xs.^2 + ys.^2;
    box{i} = exp(-rs/(2*lsig(i)^2));    
end
Ly = zeros(length(folders), 1);
Lx = zeros(length(folders), 1);

mI = cell(length(folders), 1);
for i = 1:length(folders)    
    for j = 1:length(folders{i})
        fprintf('extract mean and corr image: area %d out of %d, recording %d out of %d \n', ...
            i, length(folders), j, length(folders{j}));    
        
        fs = dir(sprintf('%s%s*.tif', root_dir, folders{i}{j}));
        fname = sprintf('%s%s%s', root_dir, folders{i}{j}, fs(1).name);
        
        Info = imfinfo(fname);
        
        Ly(i) = Info(1).Height;
        Lx(i) = Info(1).Width;
        I = zeros(Info(1).Height, Info(1).Width, nplanes(i), length(lbox));        
        
        reverseStr = '';
        for k = 1:length(fs)
            msg = sprintf('tif file %d out of %d, elapsed time %2.2f seconds \n', k, length(fs), toc);
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
            
            fname = sprintf('%s%s%s', root_dir, folders{i}{j}, fs(k).name);
            if k==length(fs)
                Info = imfinfo(fname);            
                Info = Info(1:(nplanes(i)*floor(length(Info)/nplanes(i))));
            end
            
            indx = [];
            for ip = 1:nplanes(i)
                indx = [indx; ip:skp*nplanes(i):length(Info)];
            end
            indx = indx(:);
            
            tempI = zeros(Info(1).Height, Info(1).Width, length(indx), 'uint16');
            for t = 1:length(indx)
                tempI(:,:,t) = imread(fname, 'Index', indx(t),'Info', Info);                
            end
            
            for ip = 1:nplanes(i)
               I(:,:,ip,1) = I(:,:,ip,1) + sum(single(tempI(:,:,ip:nplanes(i):end)), 3)/1e6;
            end
            
            for ip = 1:nplanes(i)                
                tempI0 = single(tempI(:,:,ip:nplanes(i):end));
                tempI0 = tempI0 - repmat(mean(tempI0,3), 1, 1, size(tempI0,3));

                for ib = 2:length(lbox)
                    tempIbox = filter2(box{ib}, tempI0(:,:), 'same');
                    tempIbox = reshape(tempIbox, size(tempI0)) - tempI0;
                    
                    I(:,:,ip, ib) = I(:,:,ip, ib) + mean(tempI0 .* tempIbox, 3)/1e6;
                end
            end
        end
        if j==1
            mI{i} = I;
        else
            mI{i} = mI{i} + I;
        end
    end        
end
%%
cd(savepath)
save('mI', 'mI')