%% parameters and user-modifiable input

% takes every skp-th frame
skp = 1; 

% replace this with a top level directory
root_dir = 'C:\Temp2p\M141003_SS027\2014-10-23\';

% replace these with the full paths to the subfolders with the data in the
% root_dir. For example, there should be tif files in the directory 
% sprintf('%s%s', root_dir, folders{1}{1})
folders = {};
% folders{1} = {'2014-12-05dftreg\1\', '2014-12-05dftreg\2\'};
% folders{2} = {'2014-12-05dftreg\4\', '2014-12-05dftreg\3\'};
% folders{3} = {'2014-12-15dftreg\2\', '2014-12-15dftreg\3\'};
% folders{4} = {'M141105_MP001dftreg\16\', 'M141105_MP001dftreg\15\'};
folders{1} = {'1\plane001\', '2\plane001\'};
folders{2} = {'1\plane002\', '2\plane002\'};
folders{3} = {'1\plane003\', '2\plane003\'};
folders{4} = {'1\plane004\', '2\plane004\'};

% get frame shifts for each area and recording
subject='M141003_SS027';
expDate = '2014-10-23';
exps=[1 2];
planes=[1 2 3 4];
for ex = 1:length(exps)
    info=ppbox.infoPopulate(subject, expDate, exps(ex));
    for iPlane = 1:length(planes)
        basename = sprintf('%s_plane%03d_registered', info.basename2p, planes(iPlane));
        filePath = fullfile(info.folderProcessed, basename);
        [~, ~, planeInfo] = loadArrInfo(filePath);
        shifts.plane(iPlane).exp(ex).xmin = planeInfo.validX(1);
        shifts.plane(iPlane).exp(ex).xmax = planeInfo.validX(end);
        shifts.plane(iPlane).exp(ex).ymin = planeInfo.validY(1);
        shifts.plane(iPlane).exp(ex).ymax = planeInfo.validY(end);
    end
end
for iPlane = 1:length(planes)
    shifts.plane(iPlane).allXmin = max([shifts.plane(iPlane).exp.xmin]);
    shifts.plane(iPlane).allXmax = min([shifts.plane(iPlane).exp.xmax]);
    shifts.plane(iPlane).allYmin = max([shifts.plane(iPlane).exp.ymin]);
    shifts.plane(iPlane).allYmax = min([shifts.plane(iPlane).exp.ymax]);
end

% how many planes per each of the experiments, should be same size as "folders"
nplanes = [1 1 1 1]; 

% replace this with an output location for variable "mI" that includes the means and
% correlation maps for all the datasets
savepath = 'C:\Temp2p\M141003_SS027\2014-10-23\result_Marius';
if exist(savepath, 'dir') == 0
    mkdir(savepath)
end

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
            i, length(folders), j, length(folders{i}));    
        
        fs = dir(sprintf('%s%s*.tiff', root_dir, folders{i}{j}));
        fname = sprintf('%s%s%s', root_dir, folders{i}{j}, fs(1).name);
        
        Info = imfinfo(fname);
        
        Ly(i) = Info(1).Height;
        Lx(i) = Info(1).Width;
        I = zeros(shifts.plane(i).allYmax-shifts.plane(i).allYmin+1, ...
            shifts.plane(i).allXmax-shifts.plane(i).allXmin+1, ...
            nplanes(i), length(lbox));        
        
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
            
            % cut frames the same size for all recordings (for each plane)
            tempI = tempI(shifts.plane(i).allYmin-shifts.plane(i).exp(j).ymin+1 : ...
                end-(shifts.plane(i).exp(j).ymax-shifts.plane(i).allYmax), :, :);
            tempI = tempI(:, shifts.plane(i).allXmin-shifts.plane(i).exp(j).xmin+1 : ...
                end-(shifts.plane(i).exp(j).xmax-shifts.plane(i).allXmax), :);
            
            for ip = 1:nplanes(i)
               I(:,:,ip,1) = I(:,:,ip,1) + sum(single(tempI(:,:,ip:nplanes(i):end)), 3)/1e6;
            end
            
            for ip = 1:nplanes(i)                
                tempI0 = single(tempI(:,:,ip:nplanes(i):end));
                tempI0 = tempI0 - repmat(mean(tempI0,3), 1, 1, size(tempI0,3));

                for ib = 2:length(lbox)
                    tempIbox = zeros(size(tempI0), 'single');
                    for pk = 1:size(tempI0,3)
                        tempIbox(:,:,pk) = filter2(box{ib}, tempI0(:,:,pk), 'same');
                    end
                    tempIbox = tempIbox - tempI0;
                    
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

cd(savepath)
save('mI', 'mI')