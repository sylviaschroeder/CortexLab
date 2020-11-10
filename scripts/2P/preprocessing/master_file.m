%%
cd('C:\STORAGE\Dropbox\Code\scripts\preprocessing') % start this code in the directory with make_db
make_db;

ops0.toolbox_path = 'C:\STORAGE\workspaces\Suite2P';
if exist(ops0.toolbox_path, 'dir')
	addpath(genpath(ops0.toolbox_path)) % add local path to the toolbox
else
	error('toolbox_path does not exist, please change toolbox_path');
end

ops0.useGPU                 = 1; % if you can use a GPU in matlab this accelerate registration approx 3 times
ops0.fig                    = 1; % turn off figure generation with 0

% root paths for files and temporary storage (ideally an SSD drive. my SSD is C:/)
% ops0.RootStorage            = '\\ZSERVER.cortexlab.net\Data\Subjects';
% ops0.RootStorage            = '//zubjects.cortexlab.net/Subjects';
ops0.RootStorage            = 'J:\DATA\raw\';
ops0.temp_tiff              = 'J:/DATA/temp.tiff'; % copy data locally first: put this on your SSD
ops0.RegFileRoot            = 'J:/DATA/tempreg'; 
ops0.DeleteBin              = 0; % set to 1 for batch processing on a limited hard drive
ops0.ResultsSavePath        = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\F';
ops0.RegFileTiffLocation    = []; %'J:/DATA/registered'; %'C:/DATA/registered'; % leave empty to NOT save registered tiffs

% registration options
ops0.doRegistration         = 1;
ops0.showTargetRegistration = 1;
ops0.PhaseCorrelation       = 1; % set to 0 for non-whitened cross-correlation
ops0.SubPixel               = Inf; % 2 is alignment by 0.5 pixel, Inf is the exact number from phase correlation
ops0.NimgFirstRegistration  = 500;
ops0.nimgbegend             = 0; % frames to average at beginning and end of blocks
ops0.dobidi                 = 0; % infer and apply bidirectional phase offset

% cell detection options
ops0.ShowCellMap            = 1;
ops0.sig                    = 0.5;  % spatial smoothing length in pixels; encourages localized clusters
ops0.nSVDforROI             = 500;
ops0.NavgFramesSVD          = 5000; % how many (pooled) frames to do the SVD based on
ops0.signalExtraction       = 'raw'; % how to extract ROI and neuropil signals: 'raw', 'regression'

ops0.getSVDcomps            = 0;
ops0.nSVD                   = 1000; % how many SVD components to keep

ops0.RegFileBinLocation     = 'J:/DATA/registered';


% ops0.clustModel  = 'neuropil';
% ops0.neuropilSub = 'model';
% ops0.Nk0                    = 600; %1300;  % how many clusters to start with
% ops0.Nk                     = 200; %650;  % how many clusters to end with
% clustrules.diameter         = 5; % expected diameter of cells (used for 0.25 * pi/4*diam^2 < npixels < 10*pi/4*diam^2)

%%
for iexp = 1:length(db)        %3:length(db)   
    gpuDevice(1);
    run_pipeline(db(iexp), ops0);
%     run_REDaddon(iexp, db, ops0);
end
%%
