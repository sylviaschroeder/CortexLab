%% Folders
% UCL
% folderData = '\\zubjects.cortexlab.net\Subjects';
% folderLocal = 'J:\Ephys';
% folderCode = 'C:\dev\workspace\CortexLab';
% folderTools = 'C:\STORAGE\workspaces';
% folderToolsNew = 'C:\dev\toolboxes';
% folderConfig = 'C:\dev\workspace\CortexLab\kilosortFiles\configFiles';
% Sussex
folderData = 'E:\Subjects';
folderLocal = 'D:\Ephys';
folderCode = 'C:\dev\workspace\CortexLab';
folderTools = 'C:\dev\tools';
folderToolsNew = 'C:\dev\tools';
folderConfig = 'C:\dev\workspace\CortexLab\kilosortFiles\configFiles';

%% Parameters
numChans = 385;
chanMap = 'neuropixPhase3A_kilosortChanMap.mat';

%% Add paths
% addpath(genpath(fullfile(folderTools, 'Kilosort2.0')))
addpath(genpath(fullfile(folderToolsNew, 'Kilosort3.0')))
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(genpath(fullfile(folderCode)))

%% Load database
db = db_ephys_task;

for k = 5
    for pr = 1 %:length(db(k).probes)
        probe = db(k).probes{pr};
        fprintf('\n%s %s %s\n', db(k).subject, db(k).date, probe)
        %% Download data
        fLocal = fullfile(folderLocal, db(k).subject, db(k).date, probe);
        if ~isfolder(fLocal)
            mkdir(fLocal)
        end
        fServer = fullfile(folderData, db(k).subject, db(k).date, ...
            ['ephys_' probe]);
        file = dir(fullfile(fServer, '*.imec.ap_CAR.bin'));
%         
%         status = copyfile(fullfile(fServer, file(1).name), fLocal);
%         if ~status
%             fprintf('  !!! NOTE: file not copied; data not sorted!\n')
%             continue
%         end
        
        %% Set parameters
        run(fullfile(folderConfig, 'configFile384.m'))
        % parameters used for ephys in task with 4 SC probes:
        ops.trange = [0 Inf]; % time range to sort
        ops.NchanTOT = numChans; % total number of channels in your recording
        ops.fproc = fullfile(fLocal, 'temp_wh.dat'); % proc file on a fast SSD
        ops.chanMap = fullfile(folderConfig, chanMap);
        ops.fbinary = fullfile(fLocal, file(1).name);
        % 20.04.21: try different sets of parameters
        ops.minFR = 0;
        
        % main parameter changes from Kilosort2 to v2.5
        ops.sig        = 20;  % spatial smoothness constant for registration
        ops.fshigh     = 300; % high-pass more aggresively
        ops.nblocks    = 5; % blocks for registration. 0 turns it off, 1 does rigid registration. Replaces "datashift" option.
        
        % main parameter changes from Kilosort2.5 to v3.0
        ops.Th       = [9 9];
        
        %% Run spike sorting algorithm: kilosort 3.0
        rez                = preprocessDataSub(ops);
        rez                = datashift2(rez, 1);
        
        [rez, st3, tF]     = extract_spikes(rez);
        
        rez                = template_learning(rez, tF, st3);
        
        [rez, st3, tF]     = trackAndSort(rez);
        
        rez                = final_clustering(rez, tF, st3);
        
        rez                = find_merges(rez, 1);
        
        rootZ = fullfile(fLocal, 'kilosort3');
        mkdir(fLocal)
        rezToPhy2(rez, fLocal);
        
        %% Run spike sorting algorithm: kilosort 2.0
%         % preprocess data to create temp_wh.dat
%         rez = preprocessDataSub(ops);
%         
%         % time-reordering as a function of drift
%         rez = clusterSingleBatches(rez);
%         
%         % saving here is a good idea, because the rest can be resumed after loading rez
%         save(fullfile(fLocal, 'rez.mat'), 'rez', '-v7.3');
%         
%         % main tracking and template matching algorithm
%         rez = learnAndSolve8b(rez);
%         
%         % final merges
%         rez = find_merges(rez, 1);
%         
%         % final splits by SVD
%         rez = splitAllClusters(rez, 1);
%         
%         % final splits by amplitudes
%         rez = splitAllClusters(rez, 0);
%         
%         % decide on cutoff
%         rez = set_cutoff(rez);
%         
%         % write to Phy
%         rezToPhy(rez, fLocal);
%         disp('  finished sorting. now copying to server')
%         
%         %% Copy results to server and delete local files
%         delete(fullfile(fLocal, file(1).name))
%         delete(fullfile(fLocal, 'temp_wh.dat'))
%         
%         fSave = fullfile(fServer, 'sorting');
%         if isfolder(fSave)
%             movefile(fSave, fullfile(fServer, 'sorting_old'));
%         end
%         if ~isfolder(fSave)
%             mkdir(fSave)
%         end
%         copyfile(fLocal, fSave)
%         rmdir(fLocal, 's')
%         
%         close all
%         clear rez
    end
end