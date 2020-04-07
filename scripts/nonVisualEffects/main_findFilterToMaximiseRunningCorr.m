%% Load database
db_sparseNoise
% db_boutons_sparseNoise

%% Folders
% loading data
folderROIData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\InfoStructs';
% saving data
folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\receptiveFields\SC neurons';
% folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\boutons\receptiveFields';

%% Parameters
% parameters to filter running trace before regressing it from calcium
% traces
runKrnlTimes = [-5 5];
runLambdas = reshape([1 5]' * logspace(-3, 0, 4), 1, []);

%% Correlate running with calcium traces

for k = 1:length(db)
    fprintf('Dataset %d: %s %s exp.: %d\n', k, db(k).subject, ...
        db(k).date, db(k).exp);
    folder = fullfile(folderROIData, db(k).subject, ...
        db(k).date, num2str(db(k).exp));
    file = [sprintf('%s_%d_%s',db(k).date, db(k).exp, ...
        db(k).subject) '_2P_plane%03d_ROI.mat'];
    
    traces = [];
    numNeurons = NaN(1, length(db(k).planes));
    for iPlane = 1:length(db(k).planes)
        % load meta
        data=load(fullfile(folder, sprintf(file,db(k).planes(iPlane))));
        meta = data.meta;
        meta.folderTL = strrep(meta.folderTL, 'ioo.ucl.ac.uk', ...
            'cortexlab.net');
        
        % get stimulus information
        if iPlane == 1
            % get trace time
            traceTimes = ppbox.getFrameTimes(meta);
            traceBin = median(diff(traceTimes));
            
            % get running speed
            ballData = nonVis.getRunningSpeed(meta);
            if isempty(ballData)
                disp('No ball data!')
                return
            end
            runSpeed = ballData.total / median(diff(ballData.t)) / 53;
            runTime = ballData.t;
            runBin = median(diff(runTime));
            % median filter running trace
            numBins = ceil(traceBin/runBin);
            runSpeed = medfilt1(runSpeed, max(5, floor(numBins/2)*2+1));
            runSpeed = smooth(runSpeed, numBins);
            runSpeed = interp1(runTime, runSpeed, traceTimes, 'pchip')';
        end
        
        %get calcium traces
        tr = meta.F_final;
        % set data of duplicate neurons or unhealthy neurons to NaN
        tr(:, meta.ROI.isDuplicate==1 | meta.ROI.isSwitchOn == 1) = NaN;
        if iPlane == 1
            traces = [traces, tr];
        else
            tr_int = NaN(length(traceTimes), size(tr,2));
            t = ppbox.getFrameTimes(meta);
            for n = 1:size(tr,2)
                if all(isnan(tr(:,n)))
                    continue
                end
                tr_int(:,n) = interp1(t, tr, traceTimes, 'pchip');
            end
            traces = [traces, tr_int];
        end
        numNeurons(iPlane) = size(tr,2);
    end

   [A, n, wins] =  krnl.getToeplitz(traceTimes, [], [], {runSpeed}, ...
       {runKrnlTimes}, true);
   A = A - mean(A);
   tracesNorm = (traces - nanmean(traces,1)) ./ nanstd(traces,0,1);
   krnls = NaN(size(A,2), size(traces,2), length(runLambdas));
   preds = NaN(size(traces,1), size(traces,2), length(runLambdas));
   lamMatrix = krnl.makeLambdaMatrix(size(A,2), 1);
   lams = sqrt(runLambdas .* size(A,1));
   for lam = 1:length(runLambdas)
       krnls(:,:,lam) = [A; lamMatrix .* lams(lam)] \ ...
           [traces; zeros(size(lamMatrix,1), size(traces,2))];
       preds(:,:,lam) = A*krnls(:,:,lam);
   end
end