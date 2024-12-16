%% Folders:
folderTools = 'C:\dev\toolboxes';
folderScript = 'C:\dev\workspaces\CortexLab';
folderBase = 'C:\Users\Sylvia\OneDrive - University of Sussex';

folderROIData_Sylvia = fullfile(folderBase, 'Lab\DATA\InfoStructs');
folderROIData_Vane = fullfile(folderBase, 'Projects\2023_OrientationColumns\RawData');
folderDataSave = fullfile(folderBase, 'Projects\2023_OrientationColumns\DataToPublish');

%% Add paths
addpath(genpath(fullfile(folderScript)));
addpath(genpath(fullfile(folderTools, 'npy-matlab', 'npy-matlab')));
addpath(fullfile(folderTools, 'Suite2P_matlab', 'cortexlab'));

%% Load database
db = db_columns_vane;

%% Convert data
experiments = {'expGratings', 'expNoise'};

for k = 19:length(db)
    fprintf('Dataset %d of %d\n', k, length(db))
    subject = db(k).subject;
    ind = strfind(subject, 'SS');
    subject = subject(ind:end);
    if str2double(subject(end-1:end)) > 50 % boutons
        folderSave = fullfile(folderDataSave, 'boutons');
    else
        folderSave = fullfile(folderDataSave, 'neurons');
    end
    folderSession = fullfile(folderSave, subject, db(k).date);
    if ~isfolder(folderSession)
        mkdir(folderSession)
    end

    planePerUnit = readNPY(fullfile(folderSession, '_ss_2pRois._ss_2pPlanes.npy'));
    cellIDs = readNPY(fullfile(folderSession, '_ss_2pRois.ids.npy'));
    masks = NaN(length(cellIDs), 500);
    maxMaskElems = 0;
    fovs = NaN(length(db(k).planes), 2);
    fovBoundaries = NaN(length(db(k).planes), 4);
    meanImages = [];
    
    for exp = 1:length(experiments)
        if ~isfield(db, experiments{exp}) || isempty(db(k).(experiments{exp}))
            continue
        end
        
        vane = false;
        folder = fullfile(folderROIData_Sylvia, db(k).subject, ...
            db(k).date, num2str(db(k).(experiments{exp})));
        if ~isfolder(folder)
            vane = true;
            folder = fullfile(folderROIData_Vane, db(k).subject, ...
                db(k).date, 'metaStructs', num2str(db(k).(experiments{exp})));
        end
        if ~isfolder(folder)
            continue
        end
        file = [sprintf('%s_%d_%s', db(k).date, db(k).(experiments{exp}), ...
            db(k).subject) '_2P_plane%03d_ROI.mat'];
        for iPlane = 1:length(db(k).planes)
            indPlane = find(planePerUnit == iPlane);
            ids = cellIDs(indPlane);
            load data
            d = load(fullfile(folder, sprintf(file, db(k).planes(iPlane))));
            meta = d.meta;
            for iCell = 1:length(ids)
                m = meta.ROI.CellMaps{ids(iCell)};
                numM = length(m);
                if numM > maxMaskElems
                    maxMaskElems = numM;
                end
                if numM > size(masks,2)
                    masks = padarray(masks, [0 numM-size(masks,2)], NaN, 'post');
                end
                masks(indPlane(iCell),1:length(m)) = m;
            end
            
            if vane
                fovs(iPlane,:) = size(meta.targetFrame);
            else
                fovs(iPlane,:) = [length(meta.validY) length(meta.validX)];
            end
            yb = meta.validY([1 end]);
            xb = meta.validX([1 end]);
            fovBoundaries(iPlane,:) = [yb(:)' xb(:)'];
            meanImages = cat(3, meanImages, meta.meanFrame);
        end
        masks(:,maxMaskElems+1:end) = [];
        writeNPY(masks, fullfile(folderSession, '_ss_2pRois.masks.npy'));
        writeNPY(fovs, fullfile(folderSession, '_ss_2pPlanes.fovSizePix.npy'));
        um2pix = infoPixUm(size(meta.targetFrame,1), meta.zoomFactor, meta.microID);
        fovsUm = fovs ./ [um2pix.yPU um2pix.xPU];
        writeNPY(fovsUm, fullfile(folderSession, '_ss_2pPlanes.fovSizeMicrons.npy'));
        writeNPY(fovBoundaries, fullfile(folderSession, '_ss_2pPlanes.fovBoundariesPix.npy'));
        meanImages = permute(meanImages, [3 1 2]);
        writeNPY(meanImages, fullfile(folderSession, '_ss_2pPlanes.meanFrame.npy'));
        break
    end
end
