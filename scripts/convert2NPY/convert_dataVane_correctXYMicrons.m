dataFolder = 'C:\Users\Sylvia\OneDrive - University of Sussex\Projects\2023_OrientationColumns\DataToPublish';

%% Load database
db = db_columns_vane;

for k = 1:length(db)
    if ~db(k).newSuite2p
        continue
    end
    name = db(k).subject;
    ind = strfind(name, 'SS');
    name = name(ind:end);
    date = db(k).date;
    f = fullfile(dataFolder, 'neurons', name, date);
    
    % load data
    data = io.getCalciumData(f);
    planes = data.planes;
    recInfo = io.getRecordingInfo(f);
    brainPos = recInfo.roiPositions;

    uniPlanes = unique(planes);
    for p = 1:length(uniPlanes)
        units = find(planes == uniPlanes(p));
        for iUnit = 1:length(units)
            m = recInfo.roiMasks(units(iUnit),:);
            m(isnan(m)) = [];
            [ySub, xSub] = ind2sub(flip(recInfo.fovPix(p,:)), m);
            xCentre = mean(xSub);
            yCentre = mean(ySub);
            brainPos(units(iUnit),1) = xCentre / recInfo.fovPix(p,2) * ...
                recInfo.fovMicrons(p,2);
            brainPos(units(iUnit),2) = yCentre / recInfo.fovPix(p,1) * ...
                recInfo.fovMicrons(p,1);
        end
    end
    brainPos(:,3) = NaN;

    writeNPY(brainPos, fullfile(f, '_ss_2pRois.xyz.npy'));
end