
fieldOfView = [2 500 520; 2.5 400 430; 3 335 370];

% animal(1).name = 'M141002_SS026';
% animal(1).expDate(1).name = '2014-10-29';
% animal(1).expDate(1).exp(1).ID = 4;
% animal(1).expDate(1).exp(1).planes = 2:3;

% animal(1).name = 'M141003_SS027';
% animal(1).expDate(1).name = '2014-10-23';
% animal(1).expDate(1).exp(1).ID = 2;
% animal(1).expDate(1).exp(1).planes = 1:4;
% 
% animal(1).expDate(1).name = '2014-10-30';
% animal(1).expDate(1).exp(1).ID = 3;
% animal(1).expDate(1).exp(1).planes = 2:4;
% 
animal(1).name = 'M141007_SS029';
animal(1).expDate(1).name = '2014-11-12';
animal(1).expDate(1).exp(1).ID = 4;
animal(1).expDate(1).exp(1).planes = 2:3;

oriDiffs = {};
dirDiffs = {};
combinedDiffs = {};
distances = {};

for an = 1:length(animal)
    for session = 1:length(animal(an).expDate);
        for exp = 1:length(animal(an).expDate(session).exp)
            info = ppbox.infoPopulate(animal(an).name, ...
                animal(an).expDate(session).name, ...
                animal(an).expDate(session).exp(exp).ID);
            for plane = animal(an).expDate(session).exp(exp).planes
                filePath = fullfile(info.folderProcessed, ...
                    sprintf('%s_plane%03d_ROI', info.basename2p, plane));
                % load meta
                load(filePath, 'meta')
                [~, stimSequence, stimMatrix, ~, samplingRate] = ...
                    ssLocal.getStimulusResponseInfo(meta);
                mapSizePixels = [length(meta.validY), length(meta.validX)];
                ind = fieldOfView(:,1) == info.zoomFactor;
                mapSizeMM = mapSizePixels .* fieldOfView(ind,2:3) ./ ...
                    size(meta.targetFrame);
                ind = strcmp(meta.ROI.CellClasses, 's');
                [oriDiffs{end+1}, dirDiffs{end+1}, combinedDiffs{end+1}] = ...
                    gratings.getOrientationDiffs(meta.F(:, ind), ...
                    samplingRate, stimMatrix, stimSequence);
                distances{end+1} = spatial.getNeuronDistances( ...
                    meta.ROI.CellMaps(ind), mapSizePixels, mapSizeMM);
            end
        end
    end
end

gratings.plotOriDiffsVsDistance(oriDiffs, distances, 'orientation');
gratings.plotOriDiffsVsDistance(dirDiffs, distances, 'direction');
gratings.plotOriDiffsVsDistance(combinedDiffs, distances, 'orientation (combined)');

% publish('plotOriDiffsVsDistance_multipleDatasets.m','outputDir','\Results_4_oriDiffsVsDist', 'maxOutputLines', 0)