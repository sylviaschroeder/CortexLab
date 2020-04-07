function orientationDifferences = getOrientationDiffs(prefOrientations, type)

% prefOrientations  [neurons x 1]
% type              'ori' or 'dir'

% tuningThreshold = 0.3;
% 
% [neuronOrientations, neuronDirections, neuronBoth] = ...
%     gratings.getPreferredOrientations(calciumTraces, samplingRate, ...
%     stimulusMatrix, directions, blank, timeAfterOnset, timeAfterOffset);
% prefOrientations = neuronOrientations(:,1);
% prefOrientations(neuronOrientations(:,2) < tuningThreshold) = NaN;
% prefDirections = neuronDirections(:,1);
% prefDirections(neuronDirections(:,2) < tuningThreshold) = NaN;
% prefOriOrDir = neuronBoth(:,1);
% prefOriOrDir(neuronBoth(:,2) < tuningThreshold) = NaN;

orientationDifferences = abs(bsxfun(@minus, prefOrientations, prefOrientations'));
if strcmp(type,'ori')
    ind = orientationDifferences > 90;
    orientationDifferences(ind) = 180 - orientationDifferences(ind);
elseif strcmp(type, 'dir')
    ind = orientationDifferences > 180;
    orientationDifferences(ind) = 360 - orientationDifferences(ind);
else
    display('type has to be ori or dir')
    orientationDifferences = [];
    return
end
% directionDifferences = abs(bsxfun(@minus, prefDirections, prefDirections'));
% ind = directionDifferences > 180;
% directionDifferences(ind) = 360 - directionDifferences(ind);
% oriOrDirDifferences = abs(bsxfun(@minus, prefOriOrDir, prefOriOrDir'));
% ind = oriOrDirDifferences > 90;
% oriOrDirDifferences(ind) = 180 - oriOrDirDifferences(ind);