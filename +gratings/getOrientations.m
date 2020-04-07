function [orientations, blanks] = getOrientations(stimulusSequence, ...
    orientationString, contrastString)

if nargin < 2 || isempty(orientationString)
    orientationString = 'ori';
end
if nargin < 3 || isempty(contrastString)
    contrastString = 'cg';
end

orientations = zeros(0,2); %1st col: orientation, 2nd col: ID of stimulus
blanks = [];
for stim = 1:length(stimulusSequence.labels)
    description = stimulusSequence.labels{stim};
    str = description(strfind(description, [orientationString ' = ']):end);
    oriPars = sscanf(str, '%*s = %d');
    str = description(strfind(description, [contrastString ' = ']):end);
    contrastPars = sscanf(str, '%*s = %d');
    if contrastPars ~= 0 % contrast is not 0
        orientations(end+1,:) = [oriPars, stim];
    else
        blanks(end+1) = stim;
    end
end
orientations = sortrows(orientations);