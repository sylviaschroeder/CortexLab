function [orientations, spatPhases] = getStaticGratingPars(stimulusSequence)

orientationString = 'ori';
phaseString = 'sph';
contrastString = 'cg';
    
orientations = zeros(length(stimulusSequence.labels),1);
spatPhases = zeros(length(stimulusSequence.labels),1);
for stim = 1:length(stimulusSequence.labels)
    description = stimulusSequence.labels{stim};
    str = description(strfind(description, [orientationString ' = ']):end);
    oriPars = sscanf(str, '%*s = %d');
    str = description(strfind(description, [phaseString ' = ']):end);
    phPars = sscanf(str, '%*s = %d');
    str = description(strfind(description, [contrastString ' = ']):end);
    contrastPars = sscanf(str, '%*s = %d');
    if contrastPars ~= 0 % contrast is not 0
        orientations(stim) = oriPars;
        spatPhases(stim) = phPars;
    else
        orientations(stim) = NaN;
        spatPhases(stim) = NaN;
    end
end