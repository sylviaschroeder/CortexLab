function [orientations, blank] = getOrientations(stimulusSequence)

orientations = zeros(0,2); %1st col: orientation, 2nd col: ID of stimulus
for stim = 1:length(stimulusSequence.labels)
    description = stimulusSequence.labels{stim};
    parameters = sscanf(description, 'dir = %d, g = %d');
    if parameters(2) == 0 % bar is black (as opposed to gray like background)
        orientations(end+1,:) = [parameters(1), stim];
    else
        blank = stim;
    end
end
orientations = sortrows(orientations);