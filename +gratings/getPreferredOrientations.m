function [neuronOrientations, neuronDirections, neuronBoth, ...
    scaledResponsesOri, scaledResponsesDir, orientations] = ...
    getPreferredOrientations(calciumTraces, samplingRate, ...
    stimulusMatrix, directions, blank, timeAfterOnset, timeAfterOffset)

limitsInFrames = round([timeAfterOnset timeAfterOffset] * samplingRate);
stimTraces = ssLocal.getTracesPerStimulus(calciumTraces, stimulusMatrix, ...
    limitsInFrames);
timeMeans = nanmean(stimTraces, 4);
trialMedians = nanmedian(timeMeans, 3);


directionVectors = directions(:,1) ./ 180 .* pi;
directionVectors = [cos(directionVectors), sin(directionVectors)];
stimOrientations = unique(mod(directions(:,1), 180));
stimOrientations = [stimOrientations, NaN(length(stimOrientations),2)];
for k = 1:size(stimOrientations,1)
    ind = find(mod(directions(:,1), 180) == stimOrientations(k,1));
    stimOrientations(k,2:length(ind)+1) = ind;
end
orientationVectors = mod(stimOrientations(:,1), 180) / 90 * pi;
orientationVectors = [cos(orientationVectors), sin(orientationVectors)];


neuronOrientations = nan(size(calciumTraces,2), 2); 
% 1st col: vector orientation, 2nd col: vector length
neuronDirections = nan(size(calciumTraces, 2), 2);
% 1st col: vector orientation, 2nd col: vector length

scaledResponsesOri = NaN(size(trialMedians,1), size(stimOrientations,1));
scaledResponsesDir = NaN(size(trialMedians,1), size(directions,1));
orientations = stimOrientations(:,1);

for roi = 1:size(calciumTraces,2)
    responses = trialMedians(roi, setdiff(1:size(trialMedians,2),blank))';
    % subtract response to blank stimulus
    responses = responses - trialMedians(roi, blank);
    % check whether neuron is mostly suppressed by visual stimulation
    if max(responses) < -min(responses)
        % if so, multiply by -1
        responses = -responses;
    end
    responsesOri = mean(responses(stimOrientations(:,2:3)), 2);
    responsesDir = responses(directions(:,2));
    
    % set negative responses to zero
    responsesOri(responsesOri < 0) = 0;
    responsesDir(responsesDir < 0) = 0;
    % scale responses to max. 1
    responsesOri = responsesOri ./ max(responsesOri);
    responsesDir = responsesDir ./ max(responsesDir);
    scaledResponsesOri(roi, :) = responsesOri;
    scaledResponsesDir(roi, :) = responsesDir;
    responseVectors = orientationVectors .* repmat(responsesOri, 1, 2);
    vectorMean = sum(responseVectors, 1) ./ ...
        sum(sqrt(sum(responseVectors .^ 2, 2)));
    neuronOrientations(roi,:) = [mod(atan2(vectorMean(2), vectorMean(1)) / ...
        pi * 90, 180), sqrt(sum(vectorMean .^ 2))];
    
    responseVectors = directionVectors .* repmat(responsesDir, 1, 2);
    vectorMean = sum(responseVectors, 1) ./ ...
        sum(sqrt(sum(responseVectors .^ 2, 2)));
    neuronDirections(roi,:) = [mod(atan2(vectorMean(2), vectorMean(1)) / ...
        pi * 180, 360), sqrt(sum(vectorMean .^ 2))];
end
neuronBoth = neuronOrientations;
ind = neuronDirections(:,2) > neuronOrientations(:,2);
neuronBoth(ind,:) = neuronDirections(ind,:);
neuronBoth(:,1) = mod(neuronBoth(:,1), 180);