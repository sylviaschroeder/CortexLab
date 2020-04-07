function [prefOri, selectOri, prefDir, selectDir] = ...
    getPreferenceAndSelectivity(responses, directions, blankResponse)

directionVectors = directions ./ 180 .* pi;
directionVectors = [cos(directionVectors), sin(directionVectors)];
stimOrientations = unique(mod(directions, 180));
stimOrientations = [stimOrientations, NaN(length(stimOrientations),2)];
for k = 1:size(stimOrientations,1)
    ind = find(mod(directions, 180) == stimOrientations(k,1));
    stimOrientations(k,2:length(ind)+1) = ind;
end
orientationVectors = mod(stimOrientations(:,1), 180) / 90 * pi;
orientationVectors = [cos(orientationVectors), sin(orientationVectors)];

% subtract response to blank stimulus or minimum stim response
responses = responses - min([responses(:);blankResponse]);
% % check whether neuron is mostly suppressed by visual stimulation
% if max(responses) < -min(responses)
%     % if so, multiply by -1
%     responses = -responses;
% end
responsesOri = mean(responses(stimOrientations(:,2:3)), 2);
responsesDir = responses(:);

% set negative responses to zero
responsesOri(responsesOri < 0) = 0;
responsesDir(responsesDir < 0) = 0;
% scale responses to max. 1
responsesOri = responsesOri ./ max(responsesOri);
responsesDir = responsesDir ./ max(responsesDir);
responseVectors = orientationVectors .* repmat(responsesOri, 1, 2);
vectorMean = sum(responseVectors, 1) ./ ...
    sum(sqrt(sum(responseVectors .^ 2, 2)));
prefOri = mod(atan2(vectorMean(2), vectorMean(1)) / pi * 90, 180);
selectOri = sqrt(sum(vectorMean .^ 2));

responseVectors = directionVectors .* repmat(responsesDir, 1, 2);
vectorMean = sum(responseVectors, 1) ./ ...
    sum(sqrt(sum(responseVectors .^ 2, 2)));
prefDir = mod(atan2(vectorMean(2), vectorMean(1)) / pi * 180, 360);
selectDir = sqrt(sum(vectorMean .^ 2));