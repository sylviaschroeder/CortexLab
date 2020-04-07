function saccades = findSaccades(eyepos,detected_saccades,step_threshold)

%findSaccades finds saccades during eye movements
%   saccades = findSaccades(eyepos,interp_pos,step_threshold)
%
%   findSaccades looks for zero-crossing of the acceleration (second
%   derivate) between periods of high positive and negative acceleration
%
%   OUTPUTS
%   saccades: time of saccades
%   
%   INPUTS
%   eyepos: 2 x t matrix, each row is the eye position in Deg along one of
%       the axis (x = 1, y =2)
%   detected_saccades: time of saccades that have already been detected
%   step_threshold: threshold 

nf = size(eyepos,2);
crossSteep = nan(2,nf-3);
for ix = 1:2
    eye_vel = diff(eyepos(ix,:)); %velocity
    eye_acc = diff(eye_vel); %acceleration
    crossSteep(ix,:) = (eye_acc(1:nf-3)).*(eye_acc(2:nf-2));
end


zeroCross_p = cell(1,2);
for ix = 1:2
    zeroCross_p{ix} = 2 + find(crossSteep(ix,:) < -step_threshold);
end
saccades = sort(unique([zeroCross_p{1} zeroCross_p{2} find(detected_saccades)]));