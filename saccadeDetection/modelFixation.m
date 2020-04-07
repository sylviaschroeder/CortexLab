function [fix_mean,fix_pos] = modelFixation(saccades,eyepos)

%modelFixation models the eye position in between the saccades
%   [fix_mean,fix_pos] = modelFixation(saccades,eyepos)
%
%   modelFixation compute the eye average position in between saccades
%
%   OUTPUTS
%   fix_mean: 2 x (nsaccades - 1) array with avearage position for each
%       period in between the saccades
%   fix_pos: 2 x t, each row is the eye position in Deg along one of
%       the axis (x = 1, y =2) averaged between saccade periods
%   
%   INPUTS
%   saccades: time of saccades
%   eyepos: 2 x t matrix, each row is the eye position in Deg along on of
%       the axis (x = 1, y =2)


nf = size(eyepos,2);
nsaccades=length(saccades) - 1;

if saccades(1) ~= 1
    saccades = [1 saccades];
end
if saccades(end) ~= nf
    saccades = [saccades nf];
end

fix_pos = nan(2,nf);
fix_mean = nan(2,nsaccades);

for istep = 1:nsaccades
    
    
    fix_period = (saccades(istep):saccades(istep+1)-1);
    if length(fix_period)> 1
        fix_mean(:,istep) = mean(eyepos(:,fix_period),2);
        fix_pos(:,fix_period) = fix_mean(:,istep)*ones(1,length(fix_period));
    end
    
end