
function [vel, acc] = computeVelocity(pos, smoothSize, Fs)
% function [vel, acc] = computeVelocity(pos, smoothSize)
%
% assumes pos is uniformly sampled in time
%
% smooth size is in units of samples of position

pos = pos(:)';

% area of this smoothing window is 1 so total values are unchanged - units
% don't change
smoothWin = gausswin(smoothSize)./sum(gausswin(smoothSize));

vel = [0 conv(diff(pos), smoothWin, 'same')];
if nargin > 2
    vel = vel .* Fs; % multiply by Fs to get cm/sec
end

if nargout>1
    % here we choose to apply the smoothing again - it's sort of
    % empirically necessary since derivatives amplify noise. 
    acc = [0 conv(diff(vel), smoothWin, 'same')]; 
    if nargin > 2
        acc = acc .* Fs; %cm/sec^2
    end
end