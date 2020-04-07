function f = orituneWrappedConditions(pars,oris,conditions) 
% ORITUNEWRAPPEDTWOCONDITIONS 	sum of two wrapped gaussians living on a 
% circle, for orientation tuning of two different curves
% 
%		orituneWrappedTwoConditions([Dp, Rp, DI, Ro, sigma, ...
%           deltaDp, deltaRp, deltaDI, deltaRo, deltaSigma], oris), where
%		Dp is the preferred direction (bet 0 and 360)
%		Rp is the response to the preferred direction;
%		DI is the direction index (Rp-Rn)/(Rp+Rn) where Rn is response to
%           null direction
%		Ro is the background response (useful only in some cases)
%		sigma is the tuning width;
%       delta[x] is the difference of parameter x for the second tuning
%           curve
% 		oris are the orientations, bet 0 and 360, the NaN entry separates
% 		orientations applied to the first curve from those applied to the
% 		second curve
%
%		orituneWrapped([Dp, Rp, Rn, Ro, sigma], 'spiel'), gives you a
%       little description of the tuning.
%
% See also: ORITUNE, FITORI, CIRCSTATS, CIRCSTATS360

% 2015 Sylvia Schroeder (based on oritune)
% part of the Matteobox toolbox


if nargin < 3
    conditions = ones(size(oris));
end
conds = unique(conditions);
numCond = length(conds);
if length(pars(:))~= 5*numCond
    error('You must specify all 5 parameters AND all deltas for each condition');
end
pars(end+1 : 5*numCond) = NaN;

Dp = pars(1) + [0 pars(6:5:end)];
Rp = pars(2) + [0 pars(7:5:end)];
Rn = [(1-pars(3))/(1+pars(3))*Rp(1) ...
    (1-sum([pars(3) pars(8:5:end)]))./(1+sum([pars(3) pars(8:5:end)])) .* Rp(2:end)];
Ro = pars(4) + [0 pars(9:5:end)];
sigma = pars(5) + [0 pars(10:5:end)];

f = NaN(size(oris));

if ischar(oris) && strcmp(oris,'spiel')
    disp(['Width at half height is ' num2str(sqrt(log(2)*2)*sigma,3) ' deg']);
    disp(['Direction index is ' num2str((Rp-Rn)/(Rp+Rn),3)]);
else
    for c = 1:numCond
        ind = conditions == conds(c);
        orisC = reshape(oris(ind),1,[]);
        gauss1 = Rp(c).*sum(exp((-bsxfun(@plus,orisC-Dp(c),(-5:5)'.*360).^2)./...
            (2*sigma(c)^2)));
        gauss2 = Rn(c).*sum(exp((-bsxfun(@plus,orisC-Dp(c)+180,(-5:5)'.*360).^2)./...
            (2*sigma(c)^2)));
        
        f(ind) = Ro(c) + gauss1 + gauss2;
    end
end
