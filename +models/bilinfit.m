function [coeffs, x0y0] = bilinfit(xx, yy)

% for some reason I think this problem is called bilinear optimization
% problem, hence the name of the function
% 
% [coeffs, x0y0] = bilinfit(xx, yy) will fit the data [xx,yy] with
% multiple linear lines, under the constraint that all of them need to
% intersect in a single point (x0y0).
% xx, yy - cell arrays of data vectors (each cell corresponds to a data for a
% single linear fit)
% coeffs - Nx2 array of coefficients [ai, bi], so that yy = ai*xx+bi for
% each group

%% reshape the data (just in case)
xx = xx(:);
yy = yy(:);

nLines = length(xx);

if nLines<3
    warning('This function is intended for use only if you need to solve more than 2 linear equations');
    warning('Otherwise the solution is trivial (any two straight lines intersect in a single point)');
    return;
end

for iLine = 1:nLines
    xx{iLine} = xx{iLine}(:);
    yy{iLine} = yy{iLine}(:);
end

%% put data into matrices and vectors
yyVector = cell2mat(yy);
nSamples = length(yyVector);
xxMatrix = zeros(nSamples, nLines);
indStart = 1;
for iLine = 1:nLines
    indEnd = indStart + length(xx{iLine}) - 1;
    idx = indStart:indEnd;
    xxMatrix(idx, iLine) = xx{iLine};
    indStart = indStart + length(xx{iLine});
end

%% Now we will solve the equation y = ai(x-x0)+y0 to find optimal ai, x0 and y0
% We will use a mixed optimization-LSE approach: optimizing for x0, and solving the
% linear regression for ai and y0.

onesVector = ones(nSamples, 1);
% A = [xxMatrix-x0, onesVector];
eps = @(x0) norm(yyVector - [xxMatrix-x0.*(xxMatrix~=0), onesVector]*pinv([xxMatrix-x0.*(xxMatrix~=0), onesVector])*yyVector);

% getting an initial guess for x0 (not perfect, just something reasonable)
x0Init = median(cell2mat(xx));

options = optimset('Display', 'off');
[x0Opt, err] = fminsearch(eps, x0Init, options);
% [x0Opt, err] = fminunc(eps, x0Init, options);

% And finally, let's get to y = ai*x+bi form
% from y = ai(x-x0)+y0 by substituting bi = y0-ai*x0
params = pinv([xxMatrix-x0Opt.*(xxMatrix~=0), onesVector])*yyVector;
y0 = params(end);
ai = params(1:end-1);
bi = y0-ai*x0Opt;

%% Finishing and exiting

coeffs = [ai(:), bi(:)];
x0y0 = [x0Opt, y0];

