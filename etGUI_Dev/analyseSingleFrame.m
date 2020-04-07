function [paramsOut, frameThresh] = analyseSingleFrame(frame, params)

gaussStd = params.gaussStd;
diskR = params.diskR;
sat = params.sat;
th = params.th;

hGauss = fspecial('gaussian', ceil(5*gaussStd), gaussStd);
hDisk = fspecial('disk', diskR);
frameFiltered = imfilter(frame, imfilter(hGauss, hDisk, 0, 'full'), 'replicate');
% frameFiltered = imfilter(frameFiltered, hDisk, 'replicate');

% normalize the frame
% c = class(frameFiltered);
frameFiltered = double(frameFiltered);
% frameFiltered = (frameFiltered-double(intmin(c)))/double(intmax(c));
frameFiltered = frameFiltered - min(frameFiltered(:));
frameFiltered = frameFiltered/max(frameFiltered(:));

% saturate the frame
frameThresh = frameFiltered;
frameThresh(frameThresh>sat)=sat;
% and then normalize once again
frameThresh = frameThresh/sat;

% threshold  = 0.5;
% frameThresh(frameThresh>=threshold) = 1;
% frameThresh(frameThresh<threshold) = 0;
% % frameThresh = logical(frameThresh);


bw = frameThresh <= th; % here the pupil is white
ebw = edge(bw);
chList = regionprops(bw, 'ConvexHull');
convex = false(size(bw));
for i = 1:length(chList)
    convex = convex | poly2mask(chList(i).ConvexHull(:, 1), chList(i).ConvexHull(:, 2),...
        size(convex, 1), size(convex, 2));
end
[convexLabeled, n] = bwlabel(convex, 8);
% if we got more than one object, select the most central
if n>1
    C = regionprops(convexLabeled, 'Centroid');
    delta = bsxfun(@minus, fliplr(size(bw)/2), cell2mat({C(:).Centroid}'));
    d2 = sum(delta.^2, 2);
    [~, iL] = min(d2);
    convexLabeled = convexLabeled==iL;
end

edgeConvexLabeled = edge(convexLabeled);
% only take the edge points that are close to original edges
edgeRes = ebw & imdilate(edgeConvexLabeled, ones(3));
[yy, xx] = ind2sub(size(bw), find(edgeRes));
nPoints = length(yy);

if nPoints >= 5

    % now let's fit the ellipse
    % we are basically solving an equation for the quadratic curve in 2-D
    % of the form: Ax^2+Bxy+Cy^2+Dx+Ey=1
    X = [xx.^2, xx.*yy, yy.^2, xx, yy];
    params = pinv(X)*ones(nPoints, 1);
    
    % now iteratively remove some of the edge points which fall inside the ellipse, but
    % leave at least  minNPoints points
    minNPoints = 10;
    minPrctile = min(minNPoints/nPoints*100, 100);
    err = (X*params-1);%.^2;
    idx = err<prctile(err, max(80, minPrctile));
    nPoints = sum(idx);
    params = pinv(X(idx, :))*ones(nPoints, 1);
    
    err = (X*params-1);%.^2;
    idx = err<prctile(err, max(60, minPrctile));
    nPoints = sum(idx);
    params = pinv(X(idx, :))*ones(nPoints, 1);
    
    err = (X*params-1);%.^2;
    idx = err<prctile(err, max(40, minPrctile));
    nPoints = sum(idx);
    params = pinv(X(idx, :))*ones(nPoints, 1);

    ell = abc2ellipse(params);
    
else
    % in case there were no contours detected at all (might happen during a
    % blink, for example) assign nans where applicable
    xx = [];
    yy = [];
    ell = abc2ellipse(NaN(5, 1));
end


paramsOut.x0 = ell.x0;
paramsOut.y0 = ell.y0;
paramsOut.aAxis = ell.a;
paramsOut.bAxis = ell.b;
paramsOut.abAxis = NaN; % don;t have this number in the canonic not rotated form
paramsOut.area = ell.area;
paramsOut.good = ell.isEllipse;
paramsOut.xx = xx;
paramsOut.yy = yy;
paramsOut.equation = ell.eq;
paramsOut.xxEll = ell.xx;
paramsOut.yyEll = ell.yy;

