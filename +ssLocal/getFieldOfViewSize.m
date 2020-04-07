function [sizeHoriz, sizeVert] = getFieldOfViewSize(zoomFactor)

zoom = [2 2.2 2.3 2.5 3];
measuredHoriz = [524 491 460 431.5 371];
measuredVert = [503.5 452.5 430 401 337.5];

ind = find(zoom == zoomFactor);
if ~isempty(ind)
    sizeHoriz = measuredHoriz(ind);
    sizeVert = measuredVert(ind);
else
    curve = fit(zoom, measuredHoriz, 'poly2');
    sizeHoriz = curve(zoomFactor);
    curve = fit(zoom, measuredVert, 'poly2');
    sizeVert = curve(zoomFactor);
end