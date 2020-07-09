function [centerY, height] = adjustCenterHeight(output, F, params)

pupilTop = output(:,2:4); % [x y likelihood]
pupilBottom = output(:,5:7);
pupilLeft = output(:,8:10);
pupilRight = output(:,11:13);
lidTop = output(:,14:16);
lidBottom = output(:,17:19);

width = pupilRight(:,1) - pupilLeft(:,1);
height = pupilBottom(:,2) - pupilTop(:,2);
centerX = pupilLeft(:,1) + 0.5 .* width;

% if either top or bottom of pupil uncertain or if top/bottom of lid is too
% close to top/bottom of pupil: (1) estimate height of pupil from width, 
% (2) find vertical center of pupil (by adding half height to bottom or top
% of pupil)
onlyBottomValid = pupilTop(:,3) < params.minCertainty & ...
    pupilBottom(:,3) > params.minCertainty;
onlyTopValid = pupilTop(:,3) > params.minCertainty & ...
    pupilBottom(:,3) < params.minCertainty;
lidNearPupilTop = pupilTop(:,2) - lidTop(:,2) < params.maxDistPupilLid;
lidNearPupilBottom = lidBottom(:,2) - pupilBottom(:,2) < params.maxDistPupilLid;
widthValid = pupilLeft(:,3) > params.minCertainty & ...
    pupilRight(:,3) > params.minCertainty;
ind = widthValid & (onlyBottomValid | onlyTopValid | ...
    lidNearPupilTop | lidNearPupilBottom);
height(ind) = F(centerX(ind), width(ind));
centerY = pupilBottom(:,2) - 0.5 .* height;
centerY(widthValid & onlyTopValid) = pupilTop(widthValid & onlyTopValid,2) + ...
    0.5 .* height(widthValid & onlyTopValid);