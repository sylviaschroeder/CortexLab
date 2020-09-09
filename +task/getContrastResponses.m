function [meanResp, semResp, cL, cR] = getContrastResponses(responses, ...
    contrastL, contrastR, pupil)

cond = unique(pupil);
cL = unique(contrastL);
cR = unique(contrastR);

meanResp = NaN(length(cL), length(cR), length(cond));
semResp = NaN(length(cL), length(cR), length(cond));

for c = 1:length(cond)
    for l = 1:length(cL)
        for r = 1:length(cR)
            ind = pupil == cond(c) & contrastL == cL(l) & contrastR == cR(r);
            meanResp(l,r,c) = mean(responses(ind));
            semResp(l,r,c) = std(responses(ind)) ./ sqrt(sum(ind));
        end
    end
end