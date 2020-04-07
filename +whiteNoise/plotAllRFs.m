function plotAllRFs(RFinfo, RFtype)

figure
hold on

n = 0;
for neuron = 1:length(RFinfo)
    if isempty(RFinfo(neuron).RFs)
        continue
    end
    if nargin > 1
        ind = find(strcmp(RFtype, {RFinfo(neuron).RFs.RFtype}));
    else
        ind = 1:length(RFinfo(neuron).RFs);
    end
    for k = ind
        RF = RFinfo(neuron).RFs(k);
        whiteNoise.plotFitReceptiveField(size(RF.RF_separated), ...
            RF.gaussianFit, RF.stimPosition, ...
            RF.RFtype, RF.RFsign);
    end
    n = n + length(ind);
end

title(sprintf('n = %d', n))