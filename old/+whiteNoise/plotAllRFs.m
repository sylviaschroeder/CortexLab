function plotAllRFs(RFinfo)

figure
hold on

for neuron = 1:length(RFinfo)
    color = 'r';
    if RFinfo(neuron).RFsign == -1
        color = 'b';
    end
    whiteNoise.plotFitReceptiveField(size(RFinfo(neuron).RF_separated), ...
        RFinfo(neuron).gaussianFit, RFinfo(neuron).stimPosition, ...
        RFinfo(neuron).RFtype, RFinfo(neuron).RFsign);
end