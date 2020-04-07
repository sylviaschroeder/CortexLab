function plotRFtypeStats(RFinfo)

types = {'Absolute', 'ON', 'OFF'};
counts = zeros(length(RFinfo), length(types));

for neuron = 1:length(RFinfo)
    if isempty(RFinfo(neuron).RFs)
        continue
    end
    for type = 1:length(types)
        if any(strcmp({RFinfo(neuron).RFs.RFtype}, types{type}))
            counts(neuron,type) = 1;
        end
    end
end

% total = [sum(counts(:,1),1), sum(all(counts(:,[1 2])==1, 2),1), ...
%     sum(all(counts(:,[1 3])==1, 2),1), sum(all(counts==1, 2),1)];
total = [sum(counts(:,1)==1 & any(counts(:,2:3)==1,2),1), ... % Absolute and at least one type
    sum(all(counts(:,[1 2])==1, 2) & counts(:,3)==0,1), ...   % only ON
    sum(all(counts(:,[1 3])==1, 2) & counts(:,2)==0,1), ...   % only OFF
    sum(all(counts==1, 2),1)];                                % both ON and OFF

figure
bar(total ./ length(RFinfo) .* 100)
title(sprintf('RF types (%d neurons)', length(RFinfo)))
ylabel('Neurons (in %)')
set(gca, 'XTickLabel', [types, {'ON+OFF'}])

fprintf('Neurons total: %d\nAbsolute RFs: %d\nOFF RFs: %d\nON RFs: %d\nOFF+ON RFs: %d\n', ...
    length(RFinfo), total(1), total(2), total(3), total(4))