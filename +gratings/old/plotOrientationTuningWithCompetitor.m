function plotOrientationTuningWithCompetitor(calciumTraces, samplingRate, ...
    stimulusMatrix, stimulusSequence, neuronIDs)

% Define time points of signal to consider (in sec)
limitsAfterOnset = [1 5];

limitsInFrames = round(limitsAfterOnset * samplingRate);
stimTraces = ssLocal.getTracesPerStimulus(calciumTraces, stimulusMatrix, ...
    limitsInFrames);
kernelSize = size(stimTraces, 4);
kernelTimes = ((0:kernelSize-1)-limitsInFrames(1)) / samplingRate;
kernelTimesPatch  = [kernelTimes, fliplr(kernelTimes)]';
stimDuration = size(stimTraces,4) - sum(limitsInFrames);

if nargin < 5
    neuronIDs = 1:size(calciumTraces, 2);
end

% Parse stimulus information
parameters = NaN(length(stimulusSequence.labels), 4);
for stim = 1:length(stimulusSequence.labels)
    description = stimulusSequence.labels{stim};
    parameters(stim,:) = sscanf(description, ...
        'ori1 = %d, cg1 = %d, ori2 = %d, cg2 = %d');
end
orientationsOut = unique(parameters(:,3))';
orientations = unique([parameters(:,1); parameters(:,3)])';
contrastsRF = setdiff(unique(parameters(:,2))', 0);
contrastsOut = setdiff(unique(parameters(:,4))', 0);
sortedParameterCombis = NaN(length(contrastsRF) + length(contrastsOut) + ...
    length(contrastsRF) * length(contrastsOut) * length(orientationsOut), ...
    3);
% first: contrast of competitor is 0
sortedParameterCombis(1:length(contrastsRF), :) = ...
    [contrastsRF(end:-1:1)', repmat([-1 0], length(contrastsRF), 1)];
k = length(contrastsRF);
% second: contrast in RF is 0
sortedParameterCombis(k+(1:length(contrastsOut)), :) = ...
    [repmat([0 -1], length(contrastsOut), 1), contrastsOut];
k = k + length(contrastsOut);
% third: all other parameter combinations
pars = [reshape(repmat(contrastsRF, ...
    length(contrastsOut) * length(orientationsOut), 1), [], 1), ...
    reshape(repmat(...
    reshape(repmat(orientationsOut, length(contrastsOut), 1), [], 1), ...
    1, length(contrastsRF)), [], 1), ...
    reshape(repmat(contrastsOut, 1, ...
    length(contrastsRF) * length(orientationsOut)), [], 1)];
sortedParameterCombis(k+1:end, :) = pars;

% match stimuli to parameter combinations
stimIDs = NaN(size(sortedParameterCombis, 1), length(orientations));
for stim = 1:size(parameters, 1)
    pars = parameters(stim,:);
    if pars(4) == 0 % contrast of competitor is 0
        pars(3) = -1; % set orientation of competitor to -1
    end
    rowInd = all(repmat(pars(2:4),size(sortedParameterCombis,1),1) == ...
        sortedParameterCombis, 2);
    if pars(2) == 0 % contrast in RF is 0
        colInd = find(orientations == pars(3)); % take orientation of competitor
    else 
        colInd = find(orientations == pars(1));
    end
    stimIDs(rowInd, colInd) = stim;
end
stimIDs(all(isnan(stimIDs), 2),:) = [];
labels = {'Contrast In: %d %%', 'Ori Out: %d deg', 'Contrast Out: %d %%'};

for roi = 1:size(calciumTraces,2)
    figure('Position', [790 190 1105 890])
    mini = min(reshape(stimTraces(roi,:,:,:),[],1));
    maxi = max(reshape(stimTraces(roi,:,:,:),[],1));
    for row = 1:size(stimIDs,1)
        for col = 1:size(stimIDs,2)
            subplot(size(stimIDs,1), size(stimIDs,2), (row-1)*size(stimIDs,2)+col)
            if row == 1
                title(sprintf('Orientation: %d deg', round(orientations(col))))
            end
            if col == 1
                ind = stimIDs(row,:);
                ind(isnan(ind)) = [];
                pars = parameters(ind,2:4);
                consistent = all(repmat(pars(1,:),size(pars,1),1) == pars, 1);
                label = '';
                for p = 1:3
                    if consistent(p)
                        label = [label sprintf([labels{p} '\n'], round(pars(1,p)))];
                    end
                end                
                ylabel(label)
            else
                set(gca, 'YTick', [])
            end
            if row == size(stimIDs,1)
                xlim(kernelTimes([1 end]))
                xlabel('Time (in s)')
            else
                set(gca, 'XTick', [])
            end
            if isnan(stimIDs(row,col))
                continue
            end
            hold on
            plot([0 0], [mini maxi], 'r--')
            plot([1 1]*stimDuration/samplingRate, [mini maxi], 'r--')
            
            resp = squeeze(stimTraces(roi, stimIDs(row,col), :, :));
            patch(repmat(kernelTimesPatch, 1, size(resp, 1)), ...
                [resp, fliplr(resp)]', 'k', 'EdgeColor', 'k', ...
                'EdgeAlpha', 0.3, 'FaceColor', 'none');
            
            kernel = nanmean(resp, 1);
            plot(kernelTimes, kernel, 'k', 'LineWidth', 2)
            
            xlim(kernelTimes([1 end]))
            ylim([mini maxi])
        end
    end
    annotation('textbox', [0 0.9 1 0.1], 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'LineStyle', 'none', ...
        'String', ['ROI ' num2str(neuronIDs(roi))])
    
    drawnow
end