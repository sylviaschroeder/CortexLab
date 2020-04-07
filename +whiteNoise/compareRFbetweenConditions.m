function RFinfo = compareRFbetweenConditions(calciumTrace, traceTimes, ...
    conditionValues, conditionTimes, condLabels, stimFrames, ...
    stimFrameTimes, stimTimes, stimPosition, RFtypes, method, titleLine)

colors = [0 0 0; 1 0 0];

condValuesResampled = interp1(conditionTimes, conditionValues, ...
    traceTimes, 'linear', NaN);
condValuesResampled = condValuesResampled >= 0.5;

[trialTraces, trialTimes] = whiteNoise.getTrialTraces(calciumTrace, traceTimes, ...
    stimTimes, stimFrameTimes, 0);
trialConds = cell(1,2);
for cond = 1:2
    trialConds{cond} = whiteNoise.getTrialTraces(condValuesResampled(:,cond), traceTimes, ...
        stimTimes, stimFrameTimes, 0);
end

figure('Position', [565 680 1320 420])
hold on
indNoClass = true(size(trialTraces));
handles = zeros(1,3);
for cond = 1:2
    trials = trialTraces;
    trials(trialConds{cond} ~= 1) = NaN;
    h = plot(trialTimes, trials, 'Color', colors(cond,:));
    handles(cond) = h(1);
    indNoClass(trialConds{cond} == 1) = false;
end
trials = trialTraces;
trials(~indNoClass) = NaN;
h = plot(trialTimes, trials, 'Color', [0.5 0.5 0.5]);
handles(3) = h(1);
axis tight
legend(handles, [condLabels, {'unclassified'}])
title(titleLine)
xlabel('Time (in s)')
ylabel('Calcium response')

condTrace = calciumTrace;
condTrace(condValuesResampled(:,1) ~= 1) = NaN;
[receptiveFields, RFtimes] = whiteNoise.getReceptiveField(...
    condTrace, traceTimes, stimFrames, stimFrameTimes, stimTimes, ...
    RFtypes, 1, method);
smoothedRFs = whiteNoise.plotReceptiveField(receptiveFields, RFtimes, ...
    stimPosition, RFtypes);

RFinds = [];
for type = 1:length(RFtypes)
    if max(abs(smoothedRFs{type}(:))) / std(smoothedRFs{type}(:)) > 5
        RFinfo(type) = whiteNoise.fitGaussianToRF(smoothedRFs{type}, ...
            RFtimes, stimPosition, 0, RFtypes{type}, 0);
        RFinds = [RFinds, type];
    end
end
if ~isempty(RFinds)
    RFinfo = RFinfo(RFinds);
    figure
    whiteNoise.plotFitReceptiveField([size(smoothedRFs{1},1) ...
        size(smoothedRFs{1},2)], ...
        {RFinfo.gaussianFit}, stimPosition, RFtypes(RFinds), [RFinfo.RFsign]);
end