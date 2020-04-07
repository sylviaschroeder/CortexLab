function plotChoiceMatrix(psychoDat, taskParams)

if nargin < 2
    taskParams = [];
end

%% Calculate psychometric data
% rwdTypes = unique(rewardSize(rewardSize>0)); %Array of unique reward types
% rwd = rewardSize';
% perf = (perf / histc(repeatNum,1))*100; %Performance = number of correct / frequency of repeatNum 1 * 100
ind = psychoDat.repeatNum == 1 & (diff(psychoDat.contrast,1,1) ~= 0 | ...
    all(psychoDat.contrast == 0, 1));
perf = sum(psychoDat.feedback(ind)==1) / sum(ind)*100;
numTrials = length(psychoDat.resp);
numFalseNogos = sum(any(psychoDat.contrast>0,1) & psychoDat.resp==3);

if isempty(taskParams)
    contrastsLeft = unique(psychoDat.contrast(1,:));
    contrastsRight = unique(psychoDat.contrast(2,:));
else
    contrastsLeft = unique(taskParams.visCueContrast(1,:));
    contrastsRight = unique(taskParams.visCueContrast(2,:));
end
respTypes = unique(psychoDat.resp(psychoDat.resp>0));

matrix = NaN(length(contrastsLeft), length(contrastsRight), length(respTypes));
for c = 1:length(respTypes)
    for l = 1:length(contrastsLeft)
        for r = 1:length(contrastsRight)
            ind = psychoDat.resp == respTypes(c) & ...
                psychoDat.contrast(1,:) == contrastsLeft(l) & ...
                psychoDat.contrast(2,:) == contrastsRight(r) & ...
                psychoDat.repeatNum == 1;
            matrix(l,r,c) = sum(ind);
        end
    end
end

% normalise by occurrences of each stimulus
matrix = matrix ./ sum(matrix,3);
matrix(isnan(matrix)) = -1/199;

cm = gray(199);
cols = lines(1);
cm = [cols; cm];

figure('Position', [180 674 1680 420])
titles = {'left', 'right', 'nogo'};
for p = 1:3
    subplot(1,3,p)
    imagesc(matrix(:,:,p), [-1/199,1])
    title(['Choice ' titles{p}])
    if p == 1
        ylabel('Contrast left')
    end
    xlabel('Contrast right')
    set(gca, 'YTick', 1:length(contrastsLeft), 'YTickLabel', contrastsLeft, ...
        'XTick', 1:length(contrastsRight), 'XTickLabel', contrastsRight, ...
        'YDir', 'normal')
end
colormap(cm)

annotation('textbox', 'String', sprintf('%s\nPerformance: %.1f%%\nNo. trials: %d\n(w/o false NoGos: %d)', ...
    psychoDat.expRef{1}, perf, numTrials, numTrials-numFalseNogos), 'Interpreter', 'none', ...
    'Position', [.01 .7 .1 .25], 'LineStyle', 'none')