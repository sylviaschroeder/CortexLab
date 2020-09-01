%% Folders
folderData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\NPY\task_2p';
folderTools = 'C:\STORAGE\workspaces';
folderCode = 'C:\dev\workspace\CortexLab';
folderPlots = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\wheelTask\pupil_neuralActivity';

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(fullfile(folderCode))

%% Get event-triggered traces of SVs
% events: stimuli of different contrasts, (beep), reward time
% differentiate between: go left + correct, go right + correct, nogo +
% correct, go left + incorrect, go right + incorrect, nogo + incorrect
stimWin = [-1 3];
beepWin = [-1 3];
fbWin = [-1 4];
leftCol = [0 0 1];
rightCol = [1 0 0];
evLabels = {'stim onset','beep','feedback'};
numSVs = 10;
numSmooth = 5; % 10

subjDirs = dir(fullfile(folderData, 'SS*'));
% for subj = 1 %:length(subjDirs)
subj = 1;
name = subjDirs(subj).name;
dateDirs = dir(fullfile(folderData, name, '2*'));
%     for dt = 1 %:length(dateDirs)
dt = 1;
date = dateDirs(dt).name;
%         if ~isfile(fullfile(folderData, name, date, 'eye.diameter.npy'))
%             continue
%         end
%         subjects{end+1,1} = name;
%         dates{end+1,1} = date;

caTraces = readNPY(fullfile(folderData, name, date, '_ss_2pCalcium.dff.npy'));
caTime = readNPY(fullfile(folderData, name, date, '_ss_2pCalcium.timestamps.npy'));
planes = readNPY(fullfile(folderData, name, date, '_ss_2pRois._ss_2pPlanes.npy'));
planeDelays = readNPY(fullfile(folderData, name, date, '_ss_2pPlanes.delay.npy'));
pupil = readNPY(fullfile(folderData, name, date, 'eye.diameter.npy'));
pupilTime = readNPY(fullfile(folderData, name, date, 'eye.timestamps.npy'));
stimT = readNPY(fullfile(folderData, name, date, '_ibl_trials.stim_intervals.npy'));
contrastL = readNPY(fullfile(folderData, name, date, '_ibl_trials.contrastLeft.npy'));
contrastR = readNPY(fullfile(folderData, name, date, '_ibl_trials.contrastRight.npy'));
beeps = readNPY(fullfile(folderData, name, date, '_ibl_trials.goCue_times.npy'));
choice = readNPY(fullfile(folderData, name, date, '_ibl_trials.choice.npy'));
fb = readNPY(fullfile(folderData, name, date, '_ibl_trials.feedbackType.npy'));
fbTimes = readNPY(fullfile(folderData, name, date, '_ibl_trials.feedback_times.npy'));

% get SVs
if numSmooth > 1
    caSmooth = smoothdata(caTraces, 1, 'movmean', numSmooth);
else
    caSmooth = caTraces;
end
[U, S, V] = svd(caSmooth - nanmean(caSmooth,1), 'econ');
% invert SVs as neccessary
signs = sign(mean(U)-median(U));
U = U .* signs;
V = V .* signs;
% plot singular values
figure
hold on
stem(cumsum(diag(S).^2./sum(diag(S).^2)), 'k')
stem(diag(S).^2./sum(diag(S).^2), 'r')
xlim([0 30])
xlabel('Singular value')
ylabel('Proportion of variance explained')

% plot trajectories of SVs
figure
hold on
plot(caTime, U(:,1:numSVs) + (0:-1:-numSVs+1).*.05)
% plot(repmat(stimT(:,1)',2,1), [.01;-.5], 'k')

% plot trajectories of first 2 SVs
figure
plot(U(:,1), U(:,2))
% plot3(U(:,1), U(:,2), U(:,3))
xlabel('SV 1')
ylabel('SV 2')
% zlabel('SV 3')

% sort stimuli according to contrast
contrasts = [contrastL, contrastR];
contr_unq = unique(contrasts, 'rows');
stimCols = zeros(size(contr_unq,1),3,2);
stimCols(:,:,1) = contr_unq(:,1).*leftCol + (1-contr_unq(:,1)).*[1 1 1];
stimCols(contr_unq(:,1)==0,:,1) = NaN;
stimCols(:,:,2) = contr_unq(:,2).*rightCol + (1-contr_unq(:,2)).*[1 1 1];
stimCols(contr_unq(:,2)==0,:,2) = NaN;
stimCols = nanmean(stimCols,3);
stimCols(all(isnan(stimCols),2),:) = 0;

%% get Toeplitz matrix for stim onset, beep, feedback time
stimOnsets = cell(1, size(contr_unq,1));
for c = 1:size(contr_unq,1)
    stimOnsets{c} = stimT(all(contrasts == contr_unq(c,:), 2),1);
end

[A, numSmpl, windows] = krnl.getToeplitz(caTime, ...
    [stimOnsets, beeps, fbTimes], [repmat({stimWin},1,size(contr_unq,1)), ...
    {beepWin}, {fbWin}]);
[ETAs, traces] = krnl.getETA(U(:,1:numSVs) * S(1:numSVs,1:numSVs), A, numSmpl);
stimETAs = cat(3, ETAs{1:size(contr_unq,1)});
[m_time,ind_time] = max(abs(stimETAs),[],1);
[m_stim, ind_stim] = max(squeeze(m_time), [], 2);

% get time range of beep after stim onset
diffs = beeps - stimT(:,1);
beepRange = [min(diffs) max(diffs)];
beepMeanDiff = mean(diffs);
% get time range of feedback after beep
diffs = fbTimes - beeps;
fbRange = [min(diffs) max(diffs)];
fbMeanDiff = mean(diffs);
% get time range of stim onset after feedback
diffs = stimT(2:end,1) - fbTimes(1:end-1);
stimRange = [min(diffs) max(diffs)];
stimMeanDiff = mean(diffs);

% plot ETAs
for sv = 1:numSVs
    sgn = sign(stimETAs(ind_time(1,sv,ind_stim(sv)),sv,ind_stim(sv)));
    ax = zeros(1,3);
    figure('position', [5 560 1914 420])
    t = tiledlayout(1,3);
    ax(1) = nexttile;
    hold on
    mini = min(reshape(stimETAs(:,sv,:).*sgn,[],1));
    maxi = max(reshape(stimETAs(:,sv,:).*sgn,[],1));
    top = maxi + 0.04 * (maxi - mini);
    for c = 1:size(contr_unq,1)
        plot(windows{c}, ETAs{c}(:,sv).*sgn, 'Color', stimCols(c,:), 'LineWidth', 2)
    end
    h = plot(beepRange, [top top], 'Color', [1 1 1].*0.5, 'LineWidth', 4);
    plot(beepMeanDiff, top, 'x', 'Color', [1 1 1].*0.5, 'LineWidth', 4, ...
        'MarkerSize', 12);
    xlim(windows{1}([1 end]))
    ylabel('F')
    xlabel(sprintf('Time from %s', evLabels{1}))
    title(sprintf('Singular vector %d', sv))
    legend(h, 'beep time')
    
    ax(2) = nexttile;
    hold on
    plot(windows{end-1}, ETAs{end-1}(:,sv), 'k', 'LineWidth', 2)
    h = plot(fbRange, [top top], 'Color', [1 1 1].*0.5, 'LineWidth', 4);
    plot(fbMeanDiff, top, 'x', 'Color', [1 1 1].*0.5, 'LineWidth', 4, ...
        'MarkerSize', 12);
    xlim(windows{end-1}([1 end]))
    xlabel(sprintf('Time from %s', evLabels{2}))
    legend(h, 'feedback time')
    
    ax(3) = nexttile;
    hold on
    plot(windows{end}, ETAs{end}(:,sv), 'k', 'LineWidth', 2)
    h = plot(stimRange, [top top], 'Color', [1 1 1].*0.5, 'LineWidth', 4);
    plot(stimMeanDiff, top, 'x', 'Color', [1 1 1].*0.5, 'LineWidth', 4, ...
        'MarkerSize', 12);
    xlim(windows{end}([1 end]))
    xlabel(sprintf('Time from %s', evLabels{3}))
    legend(h, 'stim time')
    linkaxes(ax, 'y')
end

%% Plot traces (not ETAs) sorted by time

%% get Toeplitz matrix for stim onset for different performances
% Go & correct, Go & incorrect side, Go & no stim
% NoGo & correct, NoGo & incorrect
conds = [choice<3 & [fb==1, fb==0 & sum(contrasts,2)>0, fb==0 & sum(contrasts,2)==0], ...
    choice==3 & [fb==1, fb==0]];
stimOnsets = cell(5, size(contr_unq,1));
for c = 1:size(contr_unq,1)
    indSt = all(contrasts == contr_unq(c,:), 2);
    for cond = 1:size(conds,2)
        stimOnsets{cond,c} = stimT(indSt & conds(:,cond),1);
    end
end

[A, numSmpl, windows] = krnl.getToeplitz(caTime, ...
    reshape(stimOnsets,1,[]), repmat({stimWin},1,numel(stimOnsets)));
[ETAs, traces] = krnl.getETA(U(:,1:numSVs) * S(1:numSVs,1:numSVs), A, numSmpl);
ETAs = reshape(ETAs, size(stimOnsets,1), size(stimOnsets,2));
traces = reshape(traces, size(stimOnsets,1), size(stimOnsets,2));

numCols = ceil(sqrt(size(contr_unq,1)));
numRows = ceil(size(contr_unq,1) / numCols);

% get time range of beep after stim onset
diffs = beeps - stimT(:,1);
beepRange = [min(diffs) max(diffs)];
beepMeanDiff = mean(diffs);

% plot ETAs
cols = lines(size(stimOnsets,1));
lbls = {'corr.','incorr. dir.','incorr.','corr.','incorr. no move'};
for sv = 1:numSVs
    figure('position', [5 45 1910 950])
    t = tiledlayout(numCols, numRows);
    ax = zeros(1, size(contr_unq,1));
    lb = cell(size(lbls));
    for st = 1:size(contr_unq,1)
        h = NaN(1,size(stimOnsets,1));
        ax(st) = nexttile;
        hold on
        for cond = 1:size(stimOnsets,1)
            resp = ETAs{cond,st}(:,sv);
            if all(isnan(resp))
                continue
            end
            h(cond) = plot(windows{1}, resp, 'Color', cols(cond,:), 'LineWidth', 2);
            lb{cond} = [lbls{cond} ' (' num2str(size(traces{cond,st},3))];
        end
        title(sprintf('Contrasts: %.2f vs %.2f', contr_unq(st,1), contr_unq(st,2)))
        ind = ~isnan(h);
        legend(h(ind),lb(ind))
    end
    linkaxes(ax, 'y')
    annotation('textbox', [0 .97 1 .02], 'String', ['SV ' num2str(sv)], ...
        'LineStyle', 'none', 'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold', 'FontSize', 16)
end


%% Align data to task events
pre = 2;
post = 3;
numSVs = 5;
gap = 3;
sortTime = [0 2];
% sortTime = [-1 0];

% subjects = cell(0,1);
% dates = cell(0,1);
% algnmnts = cell(0, 3);
% timeAlgn = cell(0,1);
% binSz = NaN(0,1);
% stimTimes = cell(0,1);
% contrasts = cell(0,1);
% choices = cell(0,1);
% feedbacks = cell(0,1);

subjDirs = dir(fullfile(folderData, 'SS*'));
% for subj = 1 %:length(subjDirs)
    subj = 1;
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folderData, name, '2*'));
%     for dt = 1 %:length(dateDirs)
        dt = 1;
        date = dateDirs(dt).name;
%         if ~isfile(fullfile(folderData, name, date, 'eye.diameter.npy'))
%             continue
%         end
%         subjects{end+1,1} = name;
%         dates{end+1,1} = date;
        
        caTraces = readNPY(fullfile(folderData, name, date, '_ss_2pCalcium.dff.npy'));
        caTime = readNPY(fullfile(folderData, name, date, '_ss_2pCalcium.timestamps.npy'));
        planes = readNPY(fullfile(folderData, name, date, '_ss_2pRois._ss_2pPlanes.npy'));
        planeDelays = readNPY(fullfile(folderData, name, date, '_ss_2pPlanes.delay.npy'));
        pupil = readNPY(fullfile(folderData, name, date, 'eye.diameter.npy'));
        pupilTime = readNPY(fullfile(folderData, name, date, 'eye.timestamps.npy'));
        stimT = readNPY(fullfile(folderData, name, date, '_ibl_trials.stim_intervals.npy'));
        contrastL = readNPY(fullfile(folderData, name, date, '_ibl_trials.contrastLeft.npy'));
        contrastR = readNPY(fullfile(folderData, name, date, '_ibl_trials.contrastRight.npy'));
        beeps = readNPY(fullfile(folderData, name, date, '_ibl_trials.goCue_times.npy'));
        choice = readNPY(fullfile(folderData, name, date, '_ibl_trials.choice.npy'));
        fb = readNPY(fullfile(folderData, name, date, '_ibl_trials.feedbackType.npy'));
        fbTimes = readNPY(fullfile(folderData, name, date, '_ibl_trials.feedback_times.npy'));
        
%         stimTimes{end+1,1} = stimT;
%         contrasts{end+1,1} = [contrastL, contrastR];
%         choices{end+1,1} = choice;
%         feedbacks{end+1,1} = fb;
%         feedbackTimes{end+1,1} = fbTimes;
        
        % population activity
        [U, S, V] = svd(caTraces - nanmean(caTraces,1), 'econ');

        sampleDur = diff(caTime(1:2));
        pre_ca = round(pre / sampleDur);
        post_ca = round(post / sampleDur);
        timeAlgn_ca = (-pre_ca:post_ca) .* sampleDur;
        
        sampleDur = diff(pupilTime(1:2));
        pre_ppl = round(pre / sampleDur);
        post_ppl = round(post / sampleDur);
%         binSz(end+1,1) = sampleDur;
        timeAlgn_ppl = (-pre_ppl:post_ppl) .* sampleDur;
        
        contrasts = [contrastL, contrastR];
        contr_unq = unique(contrasts, 'rows');
        
        algn_ppl = cell(size(contr_unq,1),1);
        algn_ca = cell(size(contr_unq,1),1);
        for c = 1:size(contr_unq,1)
            trials = find(all(contrasts == contr_unq(c,:), 2)); % & stimT(:,1) > 1500);
            % pupil
            algn_ppl{c} = NaN(length(trials), pre_ppl+post_ppl+1);
            ind = NaN(length(trials),1);
            for t = 1:length(trials)
                ind(t) = find(pupilTime >= stimT(trials(t),1), 1, 'first');
            end
            ind = ind + (-pre_ppl : post_ppl);
            ind2 = ind;
            invalid = ind < 1 | ind > length(pupilTime);
            ind2(invalid) = 1;
            al = pupil(ind2);
            al(invalid) = NaN;
            algn_ppl{c} = al;
            % Singular values
            algn_ca{c} = NaN(length(trials), pre_ca+post_ca+1, numSVs);
            ind = NaN(length(trials),1);
            for t = 1:length(trials)
                ind(t) = find(caTime >= stimT(trials(t),1), 1, 'first');
            end
            ind = ind + (-pre_ca : post_ca);
            ind2 = ind;
            invalid = ind < 1 | ind > length(caTime);
            ind2(invalid) = 1;
            for sv = 1:numSVs
                tr = U(:,sv);
                al = tr(ind2);
                al(invalid) = NaN;
                algn_ca{c}(:,:,sv) = al;
            end
        end
        negSVs = nanmean(nanmean(algn_ca{end}(:,pre_ca:end,:),1),2) < 0;
        inv = ones(numSVs,1);
        inv(negSVs) = -1;
        for c = 1:size(contr_unq,1)
            algn_ca{c} = algn_ca{c} .* permute(inv, [2 3 1]);
        end
        
%% Plot neural activity and pupil aligned to stim onset
        a = find(timeAlgn_ca >= sortTime(1), 1);
        o = find(timeAlgn_ca >= sortTime(2), 1);
        for sv = 1:numSVs
            ca = [];
            ppl = [];
            for c = 1:size(contr_unq,1)
                [~,order] = sort(mean(algn_ca{c}(:,a:o,sv),2));
                ca = [ca; algn_ca{c}(order,:,sv); NaN(gap,size(algn_ca{c},2))];
                ppl = [ppl; algn_ppl{c}(order,:); NaN(gap,size(algn_ppl{c},2))];
            end
            ca(end-gap+1:end,:) = [];
            ppl(end-gap+1:end,:) = [];
            figure
            subplot(1,2,1)
            imagesc(timeAlgn_ca([1 end]), [1 size(ca,1)], ca)
            
            subplot(1,2,2)
            imagesc(timeAlgn_ppl([1 end]), [1 size(ppl,1)], ppl)
        end
%     end
% end

%% Relationship between neural response and pupil at stimulus presentation
% neural measures: baseline (before stim), deviation from mean stimulus
% response (baseline corrected)
% pupil measures: baseline (before stim), during stim, diff. baseline -
% stim
timeBsln = [-.5 0];
timeStim = [0 1];

indBsln_ppl = find(timeAlgn_ppl >= timeBsln(1), 1) : find(timeAlgn_ppl >= timeBsln(2), 1);
indStim_ppl = find(timeAlgn_ppl >= timeStim(1), 1) : find(timeAlgn_ppl >= timeStim(2), 1);
indBsln_ca = find(timeAlgn_ca >= timeBsln(1), 1) : find(timeAlgn_ca >= timeBsln(2), 1);
indStim_ca = find(timeAlgn_ca >= timeStim(1), 1) : find(timeAlgn_ca >= timeStim(2), 1);

bsln_ppl = cell(size(algn_ppl));
stim_ppl = cell(size(algn_ppl));
for c = 1:size(contr_unq,1)
    bsln_ppl{c} = nanmean(algn_ppl{c}(:,indBsln_ppl),2);
    stim_ppl{c} = nanmean(algn_ppl{c}(:,indStim_ppl),2);
end
bsln_ca = cell(size(algn_ca));
stim_ca = cell(size(algn_ca));
for c = 1:size(contr_unq,1)
    for sv = 1:numSVs
        bsln_ca{c} = squeeze(nanmean(algn_ca{c}(:,indBsln_ca,:),2));
        stim_ca{c} = squeeze(nanmean(algn_ca{c}(:,indStim_ca,:),2));
    end
end

% pupil and neural responses for different stimuli
data = {bsln_ppl, stim_ppl, cellfun(@minus, stim_ppl, bsln_ppl, 'UniformOutput', false);
    bsln_ca, stim_ca, cellfun(@minus, stim_ca, bsln_ca, 'UniformOutput', false)};
dataNames = {'Pupil', 'Neural SV'};
measureNames = {'Baseline', '@Stim', '@Stim - Bsln.'};
for d = 1:2
    for sv = 1:size(data{d,1}{1},2)
        ax = zeros(1,3);
        figure('position', [5 560 1914 420])
        t = tiledlayout(1,3);
        for t = 1:3
            nexttile
            m = cellfun(@nanmean, data{d,t}, 'UniformOutput', false);
            m = cat(1, m{:});
            sem = cellfun(@nanstd, data{d,t}, 'UniformOutput', false);
            sem = cat(1, sem{:}) ./ ...
                sqrt(cellfun(@size, bsln_ppl, num2cell(ones(length(bsln_ppl),1))));
            errorbar(m(:,sv), sem(:,sv), 'ok', 'MarkerFaceColor', 'k')
            xlim([0 length(algn_ppl)+1])
            xlabel('Different stimuli')
            if d == 1
                lbl = sprintf('Pupil: %s', measureNames{t});
            else
                lbl = sprintf('SV %d: %s', sv, measureNames{t});
            end
            ylabel(lbl)
            ax(t) = gca;
        end
        linkaxes(ax(1:2), 'y')
    end
end


bc = cat(1,bsln_ca{:});
sc = cat(1,stim_ca{:});
dc = cellfun(@minus, stim_ca, bsln_ca, 'UniformOutput', false);
dc = cellfun(@minus, dc, cellfun(@mean, dc, 'UniformOutput', false), ...
    'UniformOutput', false);
dc = cat(1,dc{:});
for sv = 1:numSVs
    figure('position', [35 375 1100 600])
    t = tiledlayout(2,3);
    
    nexttile
    scatter(cat(1,bsln_ppl{:}), bc(:,sv))
    xlabel('Pupil baseline')
    ylabel('Neural baseline')
    title(sprintf('SV %d',sv))
    
    nexttile
    scatter(cat(1,stim_ppl{:}), bc(:,sv))
    xlabel('Pupil @Stim')
    ylabel('Neural baseline')
    
    nexttile
    scatter(cat(1,stim_ppl{:})-cat(1,bsln_ppl{:}), bc(:,sv))
    xlabel('Pupil @Stim - Bsln.')
    ylabel('Neural baseline')
    
    nexttile
    scatter(cat(1,bsln_ppl{:}), dc(:,sv))
    xlabel('Pupil baseline')
    ylabel('Neural dev. from mean (/wo baseline)')
    
    nexttile
    scatter(cat(1,stim_ppl{:}), dc(:,sv))
    xlabel('Pupil @Stim')
    ylabel('Neural dev. from mean (/wo baseline)')
    
    nexttile
    scatter(cat(1,stim_ppl{:})-cat(1,bsln_ppl{:}), dc(:,sv))
    xlabel('Pupil @Stim - Bsln.')
    ylabel('Neural dev. from mean (/wo baseline)')
    
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
end