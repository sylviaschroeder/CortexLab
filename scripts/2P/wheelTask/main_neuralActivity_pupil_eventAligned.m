%% Folders
folderData = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\NPY\task_2p';
folderTools = 'C:\STORAGE\workspaces';
folderCode = 'C:\dev\workspace\CortexLab';
folderPlots = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\wheelTask\pupil_neuralActivity';

%% Add paths
addpath(genpath(fullfile(folderTools, 'npy-matlab')))
addpath(fullfile(folderCode))

%% Align data to task events
pre = 2;
post = 5;
gap = 3;
% sortTime = [0 2];
sortTime = [-1 0];

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
            trials = find(all(contrasts == contr_unq(c,:), 2));
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
            % neurons
            algn_ca{c} = NaN(length(trials), pre_ca+post_ca+1, size(caTraces,2));
            pl_unq = unique(planes);
            for pl = pl_unq'
                tm = caTime + planeDelays(pl);
                ind = NaN(length(trials),1);
                for t = 1:length(trials)
                    ind(t) = find(tm >= stimT(trials(t),1), 1, 'first');
                end
                ind = ind + (-pre_ca : post_ca);
                ind2 = ind;
                invalid = ind < 1 | ind > length(tm);
                ind2(invalid) = 1;
                ns = find(planes == pl);
                for n = ns'
                    tr = caTraces(:,n);
                    al = tr(ind2);
                    al(invalid) = NaN;
                    algn_ca{c}(:,:,n) = al;
                end
            end
        end
        
        a = find(timeAlgn_ca >= sortTime(1), 1);
        o = find(timeAlgn_ca >= sortTime(2), 1);
        for n = 1:size(caTraces,2)
            ca = [];
            ppl = [];
            for c = 1:size(contr_unq,1)
                [~,order] = sort(mean(algn_ca{c}(:,a:o,n),2));
                ca = [ca; algn_ca{c}(order,:,n); NaN(gap,size(algn_ca{c},2))];
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