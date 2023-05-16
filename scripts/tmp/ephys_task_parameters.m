folder = 'Z:\UCLData\Ephys_Task\Subjects';
f_subjects = dir(folder);
subjects = {};
dates = {};
positions = [];
contrasts = [];
goStim = [];
orientations = [];
interactiveDels = [];
negPeriods = [];
posPeriods = [];
preStimDelays = [];
respWindows = [];
sigmas = [];
spatFreqs = [];
thresholds = [];
for s = 3:length(f_subjects)
    if strcmp(f_subjects(s).name, 'SS090')
        continue
    end
    f_dates = dir(fullfile(folder,f_subjects(s).name,'20*'));
    for d = 1:length(f_dates)
        subjects{end+1,1} = f_subjects(s).name;
        dates{end+1,1} = f_dates(d).name;
        pos = readNPY(fullfile(folder,f_subjects(s).name, f_dates(d).name, 'alf', '_ss_trials.stimAzimuth.npy'));
        positions(end+1,1) = pos(1);
        left = readNPY(fullfile(folder,f_subjects(s).name, f_dates(d).name, 'alf', '_ss_trials.contrastLeft.npy'));
        right = readNPY(fullfile(folder,f_subjects(s).name, f_dates(d).name, 'alf', '_ss_trials.contrastRight.npy'));
        contrasts = unique([contrasts; left]);
        goStim(end+1,1) = sum(left>0 | right>0) / length(left);
        ori = readNPY(fullfile(folder,f_subjects(s).name, f_dates(d).name, 'alf', '_ss_trials.stimOrientation.npy'));
        orientations(end+1,1) = ori(1);
        dels = readNPY(fullfile(folder,f_subjects(s).name, f_dates(d).name, 'alf', '_ss_trials.interactiveDelay_intervals.npy'));
        interactiveDels(end+1,1:2) = dels(1:2);
        negPer = readNPY(fullfile(folder,f_subjects(s).name, f_dates(d).name, 'alf', '_ss_trials.negativeFeedbackPeriod.npy'));
        negPeriods(end+1,1) = negPer(1);
        posPer = readNPY(fullfile(folder,f_subjects(s).name, f_dates(d).name, 'alf', '_ss_trials.positiveFeedbackPeriod.npy'));
        posPeriods(end+1,1) = posPer(1);
        preStimDel = readNPY(fullfile(folder,f_subjects(s).name, f_dates(d).name, 'alf', '_ss_trials.preStimDelay_intervals.npy'));
        preStimDelays(end+1,1:2) = preStimDel(1,:);
%         reps = readNPY(fullfile(folder,f_subjects(s).name, f_dates(d).name, 'alf', '_ss_trials.repNum.npy'));
        respWin = readNPY(fullfile(folder,f_subjects(s).name, f_dates(d).name, 'alf', '_ss_trials.responseWindow.npy'));
        respWindows(end+1,1) = respWin(1);
        sig = readNPY(fullfile(folder,f_subjects(s).name, f_dates(d).name, 'alf', '_ss_trials.stimSigma.npy'));
        sigmas(end+1,1) = sig(1);
        spatF = readNPY(fullfile(folder,f_subjects(s).name, f_dates(d).name, 'alf', '_ss_trials.stimSpatFreq.npy'));
        spatFreqs(end+1,1) = spatF(1);
        thresh = readNPY(fullfile(folder,f_subjects(s).name, f_dates(d).name, 'alf', '_ss_trials.stimTargetThreshold.npy'));
        thresholds(end+1,1) = thresh(1);
    end
end
