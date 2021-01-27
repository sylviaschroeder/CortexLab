function data = getTaskInfo(folder)

data.choice = readNPY(fullfile(folder, '_ss_trials.choice.npy'));
data.contrastLeft = readNPY(fullfile(folder, '_ss_trials.contrastLeft.npy'));
data.contrastRight = readNPY(fullfile(folder, '_ss_trials.contrastRight.npy'));
data.feedbackTimes = readNPY(fullfile(folder, '_ss_trials.feedback_times.npy'));
data.feedbackType = readNPY(fullfile(folder, '_ss_trials.feedbackType.npy'));
data.goCueTimes = readNPY(fullfile(folder, '_ss_trials.goCue_times.npy'));
data.goCueFreq = readNPY(fullfile(folder, '_ss_trials.goCueFreq.npy'));
data.interactiveDelayIntv = readNPY(fullfile(folder, '_ss_trials.interactiveDelay_intervals.npy'));
data.interactiveDelays = readNPY(fullfile(folder, '_ss_trials.interactiveDelays.npy'));
data.interTrialDelay = readNPY(fullfile(folder, '_ss_trials.interTrialDelay.npy'));
data.negFbPeriod = readNPY(fullfile(folder, '_ss_trials.negativeFeedbackPeriod.npy'));
data.posFbPeriod = readNPY(fullfile(folder, '_ss_trials.positiveFeedbackPeriod.npy'));
data.preStimDelayIntv = readNPY(fullfile(folder, '_ss_trials.preStimDelay_intervals.npy'));
data.preStimDelays = readNPY(fullfile(folder, '_ss_trials.preStimDelays.npy'));
data.repNum = readNPY(fullfile(folder, '_ss_trials.repNum.npy'));
data.responseWin = readNPY(fullfile(folder, '_ss_trials.responseWindow.npy'));
data.rewardVolume = readNPY(fullfile(folder, '_ss_trials.rewardVolume.npy'));
data.stimAltitude = readNPY(fullfile(folder, '_ss_trials.stimAltitude.npy'));
data.stimAzimuth = readNPY(fullfile(folder, '_ss_trials.stimAzimuth.npy'));
data.stimOnIntv = readNPY(fullfile(folder, '_ss_trials.stimOn_intervals.npy'));
data.stimOrientation = readNPY(fullfile(folder, '_ss_trials.stimOrientation.npy'));
data.stimSigma = readNPY(fullfile(folder, '_ss_trials.stimSigma.npy'));
data.stimSpatFreq = readNPY(fullfile(folder, '_ss_trials.stimSpatFreq.npy'));
data.stimTargetThresh = readNPY(fullfile(folder, '_ss_trials.stimTargetThreshold.npy'));

data.extraReward = readNPY(fullfile(folder, '_ss_extraReward.times.npy'));

data.feedbackTimes_approx = readNPY(fullfile(folder, '_ss_trials.approxFeedback_times.npy'));
data.goCueTimes_approx = readNPY(fullfile(folder, '_ss_trials.approxGoCue_times.npy'));
data.stimTimes_approx = readNPY(fullfile(folder, '_ss_trials.approxStimOn_intervals.npy'));
