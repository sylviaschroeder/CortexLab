timelineFolder = '\\ZSERVER.cortexlab.net\Data\expInfo';

db_wheelTask;

%% Load data
k = 1;
clear pars
for iSet = 1:length(db)
    for iExp = 1:length(db(iSet).exp)
        folder = fullfile(timelineFolder, db(iSet).subject, db(iSet).date, ...
            num2str(db(iSet).exp(iExp)));
        % block: stim times, conditions, beep times, choice times
        file = dir(fullfile(folder, '*_Block.mat'));
        data = load(fullfile(folder, file.name));
        block = data.block;
        pars.preStimPeriod(k,1:2) = block.parameters.preStimQuiescentPeriod;
        pars.interactiveDelay(k,1:2) = block.parameters.cueInteractiveDelay;
        pars.respWin(k) = block.parameters.responseWindow;
        pars.posPeriod(k) = block.parameters.positiveFeedbackPeriod;
        pars.negPeriod(k) = block.parameters.negativeFeedbackPeriod;
        pars.targetWidth(k) = block.parameters.targetWidth;
        pars.targetSigma(k,1:2) = block.parameters.cueSigma;
        pars.targetOri(k) = block.parameters.targetOrientation;
        pars.targetAltitude(k) = block.parameters.targetAltitude;
        pars.targetThresh(k) = block.parameters.targetThreshold;
        pars.targetDistances(k) = block.parameters.distBetweenTargets;
        k=k+1;
    end
end