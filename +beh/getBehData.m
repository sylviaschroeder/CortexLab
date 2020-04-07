function psychoDat = getBehData(allBlocks)
%  
% allBlocks         Block structure or cell array of structs
% figDestination    String; Path where fig is to be saved or 'none'
% plotFit           Logical; If true, a plot of the fitted data is made
%
% psychoDat         Struct containing expRef, contrasts, response made,
%                   feedback, repeat numbers and response times

%% Initialize & collect data
if isa(allBlocks, 'struct')
    allBlocks = {allBlocks};
elseif ~isa(allBlocks, 'struct')&&~isa(allBlocks, 'cell')
    error('allBlocks must be a single block structure or a cell array of blocks');
end

numCompletedTrials = 0;
for b = 1:length(allBlocks)
    numCompletedTrials = numCompletedTrials + allBlocks{b}.numCompletedTrials;
end

if numCompletedTrials == 0
    psychoDat = struct(...
        'expRef',{''},...
        'contrast',[],...
        'resp',[],...
        'feedback',[],...
        'repeatNum',[],...
        'trialStartTimes',[],...
        'trialStimTimes',[],...
        'trialResponseTimes',[],...
        'reactionTime',[],...
        'rewardSize',[], ...
        'respWindow', []);
    return
end

expRef = cell(1,length(allBlocks));
contrast = zeros(2,numCompletedTrials);
resp = zeros(1,numCompletedTrials);
feedback = zeros(1,numCompletedTrials);
repeatNum = zeros(1,numCompletedTrials);
trialInteractiveTimes = zeros(1,numCompletedTrials); % start of interactive period relative to stimulus onset
trialResponseTimes = zeros(1,numCompletedTrials); % time of response made relative to stimulus onset
rewardSize = zeros(1,numCompletedTrials);
respWindow = inf(1,length(allBlocks));
tInd = 1;
for b = 1:length(allBlocks)
    fOneInd = 1;
    expRef{b} = allBlocks{b}.expRef;
    respWindow(b) = allBlocks{b}.parameters.responseWindow;
%     resp = [allBlocks{b}.trial.responseMadeID];
    for t = 1:allBlocks{b}.numCompletedTrials
        contrast(:,tInd) = allBlocks{b}.trial(t).condition.visCueContrast;
        resp(tInd) = allBlocks{b}.trial(t).responseMadeID; 
        feedback(tInd) = allBlocks{b}.trial(t).feedbackType;
        repeatNum(tInd) = allBlocks{b}.trial(t).condition.repeatNum;
        trialInteractiveTimes(tInd) = allBlocks{b}.trial(t).interactiveStartedTime - ...
            allBlocks{b}.trial(t).stimulusCueStartedTime;
        trialResponseTimes(tInd) = allBlocks{b}.trial(t).responseMadeTime - ...
            allBlocks{b}.trial(t).stimulusCueStartedTime;
        if allBlocks{b}.trial(t).feedbackType==1
            rewardSize(tInd) = allBlocks{b}.rewardDeliveredSizes(fOneInd);
            fOneInd = fOneInd+1;
        end
        tInd = tInd+1;
    end
end

psychoDat.expRef = expRef;
psychoDat.contrast = contrast;
psychoDat.resp = resp;
psychoDat.feedback = feedback;
psychoDat.repeatNum = repeatNum;
psychoDat.trialInteractiveTimes = trialInteractiveTimes;
psychoDat.trialResponseTimes = trialResponseTimes;
psychoDat.rewardSize = rewardSize;
psychoDat.respWindow = respWindow;