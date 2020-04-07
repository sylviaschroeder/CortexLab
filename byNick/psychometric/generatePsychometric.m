

function psychoDat = generatePsychometric(allBlocks, figDestination)


numCompletedTrials = 0;
for b = 1:length(allBlocks)
    
    numCompletedTrials = numCompletedTrials + allBlocks(b).numCompletedTrials;
end


contrast = zeros(1,numCompletedTrials);
resp = zeros(1,numCompletedTrials);
feedback = zeros(1,numCompletedTrials);
repeatNum = zeros(1,numCompletedTrials);
rt = zeros(1,numCompletedTrials);
rewardSize = zeros(1,numCompletedTrials);
tInd = 1;
for b = 1:length(allBlocks)
    fOneInd = 1;
    for t = 1:allBlocks(b).numCompletedTrials
        
        if allBlocks(b).trial(t).condition.visCueContrast(1)>0
            contrast(tInd) = -allBlocks(b).trial(t).condition.visCueContrast(1);
        elseif allBlocks(b).trial(t).condition.visCueContrast(2)>0
            contrast(tInd) = allBlocks(b).trial(t).condition.visCueContrast(2);
        else
            contrast(tInd) = 0;
        end
        
        resp(tInd) = allBlocks(b).trial(t).responseMadeID;
        feedback(tInd) = allBlocks(b).trial(t).feedbackType;
        repeatNum(tInd) = allBlocks(b).trial(t).condition.repeatNum;
        if allBlocks(b).trial(t).feedbackType==1
            rewardSize(tInd) = allBlocks(b).rewardDeliveredSizes(fOneInd); 
            fOneInd = fOneInd+1;
        end
        
        
        rt(tInd) = allBlocks(b).trial(t).responseMadeTime- allBlocks(b).trial(t).interactiveStartedTime;
        tInd = tInd+1;
    end
end

psychoDat.contrast = contrast;
psychoDat.resp = resp;
psychoDat.feedback = feedback;
psychoDat.repeatNum = repeatNum;
psychoDat.rt = rt;
psychoDat.rewardSize = rewardSize;

respTypes = unique(resp(resp>0));
numRespTypes = numel(respTypes);

cVals = unique(contrast);

psychoM = zeros(numRespTypes,length(cVals));
psychoMCI = zeros(numRespTypes,length(cVals));
meanRTs = zeros(numRespTypes,length(cVals));
meanRTsCI = zeros(numRespTypes,length(cVals));
numTrials = zeros(1,length(cVals));
numChooseR = zeros(numRespTypes, length(cVals));
for r = 1:numRespTypes
    for c = 1:length(cVals)
        incl = repeatNum==1&contrast==cVals(c);
        numTrials(c) = sum(incl);
        numChooseR(r,c) = sum(resp==respTypes(r)&incl);
        
        psychoM(r, c) = numChooseR(r,c)/numTrials(c);
        psychoMCI(r, c) = 1.96*sqrt(psychoM(r, c)*(1-psychoM(r, c))/numTrials(c));
        
%         meanRTs(r, c) = mean(rt(resp==respTypes(r)&incl));
%         meanRTsCIlow(r, c) = std(rt(resp==respTypes(r)&incl))./sqrt(sum(resp==respTypes(r)&incl));
%         meanRTsCIup(r, c) = std(rt(resp==respTypes(r)&incl))./sqrt(sum(resp==respTypes(r)&incl));
        
        q = quantile(rt(resp==respTypes(r)&incl), 3);
        meanRTs(r, c) = q(2);
        meanRTsCIlow(r, c) = q(2)-q(1);
        meanRTsCIup(r, c) = q(3)-q(2);
    end
end
if ~isempty(figDestination)
%     colors = hsv(numRespTypes);
    colors(1,:) = [0 0.5 1];
    colors(2,:) = [1 0.5 0]; 
    colors(3,:) = [0.2 0.2 0.2];
    
    f = figure;
    
    leg = {};
    
    for r = 1:numRespTypes
        subplot(3,1,1);
        %     plot(cVals, psychoM(r, :), 'o', 'Color', colors(r,:));
        %     hold on;
        %     plot(reshape([cVals; cVals; nan(size(cVals))],1,3*length(cVals)), ...
        %         reshape([psychoM(r,:)-psychoMCI(r,:); psychoM(r,:)+psychoMCI(r,:); nan(size(psychoM(r,:)))], 1, 3*length(psychoM(r,:))), 'Color', colors(r,:));
        plotWithErr(cVals, psychoM(r,:), psychoMCI(r,:), colors(r,:));
        hold on;
        plot(cVals, psychoM(r,:), 'ko');
        leg{end+1} = num2str(respTypes(r));
        leg{end+1} = '';
        
        subplot(3,1,2);
%         plot(cVals, meanRTs(r,:), 'o', 'Color', colors(r,:));
%         hold on;
%         plot(reshape([cVals; cVals; nan(size(cVals))],1,3*length(cVals)), ...
%             reshape([meanRTs(r,:)-meanRTsCIlow(r,:); meanRTs(r,:)+meanRTsCIup(r,:); nan(size(meanRTs(r,:)))], 1, 3*length(meanRTs(r,:))), 'Color', colors(r,:));
        
        minXDiff = min(diff(cVals));
        for x = 1:length(cVals)

            
            theseRT = psychoDat.rt(psychoDat.resp==1&psychoDat.contrast==cVals(x));
            if ~isempty(theseRT)
%                 distributionPlot(theseRT', 'xValues', cVals(x), 'histOri', 'left', 'widthDiv', [2 1], 'color', colors(1,:), 'histOpt', 0, 'divFactor', [0:0.05:1.5], 'distWidth', minXDiff/3*2)
                distributionPlot(theseRT', 'xValues', cVals(x), 'histOri', 'left', 'widthDiv', [2 1], 'color', colors(1,:), 'histOpt', 0, 'distWidth', minXDiff/3*2)
            end
            hold on;
            theseRT = psychoDat.rt(psychoDat.resp==2&psychoDat.contrast==cVals(x));
            if ~isempty(theseRT)
                distributionPlot(theseRT', 'xValues', cVals(x), 'histOri', 'right', 'widthDiv', [2 2], 'color', colors(2,:), 'histOpt', 0, 'distWidth', minXDiff/3*2)
            end
        end
        
        
        xlim([cVals(1) cVals(end)]*1.1);
        xlabel('contrast');
        ylabel('reaction time (sec)');
        ylim([0 max(psychoDat.rt)]);
        makepretty;
        
    end
    
    subplot(3,1,1);
    
%     psy.plot2ADCwithAlt(obj.PsychometricAxes.Handle, obj.Block);
    
    % legend(leg, 'Location', 'EastOutside');
    plot([cVals(1) cVals(end)], [0.5 0.5], 'k:');
    ylim([0 1]);
    xlim([cVals(1) cVals(end)]*1.1);
    xlabel('contrast');
    ylabel('proportion choose R');
    expRefMod = allBlocks(1).expRef; expRefMod(expRefMod=='_') = '-';
    title([expRefMod ', numTrials = ' num2str(numCompletedTrials)]);
    makepretty;
    
    subplot(3,1,3);
    plot(cVals, numTrials, 'ko');
    xlim([cVals(1) cVals(end)]*1.1);
    ylabel('number of trials')
    xlabel('contrast');
    yl = ylim(); ylim([0 yl(2)]);
    makepretty;
    
    set(f, 'Position', [50         131         648         755]);
    
    if ~strcmp(figDestination, 'none');
        saveFig(f, fullfile(figDestination, [allBlocks(1).expRef '_psychometric']));
    end
end

psychoDat.psychoM = psychoM;
psychoDat.psychoMCI = psychoMCI;
psychoDat.meanRTs = meanRTs;
psychoDat.meanRTsCI = meanRTsCI;
psychoDat.numTrials = numTrials;
psychoDat.numChooseR = numChooseR;