function plotPsychometric(psychoDat, figDestination)
% generatePsychometric plots the psychometic data from the 
% blocks entered and optionally fits the data with phychoFit or GLM
% (for 2AFC and 2AUC respectively), as well as saving the figure to the
% path 'figDestination'.
%  
% allBlocks         Block structure or cell array of structs
% figDestination    String; Path where fig is to be saved or 'none'
% plotFit           Logical; If true, a plot of the fitted data is made
%
% psychoDat         Struct containing expRef, contrasts, response made,
%                   feedback, repeat numbers and response times

%% Calculate psychometric data
% rwdTypes = unique(rewardSize(rewardSize>0)); %Array of unique reward types
% rwd = rewardSize';
% perf = (perf / histc(repeatNum,1))*100; %Performance = number of correct / frequency of repeatNum 1 * 100
perf = sum(feedback==1&repeatNum==1)/sum(repeatNum==1)*100;
respTypes = unique(resp(resp>0));
numRespTypes = numel(respTypes);

cDiff = diff(contrast);
cVals = unique(cDiff);

psychoM = zeros(numRespTypes,length(cVals));
psychoMCI = zeros(numRespTypes,length(cVals));
meanRTs = zeros(numRespTypes,length(cVals));
meanRTsCIlow = zeros(numRespTypes,length(cVals));
meanRTsCIup = zeros(numRespTypes,length(cVals));
numTrials = zeros(1,length(cVals));
numChooseR = zeros(numRespTypes, length(cVals));
for r = 1:numRespTypes
    for c = 1:length(cVals) %For the number of unique contrasts
%         incl_lsrOn = repeatNum==1&cDiff==cVals(c)&rewardSize(2,:)>0; %Logical array of trials that aren't repeats of each contrast
%         incl_lsrOff = repeatNum==1&cDiff==cVals(c)&rewardSize(1,:)==0; %Logical array of trials that aren't repeats of each contrast
        incl = repeatNum==1&cDiff==cVals(c); %Logical array of trials that aren't repeats of each contrast
        numTrials(c) = sum(incl); %Number of trails is equal to number of non-repeat trials
        numChooseR(r,c) = sum(resp==respTypes(r)&incl);
        
        psychoM(r, c) = numChooseR(r,c)/numTrials(c);
        psychoMCI(r, c) = 1.96*sqrt(psychoM(r, c)*(1-psychoM(r, c))/numTrials(c)); % 1.96 * std-error (binomial distribution)
        
%         meanRTs(r, c) = mean(rt(resp==respTypes(r)&incl));
%         meanRTsCIlow(r, c) = std(rt(resp==respTypes(r)&incl))./sqrt(sum(resp==respTypes(r)&incl));
%         meanRTsCIup(r, c) = std(rt(resp==respTypes(r)&incl))./sqrt(sum(resp==respTypes(r)&incl));
        
        
        q = quantile(rt(resp==respTypes(r)&incl), 3);
        meanRTs(r, c) = q(2);
        meanRTsCIlow(r, c) = q(2)-q(1);
        meanRTsCIup(r, c) = q(3)-q(2);
    end
end

%% Plotting curve
if ~isempty(figDestination)
    colors([2 1],:) = lines(2); % orange (left), blue (right)
    colors(3,:) = [0.2 0.2 0.2];
    labels = {'Choose L','Choose R','NoGo'};
    
    expRefMod = allBlocks{1}.expRef; expRefMod(expRefMod=='_') = '-';
    f(1) = figure('Name', expRefMod, 'NumberTitle', 'Off', 'Position', [680 678 560 420]);
    
    trialNumber = 1:numCompletedTrials; % for reaction time plot x-axis
    
    hold on
    h = zeros(1, numRespTypes);
    for r = 1:numRespTypes
        h(r) = beh.plotWithErr(cVals, psychoM(r,:).*100, psychoMCI(r,:).*100, colors(r,:));
    end
    legend(h, labels)
    ylim([-5 105])
    set(gca, 'XTick', cVals)
    xlabel('Contrast')
    ylabel('Correct choice (%)')
    title([expRefMod ', numTrials = ' num2str(numCompletedTrials) ', perf = ' num2str(perf,3) '%']);
    
    f(2) = figure('Position', [1260 678 560 420]);
    t = .05:.1:min(max(respWindow),5);
    % reaction time to go left
    subplot(2,1,1)
    hold on
    n1 = hist(rt(repeatNum==1 & resp==1 & cDiff<0), t);
    n2 = hist(rt(repeatNum==1 & resp==1 & cDiff>=0), t);
    n1 = n1 ./ sum(n1) .* 100;
    n2 = n2 ./ sum(n2) .* 100;
    m = max([n1 n2]);
    plot(t, n1, 'Color', colors(1,:), 'LineWidth', 2)
    plot(t, n2, ':', 'Color', colors(1,:), 'LineWidth', 2)
    legend('correct', 'incorrect')
    ylabel('Trials (%)')
    title('Reaction time: choose Left')
    % reaction time to go right
    subplot(2,1,2)
    hold on
    n1 = hist(rt(repeatNum==1 & resp==2 & cDiff>0), t);
    n2 = hist(rt(repeatNum==1 & resp==2 & cDiff<=0), t);
    n1 = n1 ./ sum(n1) .* 100;
    n2 = n2 ./ sum(n2) .* 100;
    m = max([m n1 n2]);
    plot(t, n1, 'Color', colors(2,:), 'LineWidth', 2)
    plot(t, n2, ':', 'Color', colors(2,:), 'LineWidth', 2)
    ylim([0 m])
    legend('correct', 'incorrect')
    ylabel('Trials (%)')
    title('Reaction time: choose Right')
    subplot(2,1,1)
    ylim([0 m])
end 

psychoDat.expRef = expRef;
psychoDat.contrast = contrast;
psychoDat.resp = resp;
psychoDat.feedback = feedback;
psychoDat.repeatNum = repeatNum;
psychoDat.rt = rt;
psychoDat.rewardSize = rewardSize;    

% if length(unique(resp))>2&&plotFit==true
%     plotFit = 2;
% elseif length(unique(resp))<3&&plotFit==true
%     plotFit = 1;
% end
% 
% %% Plotting the fit
% switch plotFit
%     case 1 % PsychoFit
%         nn = numTrials;
%         xx = cVals*100; %Unique contrasts in percent
%         if size(psychoM,1)>1
%             pp = psychoM(2,:); %Performance
%         else
%             pp = psychoM;
%         end
% 
%         if any(xx==0)
%             intercept = psychoM(xx==0);
%         else
%             intercept=0;
% %             intercept = interp1(pp,xx, 0.5, 'spline');
%         end
% 
%         %Calculate params for line fit & pass to Psychofit
%         parstart = [intercept mean(nn) min(pp)];
%         parmin = [-40 10 min(pp)];
%         parmax = [40 max(nn) max(pp)];
%         pars = mle_fit_psycho([xx;nn;pp],'erf_psycho',parstart,parmin,parmax,25); % Matteo's function
%     %     pars = mle_fit_psycho([cVals;nn;pp],'erf_psycho_2gammas',parstart,parmin,parmax,1);
% 
%    
%         %Plot the data
%         f(2) = figure('Name', [expRefMod '_fit'], 'NumberTitle', 'Off', 'Units', 'Normalized', 'Position', [0.3 0.6 0.50 0.30]);
%         plot([0 0],[0 1],'k:');
%         hold on
%         plot(xx,pp,'ko','markerfacecolor','k');  % Plot points
%         plot(-100:100, erf_psycho( pars, -100:100 ),'b','LineWidth',1.1); % Plot the fit
%         plot([-100 100], [0.5 0.5], 'k:'); % Extend axis limit
%         xlabel('contrast (%)', 'FontWeight', 'bold');
%         ylabel('% Right', 'FontWeight', 'bold');
%         hold off
%         title(sprintf('b = %2.0f, t = %2.0f, l = %0.2f, n = %2.0f',pars, sum(nn)));
%     case 2 % GLM
%         if exist('GLM','file')==2||exist('GLM','builtin')==2
%             data.contrast_cond = contrast';
%             data.response = resp';
%             data.repeatNum = repeatNum';
%             if any(all(contrast))
%                 modelString = 'C^N';
%             else
%                 modelString = 'C^N-subset';
%             end
%             g = GLM(data);
%             g = g.setModel(modelString);
%             g = g.fit;
%             f(2) = figure('Name', [expRefMod '_' modelString '-fit'], 'NumberTitle', 'Off', 'Units', 'Normalized', 'Position', [0.3 0.6 0.50 0.30]);
%             g.plotData;g.plotFit;
%         else
%             warning('Could not plot fit; GLM toolbox not found');
%         end
% end