%% load data of animals SS061 and SS060
data1 = load('C:\DATA\electrophys\SS061\20160511\data.mat');
data2 = load('C:\DATA\electrophys\SS060\20160520\data.mat');

%% Parameters
windowStd = 0.5; % in sec; STD of Gaussian window to convolve spike rates
runningThresh = 100;
minDurPerTrial = 0.4;

%% Correlation with running (continuous traces)
% SS061
% (1) convolve neural traces with Gaussian
fs = 1/median(diff(data1.time));
winSize = round(windowStd * fs);
win = pdf('Normal', -5*winSize:5*winSize, 0, winSize);
spikeRates = zeros(size(data1.spikeRates));
for k = 1:size(spikeRates,2)
    spikeRates(:,k) = conv(data1.spikeRates(:,k), win, 'same');
end
t = data1.time;
onsets = data1.stimTimes.onset - t(1);
t = t-t(1);
results1 = nonVis.getCorrToNonVisData(spikeRates, t, ...
    data1.running', t, 'Running speed', data1.cellIDs, -1, 1, onsets);
% SS060
% (1) convolve neural traces with Gaussian
fs = 1/median(diff(data2.time));
winSize = round(windowStd * fs);
win = pdf('Normal', -5*winSize:5*winSize, 0, winSize);
spikeRates = zeros(size(data2.spikeRates));
for k = 1:size(spikeRates,2)
    spikeRates(:,k) = conv(data2.spikeRates(:,k), win, 'same');
end
% only use data up to 2300 s after first measured time point
t = data2.time;
onsets = data2.stimTimes.onset - t(1);
t = t-t(1);
ind = find(t > 2300, 1);
results2 = nonVis.getCorrToNonVisData(spikeRates(1:ind,:), t(1:ind), ...
    data2.running(1:ind)', t(1:ind), 'Running speed', data2.cellIDs, -1, 1, onsets);

% CONCLUSION: data of SS060 look very shakey (neurons appearing and
% disappearing throughout the experiment) and running correlations are
% extremely low (max at 0.07, as compared to 0.7 for SS061)
% --> abandon this dataset!

% Plot traces of most and least correlated neurons of SS061
% num neurons for each condition
n = 10;

fs = 1/median(diff(data1.time));
winSize = round(windowStd * fs);
win = pdf('Normal', -5*winSize:5*winSize, 0, winSize);
spikeRates = zeros(round(size(data1.spikeRates,1)/27),size(data1.spikeRates,2));
running = data1.running;
% convolve spike rates and downsample
for k = 1:size(spikeRates,2)
    sr = conv(data1.spikeRates(:,k), win, 'same');
    sr = decimate(sr, 3);
    sr = decimate(sr, 9);
    spikeRates(:,k) = sr;
end
spikeRates(:,sum(spikeRates,1)==0) = [];
running = decimate(running, 3);
running = decimate(running, 9);
t = data1.time;
t = t(13:27:end);
[rhos,pVals] = corr(spikeRates, running);
[rMost,mostCorr] = sort(rhos,'descend');
[rLeast,leastCorr] = sort(rhos,'ascend');
% 8th least corr. neuron has disappearing trace -> don't use
rLeast(8) = [];
leastCorr(8) = [];
dist = 5;
t = t-t(1);
figure
subplot(9,1,1)
plot(t,running,'r')
ylabel('Running')
axis tight
set(gca,'box','off','FontName','Arial','XTickLabel',[])
title('Correlation with running (gratings and optogenetics)')
subplot(9,1,2:9)
hold on
y = 0;
for j = 1:n
    tr = spikeRates(:,mostCorr(j));
    tr = (tr-mean(tr))/std(tr);
    plot(t,tr+y, 'k')
    y = y-dist;
end
y = y - dist;
for j = 1:n
    tr = spikeRates(:,leastCorr(j));
    tr = (tr-mean(tr))/std(tr);
    plot(t,tr+y, 'k')
    y = y-dist;
end
axis tight
xlabel('Time (s)')
ylabel('Firing rate (normalised)')
set(gca,'FontName','Arial','YTick',[])
% plot histogram of rhos
bins = -.5:.1:.8;
n1 = hist(rhos(pVals<0.05), bins);
n2 = hist(rhos(pVals>=0.05), bins);
figure
bar(bins, [n2', n1'], 'stacked');
colormap([1 1 1;0 0 0])
xlim(bins([1 end]))
ylim([0 max(n1+n2)*1.05])
legend('p >= 0.05', 'p < 0.05')
set(gca,'box','off')
xlabel('Correlation coeff.')
ylabel('# Neurons')
title('Gratings and V1 inactivation')

%% Compute tuning curves and linear fits for running and not running
fixedPars = [1 5];
running = double(data1.running > runningThresh);
running = ssLocal.getTracesPerStimulus(running, data1.stimMatrix, [1 0]); % [stimuli x repetitions]
running = squeeze(mean(running, 4));
running = double(running >= minDurPerTrial) + 1;
responses = ssLocal.getTracesPerStimulus(data1.spikeRates, data1.stimMatrix, [0 0]);
fs = 1/median(diff(data1.time));
stimDur = mean(sum(data1.stimMatrix,2) ./ size(responses,3)) / fs;
responses = sum(responses,4) ./ stimDur; % in Hz, [neuron x stim x rep]
blankResponses = responses(:,data1.blanks,:);
responses = responses(:,data1.directions(:,2),:);
% conditions: 1: not running + no laser, 2: running + no laser, 3: not
% running + laser on, 4: running + laser on
conditions = running + data1.laserOn * 2;
condStim = conditions(data1.directions(:,2),:);
stim = repmat(data1.directions(:,1),1,size(responses,3));
parameters = zeros(size(responses,1), 4, 5);
intercepts = NaN(size(responses,1),2);
slopes = NaN(size(responses,1),2);
maxRespDiffs = NaN(size(responses,1),2);
modulations = cell(size(responses,1),2);
for iCell = 1:size(responses,1)
    fprintf('%d ',iCell)
    resp = squeeze(responses(iCell,:,:));
    if sum(resp(:)>0) < 2
        continue
    end
    
    parsInit = fitoriWrapped(stim(:), resp(:));
    pars = NaN(4, length(parsInit));
    for cond = 1:4
        if all(resp(condStim==cond) == 0)
            continue
        end
        pars(cond,:) = fitoriWrapped(stim(condStim==cond), ...
            resp(condStim==cond), [], [parsInit(1) NaN NaN NaN parsInit(5)]);
    end
    parameters(iCell,:,:) = pars;
    
    for c = 1:2
        c1 = c*2-1;
        c2 = c*2;
        r1 = NaN(size(resp));
        r1(condStim==c1) = resp(condStim==c1);
        r2 = NaN(size(resp));
        r2(condStim==c2) = resp(condStim==c2);
        r = nansum(cat(3,r1,r2),3);
        ind = ~all(isnan([r1,r2]),2);
        s = stim(ind,:);
        run = condStim(ind,:);
        r = r(ind,:);
        r1 = r1(ind,:);
        r2 = r2(ind,:);
        pVals = anovan(r(:),[s(:),run(:)],'model','interaction','display','off');
        
        resp1 = nanmean(r1,2);
        resp2 = nanmean(r2,2);
        maxi = max([resp1; resp2]);
        resp1 = resp1 ./ maxi;
        resp2 = resp2 ./ maxi;
        if all(pVals >= 0.05) % neuron not modulated by stimulus or nonvisual signal
            modulations{iCell,c} = 'none';
        elseif pVals(2) < 0.05 && pVals(3) >= 0.05 % purely additive
            [y,gof] = fit(resp1, resp2, @(a,x) x+a, ...
                'StartPoint', median(resp1)-median(resp2));
            intercepts(iCell,c) = y.a;
            slopes(iCell,c) = 1;
            maxRespDiffs(iCell,c) = y.a;
            modulations{iCell,c} = 'additive';
        elseif pVals(2) >= 0.05 && pVals(3) < 0.05 % purely multiplicative
            meanResp = nanmean([resp1; resp2]);
            resp1 = resp1 - meanResp; % subtract mean response so that no intercept is necessary when fitting
            resp2 = resp2 - meanResp;
            [y,gof] = fit(resp1, resp2, @(a,x) x.*a, ...
                'StartPoint', nanmean(resp1./resp2));
            if gof.adjrsquare < 0.3 % check whether the interaction is multiplicative
                modulations{iCell,c} = 'none';
            else
                intercept = meanResp - y.a * meanResp; % shift fitted line by meanResp on x and y axis
                % calculate response difference;
                intercepts(iCell,c) = intercept;
                slopes(iCell,c) = y.a;
                % calculate response difference
                if abs(y.a) <= 1
                    % if the absolute slope is smaller one, then the response
                    % to the most driving stimulus during small pupil is 1
                    respDiff = y.a+intercept - 1;
                else
                    % if the absolute slope is larger one, then the response
                    % to the most driving stimulus during large pupil is 1
                    respDiff = 1 - (1-intercept)/y.a;
                end
                maxRespDiffs(iCell,c) = respDiff;
                modulations{iCell,c} = 'multiplicative';
            end
        elseif pVals(2) < 0.05 && pVals(3) < 0.05 % mixed additive and muliplicative
            y = fitlm(resp1,resp2);
            coefficients = y.Coefficients.Estimate;
            intercepts(iCell,c) = coefficients(1);
            slopes(iCell,c) = coefficients(2);
            if abs(coefficients(2)) <= 1
                respDiff = coefficients(2)+coefficients(1) - 1;
            else
                respDiff = 1 - (1-coefficients(1))/coefficients(2);
            end
            maxRespDiffs(iCell,c) = respDiff;
            modulations{iCell,c} = 'mixed';
        else % only modulated by stimulus
            intercepts(iCell,c) = 0;
            slopes(iCell,c) = 1;
            maxRespDiffs(iCell,c) = 0;
            modulations{iCell,c} = 'stimOnly';
        end
    end
end
fprintf('\n')

%% Plot tuning curves
x = 0:360;
cols = 'kr';
lin = {'-','--'};
dirs = data1.directions(:,1);
condBlank = conditions(data1.blanks,:);
cols2 = 'kc';
for iCell = 1:size(responses,1)
    resp = squeeze(responses(iCell,:,:));
    br = squeeze(blankResponses(iCell,:,:));
    figure('position',[330 655 1530 420])
    subplot(1,3,1)
    hold on
    subplot(1,3,2)
    hold on
    h = zeros(1,4);
    maxi = 0;
    mini = 0;
    for c = 1:4
        col = mod(c,2);
        if col == 0
            col = 2;
        end
        col = cols(col);
        subplot(1,3,ceil(c/2))
        r = NaN(size(resp));
        r(condStim==c) = resp(condStim==c);
        m = nanmean(r,2);
        s = nanstd(r,0,2) ./ sqrt(sum(condStim==c,2));
        errorbar(dirs, m, s, ['o' col]);
        y = orituneWrapped(squeeze(parameters(iCell,c,:)), x);
        h(c) = plot(x,y,[col,lin{ceil(c/2)}],'LineWidth',2);
        if max(m+s) > maxi
            maxi = max(m+s);
        end
        if min(m-s) < mini
            mini = min(m-s);
        end
        r = NaN(size(br));
        r(condBlank==c) = br(condBlank==c);
        m = nanmean(r(:));
        s = nanstd(r(:)) / sqrt(sum(condBlank(:)==c));
        fill(x([1 end end 1]), [[1 1].*(m+s),[1 1].*(m-s)], 'k', ...
            'FaceColor',col,'EdgeColor','none','FaceAlpha',.3)
        plot(x([1 end]),[m m],[col lin{ceil(c/2)}],'LineWidth',1)
        if max(m+s) > maxi
            maxi = max(m+s);
        end
        if min(m-s) < mini
            mini = min(m-s);
        end
    end
    rng = maxi-mini;
    maxi = maxi + 0.05*rng;
    mini = mini - 0.05*rng;
    for k = 1:2
        subplot(1,3,k)
        xlim([-5 365])
        ylim([mini maxi])
        xlabel('Direction (degrees)')
        ylabel('Firing rate (Hz)')
    end
    subplot(1,3,1)
    legend(h(1:2), {'not running','running'})
    title('V1 intact')
    subplot(1,3,2)
    legend(h(3:4), {'not running','running'})
    title('V1 inactivated')
    
    % plot not running vs running
    subplot(1,3,3)
    hold on
    h = zeros(1,2);
    for c = 1:2
        c1 = c*2-1;
        c2 = c*2;
        r1 = NaN(size(resp));
        r1(condStim==c1) = resp(condStim==c1);
        r2 = NaN(size(resp));
        r2(condStim==c2) = resp(condStim==c2);
        ind = ~all(isnan([r1,r2]),2);
        resp1 = nanmean(r1(ind,:),2);
        resp2 = nanmean(r2(ind,:),2);
        maxi = max([resp1; resp2]);
        resp1 = resp1 ./ maxi;
        resp2 = resp2 ./ maxi;
        h(c) = plot(resp1, resp2, ['o' cols2(c)]);
        mini = min([0; resp1; resp2]);
        maxi = 1;
        rng = maxi-mini;
        mini = mini-.05*rng;
        maxi = maxi+.05*rng;
        if ~isnan(intercepts(iCell,c))
            plot([mini maxi],[mini maxi].*slopes(iCell,c)+intercepts(iCell,c), ...
                'Color',cols2(c),'LineWidth',2)
        end
    end
    plot([mini maxi],[mini maxi],'k:')
    axis([mini maxi mini maxi])
    axis square
    xlabel('not running')
    ylabel('running')
    legend(h,{['V1 intact (' modulations{iCell,1} ')'], ...
        ['V1 inactivated (' modulations{iCell,2} ')']})
    
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(fullfile('C:\RESULTS\Electrophys_SC\nonVisualEffects\plots_running\tuningCurves', ...
        sprintf('tuningCurve%03d.jpg', iCell)), '-djpeg','-r0')
    close gcf
end

%% Plot max response differences between running vs not running
maxi = 1.6;
col = lines(1);
figure
plot(maxRespDiffs(:,1),maxRespDiffs(:,2),'o','color',col)
hold on
plot([-maxi maxi],[-maxi maxi],'k:')
plot([-maxi maxi],[0 0],'k:')
plot([0 0],[-maxi maxi],'k:')
ind=~any(isnan(maxRespDiffs),2);
p = signtest(maxRespDiffs(ind,1),maxRespDiffs(ind,2));
m = median(maxRespDiffs(ind,1)-maxRespDiffs(ind,2));
plot([-maxi maxi],[-maxi maxi]-m, 'color',col)
axis([-maxi maxi -maxi maxi])
axis square
set(gca,'box','off')
title(sprintf('Response diff. at pref. stim.: no running - running (n = %d, p = %.3f)',sum(ind), p))
xlabel('when V1 intact')
ylabel('when V1 inactivated')