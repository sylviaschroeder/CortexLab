folderResults = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\nonVisualEffects\modelGratingResp\';
nonVisualSignal = 'pupil';

%% Test relation between correlation to pupil and tuning modulation
data = load(fullfile(folderResults, nonVisualSignal, ...
    'corrsDuringGratingsAndGrayScreen.mat'));
corrs = data.corrs;
data = load(fullfile(folderResults, nonVisualSignal, ...
    'tuning_prefDirSigmaDIFixed_lineFit.mat'));
tuning = data.tuning;

rhosGratings = [];
rhosGrayScreen = [];
isTuned = [];
meanStimResp = [];
prefStimResp = [];
nonprefStimResp = []; % can be opposite direction or orthogonal orientation
isGad = [];
for iExp = 1:length(corrs)
    for iPlane = 1:length(corrs(iExp).plane)
        rhosGratings = [rhosGratings; corrs(iExp).plane(iPlane).gratings.rhos];
        rhosGrayScreen = [rhosGrayScreen; ...
            corrs(iExp).plane(iPlane).grayScreen.rhos];
        
        isTuned = [isTuned; tuning(iExp).plane(iPlane).isTuned];
        isGad = [isGad; tuning(iExp).plane(iPlane).isGad];
        for iCell = 1:length(tuning(iExp).plane(iPlane).cond(1).cell)
            means = NaN(1,2);
            prefs = NaN(1,2);
            nonprefs = NaN(1,2);
            for c = 1:2
                m = tuning(iExp).plane(iPlane).cond(c).cell(iCell).responses;
                if isempty(m)
                    break
                end
                pars = tuning(iExp).plane(iPlane).cond(c).cell(iCell).parameters;
                if length(pars) == 1
                    means(c) = pars;
                    continue
                else
                    means(c) = nanmean(nanmean(m,2),1);
                    oris = mod(pars(1) + [0 90 180], 360);
                    r = orituneWrappedConditions(pars, oris);
                    prefs(c) = r(1);
                    if r(1)-r(2)>0
                        nonprefs(c) = min(r(2:3));
                    else % suppressed neurons
                        nonprefs(c) = max(r(2:3));
                    end
                end
            end
            meanStimResp(end+1,:) = means;
            prefStimResp(end+1,:) = prefs;
            nonprefStimResp(end+1,:) = nonprefs;
        end
    end
end
meanStimModulation = diff(meanStimResp,1,2) ./ sum(meanStimResp,2);
prefStimModulation = diff(prefStimResp,1,2) ./ sum(prefStimResp,2);
nonprefStimModulation = diff(nonprefStimResp,1,2) ./ sum(nonprefStimResp,2);
modDepthModulation = prefStimResp - nonprefStimResp;
modDepthModulation = diff(modDepthModulation,1,2) ./ sum(modDepthModulation,2);

% Compare correlations between inhibitory and excitatory neurons
% (1) during gratings
figure
hold on
bins = -.5:.05:.8;
n = hist(rhosGratings(isGad==1), bins);
plot(bins, n./sum(n), 'b', 'LineWidth', 2)
n = hist(rhosGratings(isGad==-1), bins);
plot(bins, n./sum(n), 'r', 'LineWidth', 2)
legend('inhibitory', 'excitatory')
xlim(bins([1 end]))
xlabel('Correlation coeff.')
ylabel('Proportion of neurons')
title('Correlation with pupil during gratings')
% during gray screen
figure
hold on
bins = -.5:.05:.8;
n = hist(rhosGrayScreen(isGad==1), bins);
plot(bins, n./sum(n), 'b', 'LineWidth', 2)
n = hist(rhosGrayScreen(isGad==-1), bins);
plot(bins, n./sum(n), 'r', 'LineWidth', 2)
legend('inhibitory', 'excitatory')
xlim(bins([1 end]))
xlabel('Correlation coeff.')
ylabel('Proportion of neurons')
title('Correlation with pupil during gray screen')

% Compare correlations between non-responsive, responsive but untuned, and
% tuned neurons
% (1) during gratings
figure
hold on
bins = -.5:.05:.8;
n = hist(rhosGratings(isnan(isTuned)), bins);
plot(bins, n./sum(n), 'LineWidth', 2)
n = hist(rhosGratings(isTuned==0), bins);
plot(bins, n./sum(n), 'LineWidth', 2)
n = hist(rhosGratings(isTuned==1), bins);
plot(bins ,n./sum(n), 'LineWidth', 2)
legend('non-responsive','untuned','tuned')
xlim(bins([1 end]))
xlabel('Correlation coeff.')
ylabel('Proportion of neurons')
title('Correlation with pupil during gratings')
% (2) during gray screen
figure
hold on
bins = -.5:.05:.8;
n = hist(rhosGrayScreen(isnan(isTuned)), bins);
plot(bins, n./sum(n), 'LineWidth', 2)
n = hist(rhosGrayScreen(isTuned==0), bins);
plot(bins, n./sum(n), 'LineWidth', 2)
n = hist(rhosGrayScreen(isTuned==1), bins);
plot(bins ,n./sum(n), 'LineWidth', 2)
legend('non-responsive','untuned','tuned')
xlim(bins([1 end]))
xlabel('Correlation coeff.')
ylabel('Proportion of neurons')
title('Correlation with pupil during gray screen')

% Compare correlations to tuning modulation (pupil, gratings)
% (1) mean stimulus responses
figure
hold on
ind = isTuned == 0 & prod(sign(meanStimResp),2) > 0;
r = rhosGratings(ind);
m = meanStimModulation(ind);
rhos = r;
mods = m;
plot(r, m, '.k')
ind = isTuned == 1 & prod(sign(meanStimResp),2) > 0;
r = rhosGratings(ind);
m = meanStimModulation(ind);
rhos = [rhos; r];
mods = [mods; m];
plot(r, m, '.r')
ind = ~any(isnan([rhos, mods]),2);
[r,p] = corr(rhos(ind), mods(ind));
legend('untuned','tuned','Location','NorthWest')
xlabel('Correlations coeff.')
ylabel('Mean stim. resp. (large-small)/(large+small)')
title(sprintf('Effects of pupil during gratings (n=%d,rho=%.2f,p=%.4f)',sum(ind),r,p))
% (2) preferred stimulus responses
figure
hold on
ind = ~any(isnan(prefStimResp),2) & prod(sign(prefStimResp),2) > 0;
r = rhosGratings(ind);
m = prefStimModulation(ind);
plot(r, m, '.k')
ind = ~any(isnan([r m]),2);
[r,p] = corr(r(ind), m(ind));
xlabel('Correlations coeff.')
ylabel('Pref. stim. resp. (large-small)/(large+small)')
title(sprintf('Effects of pupil during gratings (n=%d,rho=%.2f,p=%.4f)',sum(ind),r,p))
% (3) nonpreferred stimulus responses
figure
hold on
ind = ~any(isnan(nonprefStimResp),2) & prod(sign(nonprefStimResp),2) > 0;
r = rhosGratings(ind);
m = nonprefStimModulation(ind);
plot(r, m, '.k')
ind = ~any(isnan([r m]),2);
[r,p] = corr(r(ind), m(ind));
xlabel('Correlations coeff.')
ylabel('Nonpref. stim. resp. (large-small)/(large+small)')
title(sprintf('Effects of pupil during gratings (n=%d,rho=%.2f,p=%.4f)',sum(ind),r,p))
% (4) modulation depth
figure
hold on
modDepths = prefStimResp - nonprefStimResp;
ind = ~any(isnan(modDepths),2) & prod(sign(modDepths),2) > 0;
r = rhosGratings(ind);
m = modDepthModulation(ind);
plot(r, m, '.k')
ind = ~any(isnan([r m]),2);
[r,p] = corr(r(ind), m(ind));
xlabel('Correlations coeff.')
ylabel('Modulation depth (large-small)/(large+small)')
title(sprintf('Effects of pupil during gratings (n=%d,rho=%.2f,p=%.4f)',sum(ind),r,p))
