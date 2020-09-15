%% Define data
label = 'neurons';
% label = 'boutons';

% good example dataset
% (1) SC neurons
% subject = 'M150610_SS047'; % (dataset 7)
% date = '2015-11-23';
% expGratings = 1;
% expGray = 5;
% planes = 2:4;
% (2) Boutons (for gratings and gray screen)
% subject = 'M160923_SS070'; % (dataset 5)
% date = '2016-10-18';
% expGratings = 1;
% expGray = 2;
% planes = 9:15;
% (3) Boutons (for darkness)
% subject = 'M170821_SS075'; (dataset 7)
% date = '2017-09-13';
% expGratings = 3;
% expDark = 1;
% planes = 5;
% (4) Boutons (for gratings and darkness)
% subject = 'SS078'; (dataset 16)
% date = '2017-09-28';
% expGratings = 2;
% expDark = 1;
% planes = 5;
if strcmp(label, 'neurons')
%     exSet = 2;
%     exUnits = [];
    exSet = 7;
    exUnits = [1 153; 1 171];
% OLD:
% exUnits = [1 127; 1 190; 1 153; 1 160; 1 41; ...
%     1 32; 1 171; 1 21; 1 214; 1 90]; %[1 149; 1 171]; %[1 81; 1 195];
% exSet = 2;
% exUnits = [1 56; 1 47];
else
    % exSet = 5;
    % % exUnits = [3 3; 3 43];
    % exUnits = [3 3; 5 8];
    
    % exSet = 7;
    % exUnits = [1 268; 1 325];
    % OLD:
    % % exUnits = [1 190; 1 294];
    % % exUnits = [1 333; 1 316; 1 139; 1 268; 1 271; 1 251; 1 132; 1 141; 1 16; 1 260; ...
    % %            1 56; 1 189; 1 267; 1 325; 1 292; 1 73; 1 298; 1 103; 1 122; 1 295];

    exSet = 16;
    exUnits = [1 108; 1 215];
end
numNeurons = size(exUnits,1);

%% Parameters

sigma = 1;
fields = {'expGratings','expGrayScreen','expDark'};
stimuli = {'gratings', 'grayScreen', 'dark'};
signals = {'pupil', 'running'};

xLimits = [-0.65 0.85];
% xLimits = [-0.55 0.55];

%% Folders
% file location depending on PCs
folderPC = 'C:\STORAGE\OneDrive - University College London';
% folderPC = 'C:\Users\skgtchr\OneDrive - University College London'; % ZEN notebook
% folderPC = 'C:\Storage\OneDrive - University College London'; % Lenovo laptop
% data folders
folderROIData = fullfile(folderPC, 'Lab\DATA\InfoStructs');
if strcmp(label, 'neurons')
    folderResults = fullfile(folderPC, 'Lab\RESULTS\nonvisualEffects\modelGratingResp');
    folderRF = fullfile(folderPC, 'Lab\RESULTS\receptiveFields\SC neurons');
    db_driftingGratings_blank
    
    corrections = [];
else
    folderResults = fullfile(folderPC, 'Lab\RESULTS\boutons\nonvisualEffects');
    folderRF = fullfile(folderPC, 'Lab\RESULTS\receptiveFields\boutons');
    db_boutons_driftingGratings_blanks
    
    data = load(fullfile(folderROIData, 'corrections_boutons.mat'));
    corrections = data.corrections;
    doCorrect = data.doCorrect;
end

data = load(fullfile(folderResults,'running', ...
    sprintf('corrsDuringGratingsAndGrayScreen_sigma%.2f.mat',sigma)));
corrsRunning = data.corrs;
running = data.nonVisual;
data = load(fullfile(folderResults,'pupil', ...
    sprintf('corrsDuringGratingsAndGrayScreen_sigma%.2f.mat',sigma)));
corrsPupil = data.corrs;
pupil = data.nonVisual;
% nonvisTime = data.time; %(needed only for comparison: running vs pupil
% diameter)

%% (1) Plots showing traces of behaviour (pupil+running) and most corr. neurons
% cellsMostCorr, cellsLeastCorr contain cellIDs ([plane ID]) of most
% correlated neurons

% parameters
% preprocess neural data
smoothStd = 0.25; %in sec
% num neurons for each condition
% numNeurons = 7;
minSamples = 20;
runningThreshold = 5;

refBeh = 'pupil'; % running or pupil
refStim = 'grayScreen'; % gratings, grayScreen or dark
 
% find most and least correlated neurons
if strcmp(refBeh, 'pupil')
    corrs = corrsPupil;
else
    corrs = corrsRunning;
end
cellIDs = [];
rhos = [];
isSignf = [];
for iPlane = 1:length(corrs(exSet).plane)
    n = length(corrs(exSet).plane(iPlane).(refStim).rhos);
    cellIDs = [cellIDs; [ones(n,1).*db(exSet).planes(iPlane), (1:n)']];
    r = corrs(exSet).plane(iPlane).(refStim).rhos;
    rhos = [rhos;r];
    p = prctile(corrs(exSet).plane(iPlane).(refStim).nullRhos, ...
        [2.5 97.5], 2);
    isSignf = [isSignf; r < p(:,1) | r > p(:,2)];
end
rhos(isSignf<1) = NaN;
[rhosSorted, sortedIDs] = sort(rhos,'descend');
ind = isnan(rhosSorted);
rhosSorted(ind) = [];
sortedIDs(ind) = [];
if length(rhosSorted) < 20
    display('No. of significantly correlated units < 20.')
    return
end
cellsMostCorr = cellIDs(sortedIDs(1:10),:);
cellsLeastCorr = cellIDs(sortedIDs(end:-1:end-9),:);
 
% get traces of dataset
times = cell(length(stimuli), length(db(exSet).planes));
traces = cell(length(stimuli), length(db(exSet).planes));
metas = cell(length(stimuli),1);
for st = 1:length(stimuli)
    if ~isfield(db, fields{st})
        continue
    end
    folder = fullfile(folderROIData, db(exSet).subject, ...
        db(exSet).date, num2str(db(exSet).(fields{st})));
    file = [sprintf('%s_%d_%s', db(exSet).date, ...
        db(exSet).(fields{st}), db(exSet).subject), ...
        '_2P_plane%03d_ROI.mat'];
    for iPlane=1:length(db(exSet).planes)
        if isempty(db(exSet).(fields{st}))
            continue
        end
        % load meta
        data=load(fullfile(folder, sprintf(file,db(exSet).planes(iPlane))));
        meta = data.meta;
        meta.folderTL = strrep(meta.folderTL, 'ioo.ucl.ac.uk', 'cortexlab.net');
        metas{st} = meta;
        if strcmp(stimuli{st},'gratings') && iPlane == 1
            stimTimes = ssLocal.getStimulusResponseInfo(meta);
        end
        times{st,iPlane} = ppbox.getFrameTimes(meta);
        traces{st,iPlane} = metas{st}.F_final;
    end
end
% interpolate neural responses to get a single series of time points
for st = 1:length(stimuli)
    if isempty(times{st,1})
        continue
    end
    for iPlane = 2:length(db(exSet).planes)
        tr = traces{st,iPlane};
        ind = ~all(isnan(tr),1);
        traces{st,iPlane} = interp1(times{st,iPlane}, ...
            tr(:,ind), times{st,1}, 'pchip');
    end
end
times(:,2:end) = [];
% select traces of example cells and convolve them
tracesMost = cell(length(stimuli),1);
tracesLeast = cell(length(stimuli),1);
stdSamples = round(smoothStd / median(diff(times{1,1})));
convWindow = normpdf(-3*stdSamples:3*stdSamples, 0, stdSamples);
for st = 1:length(stimuli)
    if all(cellfun(@isempty, traces(st,:)))
        continue
    end
    tracesMost{st} = zeros(size(traces{st,1},1), numNeurons);
    tracesLeast{st} = zeros(size(traces{st,1},1), numNeurons);
    for n = 1:numNeurons
        pl = cellsMostCorr(n,1) == db(exSet).planes;
        tr = traces{st,pl}(:,cellsMostCorr(n,2));
        tracesMost{st}(:,n) = conv(tr, convWindow, 'same');
        pl = cellsLeastCorr(n,1) == db(exSet).planes;
        tr = traces{st,pl}(:,cellsLeastCorr(n,2));
        tracesLeast{st}(:,n) = conv(tr, convWindow, 'same');
    end
end
    
% load behavioural data
running = cell(length(stimuli),1);
runningTime = cell(length(stimuli),1);
pupil = cell(length(stimuli),1);
pupilTime = cell(length(stimuli),1);
starts = cell(length(stimuli),1);
stops = cell(length(stimuli),1);
for d = 1:length(stimuli)
    if isempty(metas{d})
        continue
    end
    ballData = nonVis.getRunningSpeed(metas{d});
    running{d} = ballData.total / median(diff(ballData.t)) / 53;
    running{d} = interp1(ballData.t, running{d}, times{d});
    running{d} = conv(running{d}, convWindow, 'same');
    
    [p, t] = nonVis.loadPupilData(metas{d});
    if ~isempty(p)
        t(length(p.x)+1:end) = [];
        pupil{d} = nonVis.getPupilDiam(p);
        %     ind = isnan(pupil{d});
        %     indInterp = hist(t(ind), times{d}) > 0;
        warning off
        pupil{d} = interp1(t, pupil{d}, times{d}, 'pchip');
        pupil{d} = conv([ones(1,ceil((length(convWindow)-1)/2)) .* ...
            mean(pupil{d}(1:100)), pupil{d}, ...
            ones(1,floor((length(convWindow)-1)/2)) .* ...
            mean(pupil{d}(end-99:end))], convWindow, 'valid');
        %     pupil{d}(indInterp) = NaN;
    end
    
    if strcmp(refBeh, 'pupil')
        beh = pupil{d};
        ind = pupil{d} > nanmedian(pupil{d});
    else
        beh = running{d};
        ind = running{d} > runningThreshold;
    end
    sta = find(diff(ind)==1);
    sto = find(diff(ind)==-1);
    if sta(1)>sto(1)
        sta = [1, sta];
    end
    if sta(end)>sto(end)
        sto(end+1) = length(beh);
    end
    ind = find(sta(2:end) - sto(1:end-1) < minSamples);
    sta(ind+1) = [];
    sto(ind) = [];
    ind = (sto - sta) >= minSamples;
    starts{d} = sta(ind);
    stops{d} = sto(ind);
end
    
% make figures
% (1) traces of behaviour and traces of most and least correlated neurons
dist = 5;
cols = [.7 0 0; 0 0 .7];
for d = 1:length(stimuli)
    if isempty(traces{d})
        continue
    end
    ax = [];
    figure('Position', [1 41 1920 1083])
    
    subplot(10,1,1)
    if ~isempty(pupil{d})
        hold on
        mini = min(pupil{d});
        maxi = max(pupil{d});
        fill(times{d}([starts{d};stops{d};stops{d};starts{d}]), ...
            [mini mini maxi maxi]', 'k', 'EdgeColor', 'none', 'FaceColor', ...
            [1 1 1].*.9)
        plot(times{d},pupil{d},'k')
        ylabel({'Pupil';'diam.'})
        axis tight
        %     ylim([39 75])
        %     set(gca,'box','off','FontName','Arial','YTick',[],'XTickLabel',[])
        set(gca,'box','off','XTick',[])
        ax(end+1) = gca;
    end
    title(sprintf('Correlation with behaviour - %s (reference: %s, %s)', ...
        stimuli{d}, refBeh, refStim))
    
    subplot(10,1,2)
    hold on
    mini = min(running{d});
    maxi = max(running{d});
    fill(times{d}([starts{d};stops{d};stops{d};starts{d}]), ...
        [mini mini maxi maxi]', 'k', 'EdgeColor', 'none', 'FaceColor', ...
        [1 1 1].*.9)
    plot(times{d},running{d},'k')
    ylabel({'Running';'(cm/s)'})
    axis tight
    %     ylim([-1.5 75])
    set(gca,'box','off','XTick',[])
    ax(end+1) = gca;
    
    subplot(10,1,3:10)
    hold on
    mini = -(2*numNeurons + .7) * dist;
    maxi = 1.5*dist;
    fill(times{d}([starts{d};stops{d};stops{d};starts{d}]), ...
        [mini mini maxi maxi]', 'k', 'EdgeColor', 'none', 'FaceColor', ...
        [1 1 1].*.9)
    y = 0;
    yTicks = [];
    for j = 1:numNeurons
        tr = tracesMost{d}(:,j);
        tr = (tr-nanmean(tr))/nanstd(tr);
        ind = ~isnan(tr);
        plot(times{d}(ind),tr(ind)+y, 'Color', cols(1,:))
        yTicks = [yTicks y];
        y = y-dist;
    end
    y = y - dist;
    for j = 1:numNeurons
        tr = tracesLeast{d}(:,j);
        tr = (tr-nanmean(tr))/nanstd(tr);
        ind = ~isnan(tr);
        plot(times{d}(ind),tr(ind)+y, 'Color', cols(2,:))
        yTicks = [yTicks y];
        y = y-dist;
    end
    if d == 1 %gratings
        stims = [stimTimes.onset'; stimTimes.offset'];
%         plot(stims, ones(2, size(stims,2)) .* 1.4*dist, 'k')
        for j = 1:size(stims,2)
            fill(stims([1 2 2 1],j), [dist dist 1.4*dist 1.4*dist]', ...
                'k', 'EdgeColor', 'none', 'FaceColor', 'k')
        end
    end
    axis tight
    ylim([mini-.2 maxi])
    xlabel('Time (s)')
    ylabel('\DeltaF/F (normalised)')
    set(gca,'FontName','Arial','YTick',flip(yTicks),'YTickLabel', ...
        cellstr(num2str(flip([cellsMostCorr; cellsLeastCorr],1))))
    ax(end+1) = gca;
    
    linkaxes(ax, 'x')
end

%% (2) Get convolved data

rhosPupil = cell(length(corrsPupil), length(stimuli));
nullRhosPupil = cell(length(corrsPupil), length(stimuli));
rhosRunning = cell(length(corrsRunning), length(stimuli));
nullRhosRunning = cell(length(corrsRunning), length(stimuli));
resp = cell(length(corrsPupil),1);
isGad = cell(length(corrsPupil),1);
for k = 1:length(corrsRunning)
    for iPlane = 1:length(corrsRunning(k).plane)
        for st = 1:length(stimuli)
            if ~isfield(corrsPupil(k).plane(iPlane), stimuli{st})
                continue
            end
            rhosPupil{k,st} = [rhosPupil{k,st}; ...
                corrsPupil(k).plane(iPlane).(stimuli{st}).rhos];
            nullRhosPupil{k,st} = [nullRhosPupil{k,st}; ...
                corrsPupil(k).plane(iPlane).(stimuli{st}).nullRhos];
            rhosRunning{k,st} = [rhosRunning{k,st}; ...
                corrsRunning(k).plane(iPlane).(stimuli{st}).rhos];
            nullRhosRunning{k,st} = [nullRhosRunning{k,st}; ...
                corrsRunning(k).plane(iPlane).(stimuli{st}).nullRhos];
        end
        if isfield(corrsRunning(k).plane, 'responsive')
            if isempty(corrsRunning(k).plane(iPlane).responsive)
                r = NaN(size(corrsRunning(k).plane(iPlane).gratings.rhos));
            else
                r = corrsRunning(k).plane(iPlane).responsive;
            end
            resp{k} = [resp{k}; r];
        end
        if isfield(corrsRunning(k).plane, 'isGad')
            if isempty(corrsRunning(k).plane(iPlane).isGad)
                g = NaN(size(corrsRunning(k).plane(iPlane).gratings.rhos));
            else
                g = corrsRunning(k).plane(iPlane).isGad;
            end
            isGad{k} = [isGad{k}; g];
        end
    end 
end

rhos_total = [cat(1,rhosPupil{:}); cat(1,rhosRunning{:})]; % [cat(1,rhosPupil{:}); cat(1,rhosRunning{:})];
% bins = floor((min(rhos_total)+.05)*10)/10 : 0.1 : ceil((max(rhos_total)-.05)*10)/10;

%% (3) Histograms of correlation coefficients
bins = ceil(xLimits(1)*10)/10 : 0.1 : floor(xLimits(2)*10)/10;

condLabels = {'resp.','non-resp.','all'};
cols = lines(2);
for s = 1:2
    maxi = zeros(1, length(condLabels));
    figs = zeros(length(stimuli), length(condLabels));
    for st = 1:length(stimuli)
        if s == 1 % pupil
            sets = ~cellfun(@isempty,rhosPupil(:,st)) & ~cellfun(@isempty,resp);
            rhos = cat(1,rhosPupil{sets,st});
            nulls = cat(1,nullRhosPupil{sets,st});
            responsive = cat(1,resp{sets});
            if isfield(corrsPupil, 'sigma')
                sig = corrsPupil.sigma;
            else
                sig = sigma;
            end
        else % running
%             sets = ~cellfun(@isempty,rhosRunning(:,st)) & ~cellfun(@isempty,resp);
            sets = true(size(rhosRunning,1),1);
            rhos = cat(1,rhosRunning{sets,st});
            nulls = cat(1,nullRhosRunning{sets,st});
%             responsive = cat(1,resp{sets});
            if isfield(corrsRunning, 'sigma')
                sig = corrsRunning.sigma;
            else
                sig = sigma;
            end
        end
        valid = ~isnan(rhos);
%         conds = [responsive==1, responsive==0, ~isnan(responsive)];
        conds = true(length(rhos),3);
        
        prctls = prctile(nulls, [2.5 97.5], 2);
        isSgnf = rhos < prctls(:,1) | rhos > prctls(:,2);
        pVals = sum(nulls < rhos, 2) ./ size(nulls,2);
        ind = pVals > 0.5;
        pVals(ind) = 1 - pVals(ind);
        pVals = 2 .* pVals;
        pVals(pVals==0) = 1/size(nulls,2);
        for c = 3 %1:length(condLabels)
            % Fisher's method to combine p-values
            chi_vals = -2.*log(pVals(valid & conds(:,c)));
            group_pval = 1 - chi2cdf(sum(chi_vals),2*length(chi_vals));
            
            n1 = hist(rhos(isSgnf & valid & conds(:,c)),bins)';
            n2 = hist(rhos(~isSgnf & valid & conds(:,c)),bins)';
            maxi(c) = max(maxi(c), ceil(max(n1+n2)/10)*10);
            % gratings
            figs(st,c) = figure;
            b = bar(bins, [n1,n2], 'stacked');
            b(1).FaceColor = 'k';
            b(2).FaceColor = 'w';
            xlabel('Correlation coeff.')
            ylabel(['#' label])
            n = sum(valid & conds(:,c));
            title(sprintf('%s - %s - %s (n = %d, %d%% signif., p = %.2e, Fisher''s test)', ...
                stimuli{st}, signals{s}, condLabels{c}, n, round(sum(isSgnf & ...
                valid & conds(:,c)) / n*100), group_pval))
            xlim([bins(1)-.05 bins(end)+.05])
            ax = gca;
            ax.Box = 'off';
            legend('p < 0.05', 'p >= 0.05')
        end
    end
%     for c = 1:length(condLabels)
%         for st = 1:length(stimuli)
%             figure(figs(st,c))
%             ylim([0 maxi(c)])
%         end
%     end
end

%% (4) Cumulative distribution of correlation coeffs. (pupil, examples)
leg = {'data', 'null distribution'};

% % for boutons, gratings+gray screen: get data for neurons for comparison --
% leg = {'data', 'null distribution', 'neurons'};
% fR = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\nonVisualEffects\modelGratingResp\';
% data = load(fullfile(fR,'pupil', ...
%     sprintf('corrsDuringGratingsAndGrayScreen_sigma%.2f.mat',sigma)));
% corrsPupil_neurons = data.corrs;
% rhosPupil_neurons = cell(length(corrsPupil_neurons), length(stimuli));
% for k = 1:length(corrsPupil_neurons)
%     for iPlane = 1:length(corrsPupil_neurons(k).plane)
%         for st = 1:length(stimuli)
%             if ~isfield(corrsPupil_neurons(k).plane(iPlane), stimuli{st})
%                 continue
%             end
%             rhosPupil_neurons{k,st} = [rhosPupil_neurons{k,st}; ...
%                 corrsPupil_neurons(k).plane(iPlane).(stimuli{st}).rhos];
%         end
%     end 
% end
% % -------------------------------------------------------------------------
% for neurons, gratings+gray screen: get data for boutons for comparison --
if strcmp(label, 'neurons')
    leg = {'data', 'null distribution', 'boutons'};
    fR = fullfile(folderPC, 'Lab\RESULTS\boutons\nonvisualEffects');
    data = load(fullfile(fR,'pupil', ...
        sprintf('corrsDuringGratingsAndGrayScreen_sigma%.2f.mat',sigma)));
    corrsPupil_neurons = data.corrs;
    rhosPupil_neurons = cell(length(corrsPupil_neurons), length(stimuli));
    for k = 1:length(corrsPupil_neurons)
        for iPlane = 1:length(corrsPupil_neurons(k).plane)
            for st = 1:length(stimuli)
                if ~isfield(corrsPupil_neurons(k).plane(iPlane), stimuli{st})
                    continue
                end
                rhosPupil_neurons{k,st} = [rhosPupil_neurons{k,st}; ...
                    corrsPupil_neurons(k).plane(iPlane).(stimuli{st}).rhos];
            end
        end
    end
end
% -------------------------------------------------------------------------

rhosEx = zeros(numNeurons,length(stimuli),2); % [neurons x stimulus x behaviour]

for sig = 1:2 %2 %1: pupil, 2: running
    if sig == 1
        rhos_ = rhosPupil;
        nulls_ = nullRhosPupil;
        corrs = corrsPupil;
    else
        rhos_ = rhosRunning;
        nulls_ = nullRhosRunning;
        corrs = corrsRunning;
    end
    for st = 1:length(stimuli)
        rhos = cat(1,rhos_{:,st});
        ind = isnan(rhos);
        rhos(ind) = [];
        x = sort(rhos, 'ascend');
        y = (1:length(x))' ./ length(x);
        x = [-1; x; 1];
        y = [0; y; 1];
        
        nulls = cat(1,nulls_{:,st});
        nulls(ind,:) = [];
        
        xNull = sort(nulls, 1, 'ascend');
        xNull = sort(xNull, 2, 'ascend');
        limNull = prctile(xNull, [2.5 97.5], 2);
        limNull = [[-1 -1]; limNull; [1 1]];
        xNull = median(xNull, 2);       
%         xNull = sort(nulls(:), 'ascend');

        yNull = (1:length(xNull))' ./ (length(xNull));
        xNull = [-1; xNull; 1];
        yNull = [0; yNull; 1];
        for n = 1:numNeurons
            rhosEx(n,st,sig) = corrs(exSet).plane(exUnits(n,1)) ...
                .(stimuli{st}).rhos(exUnits(n,2));
        end
        figure
        hold on
        h = [0 0];
        % if neurons, gratings+grayScreen ---------------------------------
        if strcmp(label, 'neurons')
            h = [0 0 0];
            rhos_n = cat(1,rhosPupil_neurons{:,st});
            ind = isnan(rhos_n);
            rhos_n(ind) = [];
            x_n = sort(rhos_n, 'ascend');
            y_n = (1:length(x_n))' ./ length(x_n);
            x_n = [-1; x_n; 1];
            y_n = [0; y_n; 1];
            if ~isempty(x_n)
                h(3) = plot(x_n, y_n, 'Color', [1 1 1].*0.7, 'LineWidth', 2);
            end
        end
        % -----------------------------------------------------------------
        fill([limNull(:,1);flip(limNull(:,2))], [yNull; flip(yNull)], ...
            'k', 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.2)
        h(1) = plot(x, y, 'k', 'LineWidth', 2);
        h(2) = plot(xNull, yNull, 'k:', 'LineWidth', 2);
        heights = interp1(x, y, rhosEx(:,st,sig));
        plot(rhosEx(:,st,sig), heights, 'o', 'MarkerEdgeColor', 'none', ...
            'MarkerFaceColor', 'r')
        xlim(xLimits)
        xlabel(['Correlation with ' signals{sig}])
        ylabel(['Proportion of ' label])
        title(sprintf('Corr. with %s during %s (n = %d)', signals{sig}, ...
            stimuli{st}, length(rhos)))
        legend(h, leg, 'Location', 'NorthWest')
    end
end

%% (4b) Histograms and cumulative distributions of excitatory and inhibitory neurons
bins = ceil(xLimits(1)*10)/10 : 0.1 : floor(xLimits(2)*10)/10;

for st = 1:2
    % pupil
    sets = ~cellfun(@isempty,rhosPupil(:,st)) & ~cellfun(@isempty,resp);
    rhos = cat(1,rhosPupil{sets,st});
    nulls = cat(1,nullRhosPupil{sets,st});
    responsive = cat(1,resp{sets});
    gad = cat(1,isGad{sets});
    if isfield(corrsPupil, 'sigma')
        sig = corrsPupil.sigma;
    else
        sig = sigma;
    end
    valid = ~isnan(rhos);
    prctls = prctile(nulls, [2.5 97.5], 2);
    isSgnf = rhos < prctls(:,1) | rhos > prctls(:,2);
    % excitatory
    n1 = hist(rhos(isSgnf & valid & gad==-1),bins)';
    n2 = hist(rhos(~isSgnf & valid & gad==-1),bins)';
    figure;
    b = bar(bins, [n1,n2], 'stacked');
    b(1).FaceColor = 'r';
    b(2).FaceColor = 'w';
    xlabel('Correlation coeff.')
    ylabel(['#' label])
    n = sum(valid & gad==-1);
    title(sprintf('excitatory - %s - pupil (n = %d, %d%% signif., %.2f s filter)', ...
        stimuli{st}, n, round(sum(isSgnf & valid & gad==-1) / n*100), sig))
    xlim([bins(1)-.05 bins(end)+.05])
    ax = gca;
    ax.Box = 'off';
    legend('p < 0.05', 'p >= 0.05')
    % inhibitory
    n1 = hist(rhos(isSgnf & valid & gad==1),bins)';
    n2 = hist(rhos(~isSgnf & valid & gad==1),bins)';
    figure;
    b = bar(bins, [n1,n2], 'stacked');
    b(1).FaceColor = 'b';
    b(2).FaceColor = 'w';
    xlabel('Correlation coeff.')
    ylabel(['#' label])
    n = sum(valid & gad==1);
    title(sprintf('inhibitory - %s - pupil (n = %d, %d%% signif., %.2f s filter)', ...
        stimuli{st}, n, round(sum(isSgnf & valid & gad==1) / n*100), sig))
    xlim([bins(1)-.05 bins(end)+.05])
    ax = gca;
    ax.Box = 'off';
    legend('p < 0.05', 'p >= 0.05')
    % cumulative
    xEx = sort(rhos(valid & gad==-1 & responsive==1), 'ascend');
    yEx = (1:length(xEx))' ./ length(xEx);
    xEx = [-1; xEx; 1];
    yEx = [0; yEx; 1];
    xIn = sort(rhos(valid & gad==1 & responsive==1), 'ascend');
    yIn = (1:length(xIn))' ./ length(xIn);
    xIn = [-1; xIn; 1];
    yIn = [0; yIn; 1];
    figure
    hold on
    h = [0 0];
    h(1) = plot(xEx, yEx, 'r', 'LineWidth', 2);
    h(2) = plot(xIn, yIn, 'b', 'LineWidth', 2);
    xlim(xLimits)
    xlabel('Correlation with pupil')
    ylabel(['Proportion of ' label])
    [~,p] = kstest2(rhos(valid & gad==-1), rhos(valid & gad==1));
    title(sprintf('Corr. with pupil during %s (p = %.2e, KS-test)', stimuli{st}, p))
    legend(h, {'excitatory','inhibitory'}, 'Location', 'NorthWest')
end

%% (5)  Scatterplot: corr. with pupil vs. with running (during gray screen)
if strcmp(label, 'boutons')
    stim = 1;
else
    stim = 2;
end
ind = ~cellfun(@isempty,rhosRunning(:,stim)) & ~cellfun(@isempty,rhosPupil(:,stim));
r = cat(1,rhosRunning{ind,stim});
p = cat(1,rhosPupil{ind,stim});
figure
plot(r,p,'k.', 'MarkerSize', 3)
ind = ~any(isnan([r,p]),2);
[rho,pVal] = corr(r(ind),p(ind));
xlabel('Corr. with running')
ylabel('Corr. with pupil diameter')
title(sprintf('Correlation of %s during %s (n = %d, rho = %.3f, p = %.2e)', ...
    label, stimuli{stim}, sum(ind), rho, pVal))
axis square equal
hold on
plot(rhosEx(:,stim,2), rhosEx(:,stim,1), '.r', 'MarkerSize', 10)
% plot(rhosExLeast(:,2,2), rhosExLeast(:,2,1), '.b', 'MarkerSize', 10)
axis(repmat(xLimits,1,2))
ax = gca;
ax.Box = 'off';

%% (6)  Scatterplot: corr. with pupil during gray screen vs. during gratings (neurons and boutons)
% OR corr. with running during darkness vs. during gratings (boutons)
stim = 3; % 2 (gray) or 3 (dark)
beh = 2; % 1(pupil) or 2 (running) or 3 (pupil for gratings, running for spont)
switch beh
    case 1
        p_grat = cat(1,rhosPupil{:,1});
        p_spont = cat(1,rhosPupil{:,stim});
    case 2
        p_grat = cat(1,rhosRunning{:,1});
        p_spont = cat(1,rhosRunning{:,stim});
        behEx = 2;
    case 3
        ind = ~cellfun(@isempty,rhosPupil(:,1)) & ...
            ~cellfun(@isempty,rhosRunning(:,stim));
        p_grat = cat(1,rhosPupil{ind,1});
        p_spont = cat(1,rhosRunning{ind,stim});
        behEx = 2;
end
figure
plot(p_spont,p_grat,'k.', 'MarkerSize', 3)
ind = ~any(isnan([p_spont,p_grat]),2);
[rho,p] = corr(p_spont(ind),p_grat(ind));
xlabel(sprintf('Corr. with %s during %s', signals{behEx}, stimuli{stim}))
ylabel(sprintf('Corr. with %s during gratings', signals{behEx}))
title(sprintf('Correlation (n = %d, rho = %.3f, p = %.2e)', sum(ind), rho, p))
axis square
hold on
plot(rhosEx(:,stim,behEx), rhosEx(:,1,1), '.r', 'MarkerSize', 10)
axis(repmat(xLimits,1,2))
ax = gca;
ax.Box = 'off';
axis square

%% (7) Compare running with pupil
% cross-correlations
maxLag = 100; % in sec
binSize = 0.2;
lags = -maxLag : binSize : maxLag;
cc = [];
pupilAll = [];
runningAll = [];
for k = 1:size(running,1)
    for st = 1:size(running,2)
        if isempty(running{k,st}) || isempty(pupil{k,st})
            continue
        end
        t = nonvisTime{k,st}(1) : binSize : nonvisTime{k,st}(end);
        p = interp1(nonvisTime{k,st}, pupil{k,st}, t, 'pchip');
        p = (p - min(p)) ./ (max(p) - min(p));
        r = interp1(nonvisTime{k,st}, running{k,st}, t, 'pchip');
        r = (r - min(r)) ./ (max(r) - min(r));
        cc(:,end+1) = xcorr(zscore(p),zscore(r),(length(lags)-1)/2, 'unbiased');
        pupilAll = [pupilAll; p'];
        runningAll = [runningAll; r'];
    end
end

% cross-correlation
figure
plot(lags, cc, 'Color', [1 1 1].* 0.5)
hold on
plot(lags, mean(cc,2), 'k', 'LineWidth', 2)
xlabel('Time lag (of pupil, in sec)')
ylabel('Cross-correlation')
title('Pupil diameter vs. running speed')

% scatter plot
figure
plot(pupilAll, runningAll, 'k.')
xlabel('Pupil diameter')
ylabel('Running speed')
rho = corr(pupilAll,runningAll);
title(sprintf('Rho = %.3f', rho))

% make conditional probability plot
pBins = prctile(pupilAll,0:10:100);
pBins(end) = pBins(end) + 1;
rBins = unique(prctile(runningAll,0:10:100));
prob = NaN(length(rBins)-1, length(pBins)-1);
for b = 1:length(pBins)-1
    ind = pupilAll>=pBins(b) & pupilAll<pBins(b+1);
    n = histcounts(runningAll(ind), rBins);
    prob(:,b) = n ./ sum(ind);
end
figure
imagesc([5 95], [5 95], prob)
h = colorbar;
h.Label.String = 'P(running|pupil)';
colormap gray
ax = gca;
ax.YDir = 'normal';
axis square
xlabel('Percentile of pupil diameter')
ylabel('Percentile of running speed')

%% (8) Plot traces of all neurons and first PC (one dataset)
exps = [1 2]; % 1: gratings, 2: gray screen, 3: darkness
corrs = corrsPupil;
% corrs = corrsRunning;

valid = cell(1,length(exps));
cellIDs = [];
calciumTraces = cell(1,length(exps));
singleNeurons = cell(1,length(exps));
nonVisual = cell(2,length(exps));
pcs = cell(1,length(exps));
time = cell(1,length(exps));
rhos_b = cell(1,length(exps));
for iPlane = 1:length(db(exSet).planes)
    for exp = 1:length(exps)
        folder = fullfile(folderROIData, db(exSet).subject, ...
            db(exSet).date, num2str(db(exSet).(fields{exps(exp)})));
        % load meta
        data=load(fullfile(folder, sprintf('%s_%d_%s_2P_plane%03d_ROI.mat', ...
            db(exSet).date, db(exSet).(fields{exps(exp)}), db(exSet).subject, ...
            db(exSet).planes(iPlane))));
        meta = data.meta;
        meta.folderTL = strrep(meta.folderTL, 'ioo.ucl.ac.uk', ...
            'cortexlab.net');
        frameTimes = ppbox.getFrameTimes(meta);
        if iPlane == 1
            t_nv = cell(1,2);
            stdSamples = round(sigma / median(diff(frameTimes)));
            convWindow = normpdf(-4*stdSamples:4*stdSamples, 0, stdSamples);
            time{exp} = frameTimes;
            frameDur = median(diff(frameTimes));
            % load running
            ballData = nonVis.getRunningSpeed(meta);
            if ~isempty(ballData)
                nonVisual{1,exp} = ballData.total / median(diff(ballData.t)) / 53; % cm/s
                t_nv{1} = ballData.t;
            end
            % load pupil
            [pupilData, t] = nonVis.loadPupilData(meta);
            if ~isempty(pupilData)
                t(length(pupilData.x)+1:end) = [];
                t_nv{2} = t;
                nonVisual{2,exp} = nonVis.getPupilDiam(pupilData);
            end
            for b = 1:2
                if isempty(nonVisual{b,exp})
                    continue
                end
                ind = isnan(nonVisual{b,exp});
                indInterp = hist(t_nv{b}(ind), time{exp}) > 0;
                nonVisual{b,exp} = interp1(t_nv{b}(~ind), ...
                    nonVisual{b,exp}(~ind), time{exp}, 'pchip')';
                nonVisual{b,exp} = conv([ones(ceil((length(convWindow)-1)/2),1) .* ...
                    mean(nonVisual{b,exp}(1:round(length(convWindow)/2))); ...
                    nonVisual{b,exp}; ones(floor((length(convWindow)-1)/2),1) .* ...
                    mean(nonVisual{b,exp}(end-round(length(convWindow)/2):end))], ...
                    convWindow, 'valid');
                nonVisual{b,exp}(indInterp) = NaN;
            end
            if exps(exp) == 1
                stimTimes = ssLocal.getStimulusResponseInfo(meta);
            end
        end
        % only consider ROIs that are unique and not switch-on
        val = ~all(isnan(meta.F_final),1)';
        if isfield(meta.ROI, 'isDuplicate')
            val = val & meta.ROI.isDuplicate == 0;
        end
        if isfield(meta.ROI, 'isSwitchOn')
            val = val & meta.ROI.isSwitchOn == 0;
        end
        valid{exp} = [valid{exp}; val];
        traces = NaN(length(time{exp}),size(meta.F_final,2));
        for n = 1:size(meta.F_final,2)
            if ~val(n)
                continue
            end
            ind = isnan(meta.F_final(:,n));
            indInterp = hist(frameTimes(ind), time{exp}) > 0;
            badEpisodes = [0 diff(indInterp)];
            starts = find(badEpisodes==1);
            ends = find(badEpisodes==-1);
            if ~isempty([starts ends])
                if isempty(ends) || starts(1)>ends(1), starts = [1 starts]; end
                if isempty(starts) || ends(end)<starts(end), ends = [ends length(time{exp})]; end
            end
            bad = find(ends - starts > 5 / frameDur);
            indInterp = false(1, length(time{exp}));
            for b = bad, indInterp(starts(b):ends(b)) = true; end
            traces(:,n) = interp1(frameTimes(~ind), meta.F_final(~ind, n), ...
                time{exp}, 'pchip');
            traces(indInterp,n) = NaN;
        end
        calciumTraces{exp} = [calciumTraces{exp}, traces];
        if exp == 1
            cellIDs = [cellIDs; [ones(size(meta.F_final,2),1).*iPlane, ...
                (1:size(meta.F_final,2))']];
        end
        
        % for retinal boutons
        rhos_b{exp} = [rhos_b{exp}; corrs(exSet).plane(iPlane) ...
            .(stimuli{exps(exp)}).rhos];
    end
end

filteredTraces = cell(size(calciumTraces));
centerIndex = NaN(1,length(exps));
for exp = 1:length(exps)
    filteredTraces{exp} = NaN(size(calciumTraces{exp}));
    for n = 1:size(filteredTraces{exp},2)
        tr = calciumTraces{exp}(:,n);
        filteredTraces{exp}(:,n) = ...
            conv([ones(ceil((length(convWindow)-1)/2),1) .* ...
            mean(tr(1:round(length(convWindow)/2))); ...
            tr; ...
            ones(floor((length(convWindow)-1)/2),1) .* ...
            mean(tr(end-round(length(convWindow)/2)))], ...
            convWindow, 'valid');
    end
    % for sSC neurons
    centerIndex(exp) = floor(size(calciumTraces{exp},1)/2);
    ind = true(size(calciumTraces{exp},1),1);
    ind(centerIndex(exp)+1:end) = false;
    ind(isnan(nonVisual{2,exp})) = false;
    rhos_b{exp} = corr(nonVisual{2,exp}(ind), ...
        filteredTraces{exp}(ind,:))';
    % for boutons
%     centerIndex(exp) = 0;
end

valid = all(cat(2,valid{:}),2);
cellIDs = cellIDs(valid,:);
% r_p = rhos_b{1}; % sort on gratin experiment
r_p = mean(cat(2,rhos_b{:}),2); % sort by mean (for boutons: gratings and darkness)
r_p = r_p(valid);
[r_p_sorted,order] = sort(r_p,'descend');
cellIDs2 = cellIDs(order,:);
% cols = 'rrrrrrrrrrbbbbbbbbbb';
cols = 'rb';
for exp = 1:length(exps)
    traces = calciumTraces{exp}(:,valid);
    filtered = filteredTraces{exp}(:,valid);
%     for n = 1:size(filtered,2)
%         tr = traces(:,n);
%         filtered(:,n) = ...
%             conv([ones(ceil((length(convWindow)-1)/2),1) .* ...
%             mean(tr(1:round(length(convWindow)/2))); ...
%             tr; ...
%             ones(floor((length(convWindow)-1)/2),1) .* ...
%             mean(tr(end-round(length(convWindow)/2)))], ...
%             convWindow, 'valid');
%     end
    
    zscored = (filtered - nanmean(filtered,1)) ./ nanstd(filtered,0,1) ./ ...
        size(filtered,1).^0.5;
    nInd = any(isnan(zscored),1); % neurons with NaN values
    if sum(nInd) < 0.1 * size(zscored,2)
        ind = true(size(zscored,1),1);
    else
        nInd = false(1, size(zscored,2));
        ind = ~any(isnan(zscored),2);
    end
    [U,S,V] = svdecon(zscored(ind,~nInd)');
    zscored(:,nInd) = NaN;
    zscored = zscored(:,order);
    zscored(:,all(isnan(zscored),1)) = [];
    pc = NaN(size(zscored,1),10);
    pc(ind,:) = zscore(V(:,1:10));
    
    r = nonVisual{1,exp};
%     ind = ~isnan(r);
%     r(ind) = zscore(r(ind));
    ind = ~isnan(r) & ~any(isnan(pc),2);
    r_corr = corr(pc(ind,:), r(ind));
    p = nonVisual{2,exp};
    if ~isempty(p)
%         ind = ~isnan(p);
%         p(ind) = zscore(p(ind));
        ind = ~isnan(p) & ~any(isnan(pc),2);
        p_corr = corr(pc(ind,:), p(ind));
        [~,best] = max(abs(p_corr+r_corr));
        r_corr = r_corr(best);
        p_corr = p_corr(best);
        pc = pc(:,best);
        if p_corr < 0
            pc = -pc;
            r_corr = -r_corr;
            p_corr = -p_corr;
        end
        ind = ~isnan(p);
        p = interp1(time{exp}(ind), p(ind), time{exp}, 'pchip');
    else
        p = NaN(1,length(time{exp}));
        p_corr = NaN;
        [~,best] = max(abs(r_corr));
        r_corr = r_corr(best);
        pc = pc(:,best);
        if r_corr < 0
            pc = -pc;
            r_corr = -r_corr;
        end
    end
    
    cm = flip(gray);
    ax = [0 0];
    figure('Position',[1921 1 1920 1123])
    subplot(5,1,1)
    hold on
    plot(time{exp}, p, 'Color', [0 0.7 0.5])
    plot(time{exp}, r-3, 'Color', [0 0 .7])
    plot(time{exp}, pc-6, 'Color', [0.7 0.2 0.2])
    for ex = 1:size(exUnits,1)
        ind = find(all(cellIDs==exUnits(ex,:),2));
        plot(time{exp}, smooth(traces(:,ind),5)./2-6-4*ex, cols(ex))
%         plot(time{exp}, zscore(smooth(traces(:,ind),5))./2-6-4*ex, cols(ex))
        fprintf('rho = %.3f\n',corr(smooth(traces(:,ind),5),r))
    end
    if exps(exp) == 1 % gratings
        t = [repmat(stimTimes.onset',2,1);NaN(1,length(stimTimes.onset))];
        s = [zeros(1,length(stimTimes.onset)); ...
            ones(1,length(stimTimes.onset));NaN(1,length(stimTimes.onset))];
%         t = [0; reshape([repmat(stimTimes.onset',2,1); ...
%             repmat(stimTimes.offset',2,1)],[],1); time{exp}(end)];
%         s = [0; repmat([0 1 1 0]',length(stimTimes.onset),1); 0];
        plot(t,s-20,'k')
    else % gray screen
        plot(time{exp}([1 end]),[0 0]-20,'k')
    end
%     ylim([-22 max([p 1])])
    % for sSC neurons
    title(sprintf('%s, %s %s, %dth PC - running: %.2f, - pupil: %.2f (sorted on first half)', ...
        stimuli{exps(exp)}, db(exSet).subject, db(exSet).date, best, r_corr, p_corr), ...
        'Interpreter', 'none')
    % for boutons
%     title(sprintf('%s, %s %s, %dth PC - running: %.2f, - pupil: %.2f (sorted by pupil in gratings)', ...
%         stimuli{exps(exp)}, db(exSet).subject, db(exSet).date, best, r_corr, p_corr), ...
%         'Interpreter', 'none')
    ax(1) = gca;
    leg = legend('Pupil','Running','PC',num2str(exUnits(1,:)), ...
        num2str(exUnits(2,:)),'stimulus');
    leg.Position = [0.92,0.83,0.05,0.09];
    subplot(5,1,2:5)
    imagesc(time{exp}([1 end]),[1 size(zscored,2)], zscored', ...
        prctile(zscored(:),[5 95]))
    colormap(cm)
    ax(2) = gca;
    linkaxes(ax,'x')
    xlim(time{exp}([centerIndex(exp)+1 end]))
    xlabel('Time (s)')
    set(gca,'box','off')
end