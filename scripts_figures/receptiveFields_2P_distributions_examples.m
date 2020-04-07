%% Define data

% label = 'neurons';
label = 'boutons';

if strcmp(label, 'neurons')
    exampleSets(1).subject = 'M150410_SS044'; % corr gratings: -0.2212, p = 0.024
    exampleSets(1).date = '2015-05-29';       % no gray screen
    exampleSets(1).plane = 2;                 % DI pref resp: -0.0776, p = 0.02
    exampleSets(1).neuron = 72;
    exampleSets(1).RFtype = 'ON';
    exampleSets(1).corrDirection = -1; 
    exampleSets(2).subject = 'M150410_SS044'; % corr gratings: 0.2021, p = 0.036
    exampleSets(2).date = '2015-04-28';       % corr gray screen: -0.0144, p = 0.436
    exampleSets(2).plane = 2;                 % DI pref resp: 0.1324, p = 0.01
    exampleSets(2).neuron = 204;
    exampleSets(2).RFtype = 'ON';
    exampleSets(2).corrDirection = 1;
    exampleSets(3).subject = 'M150305_SS041'; % corr gratings: -0.2171, p = 0.034
    exampleSets(3).date = '2015-04-23';       % corr gray screen: -0.2397, p = 0.002
    exampleSets(3).plane = 4;                 % DI pref resp: -0.1449, p = 0.01
    exampleSets(3).neuron = 58;
    exampleSets(3).RFtype = 'OFF';
    exampleSets(3).corrDirection = -1;
    exampleSets(4).subject = 'M150611_SS048'; % corr gratings: 0.3802, p = 0.002
    exampleSets(4).date = '2015-12-02';       % corr gray screen: 0.5246, p = 0.01
    exampleSets(4).plane = 3;                 % DI pref resp: 0.1956, p = 0.02
    exampleSets(4).neuron = 62;               % alternatives: set 10, pl 3, neuron 15;
    exampleSets(4).RFtype = 'OFF';            % set 3, pl 3, neuron 159
    exampleSets(4).corrDirection = 1;         % set 3, pl 4, neuron 60
    exampleSets(5).subject = 'M150305_SS041'; % corr gratings: -0.3017, p = 0.0000
    exampleSets(5).date = '2015-04-23';       % corr gray screen: -0.5497, p = 0.000
    exampleSets(5).plane = 4;                 % DI pref resp: -0.2956, p = 0.0000
    exampleSets(5).neuron = 43;
    exampleSets(5).RFtype = 'ON+OFF';
    exampleSets(5).corrDirection = -1;
    exampleSets(6).subject = 'M150410_SS044';
    exampleSets(6).date = '2015-04-28';
    exampleSets(6).plane = 3;
    exampleSets(6).neuron = 188;
    exampleSets(6).RFtype = 'ON+OFF';
    exampleSets(6).corrDirection = 1;
else % boutons in sSC
    exampleSets(1).subject = 'SS078';
    exampleSets(1).date = '2017-10-05';
    exampleSets(1).plane = 1;
    exampleSets(1).neuron = 156;
    exampleSets(1).RFtype = 'ON';
    exampleSets(1).corrDirection = -1;
    exampleSets(2).subject = 'SS078';
    exampleSets(2).date = '2017-10-05';
    exampleSets(2).plane = 1;
    exampleSets(2).neuron = 225;
    exampleSets(2).RFtype = 'ON';
    exampleSets(2).corrDirection = 1;
%     exampleSets(3).subject = 'SS076';
%     exampleSets(3).date = '2017-10-02';
%     exampleSets(3).plane = 1;
%     exampleSets(3).neuron = 100;
%     exampleSets(3).RFtype = 'OFF';
%     exampleSets(3).corrDirection = -1;
    exampleSets(3).subject = 'SS078';
    exampleSets(3).date = '2017-10-05';
    exampleSets(3).plane = 1;
    exampleSets(3).neuron = 22;
    exampleSets(3).RFtype = 'OFF';
    exampleSets(3).corrDirection = -1;
    exampleSets(4).subject = 'SS078';
    exampleSets(4).date = '2017-09-28';
    exampleSets(4).plane = 1;
    exampleSets(4).neuron = 243;
    exampleSets(4).RFtype = 'OFF';
    exampleSets(4).corrDirection = 1;
    exampleSets(5).subject = 'M160923_SS069';
    exampleSets(5).date = '2016-10-21';
    exampleSets(5).plane = 1;
    exampleSets(5).neuron = 223;
    exampleSets(5).RFtype = 'ON+OFF';
    exampleSets(5).corrDirection = -1;
    exampleSets(6).subject = 'SS078';
    exampleSets(6).date = '2017-09-28';
    exampleSets(6).plane = 1;
    exampleSets(6).neuron = 263;
    exampleSets(6).RFtype = 'ON+OFF';
    exampleSets(6).corrDirection = 1;
end

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
    
    corrections = [];
else
    folderResults = fullfile(folderPC, 'Lab\RESULTS\boutons\nonvisualEffects');
    folderRF = fullfile(folderPC, 'Lab\RESULTS\receptiveFields\boutons');
    
    data = load(fullfile(folderROIData, 'corrections_boutons.mat'));
    corrections = data.corrections;
    doCorrect = data.doCorrect;
end

%% Parameters
RFtypes = {'ON','OFF','ON+OFF'};

% for evaluation of receptive fields (significance/goodness)
minExplainedVarianceStim = 0.01;
minPVal = 0.05;
maxLambda = 1;

% to group into ON, OFF and ON+OFF
onThr = 0.5; % ON field is at least three times stronger than OFF field
offThr = -0.5; % OFF field is at least three times stronger than ON field

% colormaps
red = [1 0 .5];
% green = [0 1 .5];
blue = [0 .5 1];
% orange = [1 .5 0];
black = [1 1 1].*0.5;
grad = linspace(0,1,100)';
reds = red.*flip(grad) + [1 1 1].*grad;
% greens = green.*flip(grad) + [1 1 1].*grad;
blacks = black.*flip(grad) + [1 1 1].*grad;
% cm_ON = [greens; flip(reds(1:end-1,:),1)];
cm_ON = [blacks; flip(reds(1:end-1,:),1)];
blues = blue.*flip(grad) + [1 1 1].*grad;
% oranges = orange.*flip(grad) + [1 1 1].*grad;
% cm_OFF = [oranges; flip(blues(1:end-1,:),1)];
cm_OFF = [blacks; flip(blues(1:end-1,:),1)];
colormaps = {cm_ON, cm_OFF};

%% Load data
data = load(fullfile(folderRF, 'receptiveFields.mat'));
RFs = data.RFs;

%% Plot RFs of example units
titles = {'ON field','OFF field'};
for ex = 1:length(exampleSets)
    s = strcmp({RFs.subject}, exampleSets(ex).subject) & ...
        strcmp({RFs.date}, exampleSets(ex).date);
    rf = RFs(s).plane(exampleSets(ex).plane).receptiveFields(:,:,:,:, ...
        exampleSets(ex).neuron);
    rf(:,:,:,2) = -rf(:,:,:,2);
    pos = RFs(s).stimPosition;
    squW = diff(pos(1:2)) / size(rf,2);
    squH = diff(pos(3:4)) / size(rf,1);
    [mx,mxTime] = max(max(reshape(permute(abs(rf),[1 2 4 3]),[],size(rf,3)),[],1));
    % fit Gaussian
    fitPars = whiteNoise.fit2dGaussRF(mean(rf(:,:,mxTime,:),4), false);
    [x0, y0] = meshgrid(linspace(0.5, size(rf,2)+0.5, 100), ...
        linspace(0.5, size(rf,1)+0.5, 100));
    fitRF = whiteNoise.D2GaussFunctionRot(fitPars, cat(3, x0, y0));
    [x1,y1] = meshgrid(linspace(pos(1), pos(2), 100), ...
        linspace(pos(3), pos(4), 100));
    figure
    for f = 1:2
        subplot(2,1,f)
        imagesc([pos(1)+squW/2 pos(2)-squW/2], [pos(3)+squH/2 pos(4)-squH/2], ...
            rf(:,:,mxTime,f),[-mx mx])
        hold on
        contour(x1, y1, fitRF, [1 1].*fitPars(1)/2, 'k', 'LineWidth', 1)
        axis image
        set(gca, 'box', 'off', 'XTick', pos(1:2), 'YTick', [pos(3) 0 pos(4)], ...
            'YTickLabel', [-pos(3) 0 -pos(4)])
        colormap(gca, colormaps{f})
        ylabel(titles{f})
        colorbar
        if pos(2) > 0 % RF map crosses vertical midline
            xlim([pos(1) 0])
            set(gca, 'XTick', [pos(1) 0])
        end
        if strcmp(label, 'neurons') && pos(2) > -75 % RF map crosses vertical midline
            xlim([pos(1) -75])
            set(gca, 'XTick', [pos(1) -75])
        end
        if f == 1
            title(sprintf('%s, %s, plane %d, unit %d (corr with pupil: %d)', ...
                exampleSets(ex).subject, exampleSets(ex).date, ...
                exampleSets(ex).plane, exampleSets(ex).neuron, ...
                exampleSets(ex).corrDirection), 'interpreter', 'none')
        end
    end
end

%% Collect relevant variables from visual noise data
% for running + stimulus, select best stimulus model (linear, absolute,
% white, or black)

datasets = [];
planes = [];
neurons = [];
evTotal = [];
evStim = [];
evRun = [];
lambdasStim = [];
pValues = [];
OnOffValues = [];

for k = 1:length(RFs)
    for iPlane = 1:length(RFs(k).plane)
        n = length(RFs(k).plane(iPlane).explainedVariances);
        datasets = [datasets; ones(n,1).*k];
        planes = [planes; ones(n,1).*iPlane];
        neurons = [neurons; (1:n)'];
        evTotal = [evTotal; RFs(k).plane(iPlane).explainedVariances];
        evStim = [evStim; RFs(k).plane(iPlane).explainedVariances_stimOnly];
        evRun = [evRun; RFs(k).plane(iPlane).explainedVariances_runOnly];
        lambdasStim = [lambdasStim; RFs(k).plane(iPlane).lambdasStim'];
        pValues = [pValues; RFs(k).plane(iPlane).pVal_RFonly];
        
        oov = NaN(n,2);
        for iCell = 1:n
            rf = RFs(k).plane(iPlane).receptiveFields(:,:,:,:,iCell);
            [~,t] = max(max(reshape(permute(abs(rf),[1 2 4 3]), [], ...
                size(rf,3)), [], 1));
            rf = squeeze(rf(:,:,t,:));
            [mx,row] = max(max(abs(rf),[],3),[],1);
            [~,col] = max(mx);
            row = row(col);
            rf = squeeze(rf(row,col,:));
            rf(2) = -rf(2);
            oov(iCell,:) = rf;
        end
        OnOffValues = [OnOffValues; oov];
    end
end

[~,type] = max(abs(OnOffValues),[],2);
ind = sub2ind(size(OnOffValues), (1:size(OnOffValues,1))', type);
signs = sign(OnOffValues(ind));
OnOffRatios = OnOffValues;
OnOffRatios(signs<0,:) = -OnOffRatios(signs<0,:);
OnOffRatios(OnOffRatios<0) = 0;
OnOffRatios = (OnOffRatios(:,1)-OnOffRatios(:,2))./sum(OnOffRatios,2);

validRF = pValues < minPVal & evStim > minExplainedVarianceStim & ...
    lambdasStim < maxLambda;

%% Plot distribution of ON-OFF-ratios
exSets = zeros(1,length(exampleSets));
for ex = 1:length(exampleSets)
    exSets(ex) = find(strcmp({RFs.subject}, exampleSets(ex).subject) & ...
        strcmp({RFs.date}, exampleSets(ex).date));
end

binSize = 0.05;
edges = -1:binSize:1;
bins = edges(1:end-1)+binSize/2;
figure
histogram(OnOffRatios(validRF), edges, 'FaceColor', 'k')
n = histcounts(OnOffRatios(validRF), edges);
mx = max(n) * 1.05;
hold on
plot([1 1].*offThr, [0 mx], 'k')
plot([1 1].*onThr, [0 mx], 'k')
for ex = 1:length(exampleSets)
    switch exampleSets(ex).RFtype
        case 'ON'
            color = red;
        case 'OFF'
            color = blue;
        otherwise
            color = (red+blue)./2;
    end
    if exampleSets(ex).corrDirection > 0
        color = color .* 0.7 + [1 1 1].*0.3;
    else
        color = color .* 0.7 + [0 0 0].*0.3;
    end
    
    plot(OnOffRatios(datasets==exSets(ex) & planes==exampleSets(ex).plane & ...
        neurons==exampleSets(ex).neuron), mx, 'v', ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'none')        
end
set(gca, 'box', 'off', 'XTick', [-1 offThr 0 onThr 1])
xlim([-1 1])
ylim([0 mx])
xlabel('ON/OFF-ratio ((ON-OFF)/(ON+OFF))')
ylabel(sprintf('#%s',label))
title(sprintf('ON/OFF ratios (n=%d)', sum(validRF)))

%% Collect correlation values from analyses of continuous traces
data = load(fullfile(folderCorrResults, nonvis, ...
    'corrsDuringGratingsAndGrayScreen_sigma1.00.mat'));
corrs = data.corrs;

rhos = [];
nullRhos = [];
resp = [];
isGad = [];
n = 0;
for k = 1:length(RFs)
    if isempty(RFs(k).RFTimes)
        continue
    end
    kcorr = find(strcmp(RFs(k).subject, {corrs.subject}) & ...
        strcmp(RFs(k).date, {corrs.date}));
    if isempty(kcorr)
        numCells = sum(cellfun(@length, {RFs(k).plane.explainedVariances_runOnly}));
        rhos = [rhos; NaN(numCells, length(stimuli))];
        nullRhos = [nullRhos; NaN(numCells, length(stimuli), ...
            size(corrs(1).plane(1).gratings.nullRhos,2))];
        resp = [resp; NaN(numCells, 1)];
        if isfield(corrs(1).plane, 'isGad')
            isGad = [isGad; NaN(numCells, 1)];
        end
        n = n + numCells;
        continue
    end
    for iPlane = 1:length(corrs(kcorr).plane)
        fields = fieldnames(corrs(kcorr).plane(iPlane));
        st = find(~structfun(@isempty, corrs(kcorr).plane(iPlane)),1);
        numCells = length(corrs(kcorr).plane(iPlane).(fields{st}).rhos);
        rhos = [rhos; NaN(numCells, length(stimuli))];
        nullRhos = [nullRhos; NaN(numCells, length(stimuli), ...
            size(corrs(1).plane(1).gratings.nullRhos,2))];
        for st = 1:length(stimuli)
            if ~isfield(corrs(kcorr).plane(iPlane), stimuli{st})
                continue
            end
            rhos((1:numCells) + n, st) = ...
                corrs(kcorr).plane(iPlane).(stimuli{st}).rhos;
            nullRhos((1:numCells) + n, st, :) = ...
                corrs(kcorr).plane(iPlane).(stimuli{st}).nullRhos;
        end
        if isfield(corrs(kcorr).plane(iPlane), 'responsive')
            resp = [resp; corrs(kcorr).plane(iPlane).responsive];
        else resp = [resp; NaN(numCells,1)];
        end
        if isfield(corrs(kcorr).plane, 'isGad')
            g = corrs(kcorr).plane(iPlane).isGad;
            isGad = [isGad; g];
        end
        
        n = n + numCells;
    end 
end

%% Plot running/pupil correlation for different RF types
indsEnh = {maxSignStim > 0, maxSignStim < 0};

% Cum. histogram of correlation coefficients (data and null distribution)
% for each RF type (On enh., On supp., Off enh., Off supp., ...)
RFTypes = {'Linear', 'Absolute', 'On', 'Off'};
respTypes = {'enhanced', 'suppressed'};
exPrctInDistrib = [0.1 0.5 0.9];
% xLimits = [floor(min(rhos(:))*10)/10, ceil(max(rhos(:))*10)/10];
xLimits = [-0.65 0.85];
colors = lines(3);
for stim = 1:length(stimuli)
    if all(isnan(rhos(:,stim)))
        continue
    end
    
    figure
    h = [0 0 0];
    ind = model==0 & ~isnan(rhos(:,stim));
    h0 = general.plotCumHist(gca, rhos(ind,stim), squeeze(nullRhos(ind,stim,:)), ...
        xLimits, [0 0 0; 0.5 0.5 0.5]);
    h(1:2) = h0;
    h0 = general.plotCumHist(gca, rhos(~isnan(rhos(:,stim)),stim), [], ...
        xLimits, [0 0 0], {'--'});
    h(3) = h0;
    xlim(xLimits)
    xlabel(['Correlation with ' nonvis])
    ylabel(['Proportion of ' label])
    title(sprintf('No RF, correlation during %s (n = %d)', stimuli{stim}, sum(ind)))
    legend(h, {'no RF', 'shifted', 'all'})
%     fprintf('%s at 10%%, 50%%, and 90%% of distribution, for corr during %s\n', ...
%         label, stimuli{stim})
    for m = 1:4
        for r = 1:2
            figure
            h = [0 0 0];
            ind = find(model==m & indsEnh{r} & ~isnan(rhos(:,stim)));
            h0 = general.plotCumHist(gca, rhos(ind,stim), squeeze(nullRhos(ind,stim,:)), ...
                xLimits, [0 0 0; 0.5 0.5 0.5]);
            h(1:2) = h0;
            h0 = general.plotCumHist(gca, rhos(~isnan(rhos(:,stim)),stim), [], ...
                xLimits, [0 0 0], {'--'});
            h(3) = h0;
            xlim(xLimits)
            xlabel(['Correlation with ' nonvis])
            ylabel(['Proportion of ' label])
            title(sprintf('%s %s, correlation during %s (n = %d)', RFTypes{m}, ...
                respTypes{r}, stimuli{stim}, length(ind)))
            legend(h, {[RFTypes{m} ' ' respTypes{r}], 'shifted', 'all'})
            
            % select examples
            [~,j] = sort(rhos(ind,stim));
            jind = ceil(exPrctInDistrib .* length(ind));
            if m==3 && r==1 && strcmp(label,'neurons') 
                jind(1) = 3;
            end
            ex = ind(j(jind));
            hold on
            for n = 1:3
                plot(rhos(ex(n),stim), jind(n)./length(ind), '.', ...
                    'Color', colors(n,:), 'MarkerSize', 30)
            end
%             fprintf('  %s RF, %s (Dataset Plane NeuronID correlation)\n', RFTypes{m}, respTypes{r})
%             disp([dataset(ex) plane(ex) neuron(ex) rhos(ex,stim)])
            % plot RFs of examples
            for n = 1:3
                figure
                rf = RFs(dataset(ex(n))).plane(plane(ex(n))) ...
                    .receptiveFields(:,:,:,neuron(ex(n)),m);
                cnf = squeeze(prctile(nullRhos(ex(n),stim,:), [2.5 97.5], 3));
                pos = RFs(dataset(ex(n))).stimPosition;
                [maxVal,maxFrame] = max(max(max(abs(rf),[],1),[],2),[],3);
                imagesc(rf(:,:,maxFrame),[-maxVal maxVal])
                colormap(cm)
                axis image
                set(gca, 'XTick', [.5 size(rf,2)+.5], 'XTickLabel', pos(1:2), ...
                    'YTick', [.5 size(rf,1)+.5], 'YTickLabel', pos(3:4))
                title(sprintf('%s RF (%s) with corr at %d%% (%.3f [%.3f %.3f])\n(set %d, plane %d, neuron %d)', ...
                    RFTypes{m}, respTypes{r}, exPrctInDistrib(n)*100, rhos(ex(n),stim), ...
                    cnf(1), cnf(2), dataset(ex(n)), plane(ex(n)), neuron(ex(n))))
            end
        end
    end
end