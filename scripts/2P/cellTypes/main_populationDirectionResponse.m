%% Define data

% label = 'neurons';
label = 'boutons';

% file location depending on PCs
folderPC = 'C:\STORAGE\OneDrive - University College London';
% folderPC = 'C:\Users\skgtchr\OneDrive - University College London'; % ZEN notebook
% folderPC = 'C:\Storage\OneDrive - University College London'; % Lenovo laptop
% data folders
folderROIData = fullfile(folderPC, 'Lab\DATA\InfoStructs');
if strcmp(label, 'neurons')
    folderResults = fullfile(folderPC, 'Lab\RESULTS\nonvisualEffects\modelGratingResp');
    
    corrections = [];
else
    folderResults = fullfile(folderPC, 'Lab\RESULTS\boutons\nonvisualEffects');
    
    data = load(fullfile(folderROIData, 'corrections_boutons.mat'));
    corrections = data.corrections;
    doCorrect = data.doCorrect;
end

data = load(fullfile(folderResults, 'pupil', ...
    'tuning_prefDirSigmaDIFixed_lineFit.mat'));
tuning = data.tuning;
degrees = data.x;
data = load(fullfile(folderResults, 'pupil', ...
    'nullTuning_prefDirSigmaDIFixed.mat'));
null = data.null;

%% Collect variables of tuning curves

prefDirs = [];
responses = cell(1, length(tuning));
isSuppr = [];
dataset = [];
for iExp = 1:length(tuning)
    fprintf('Dataset %d: %s %s exp.: %d\n', iExp, tuning(iExp).subject, ...
        tuning(iExp).date, tuning(iExp).exp);
    respD = [];
    for iPlane = 1:length(tuning(iExp).plane)
        data = tuning(iExp).plane(iPlane);
        if isempty(data.cellIDs)
            continue
        end
        neurons = find(~cellfun(@isempty, {data.cond(1).cell.parameters}));
        
        for j = 1:length(neurons)
            iCell = neurons(j);
            prefDir = NaN;
            resp = NaN(1,12,size(data.cond(1).cell(iCell).responses,2),2);
            for c = 1:2
                resp(1,:,:,c) = data.cond(c).cell(iCell).responses;
                
                pars = data.cond(c).cell(iCell).parameters;
                curve = data.cond(c).cell(iCell).curve;
                if length(pars) > 1 % tuned
                    prefDir = pars(1);
                end
            end
            
            if ~isempty(corrections)
                a = corrections(iExp).plane(iPlane).a{tuning(iExp).exp} ...
                    (data.cellIDs(iCell));
                b = corrections(iExp).plane(iPlane).b{tuning(iExp).exp} ...
                    (data.cellIDs(iCell));
                resp = doCorrect(a,b,resp);
            end
            
            prefDirs(end+1,1) = prefDir;
            respD(end+1,:,:,:) = resp;
            isSuppr(end+1,:) = data.isSuppressed(iCell);
        end
        dataset = [dataset; ones(length(neurons),1).*iExp; neurons'];
        responses{iExp} = respD;
    end
end

%% Determine mean responses and DSIs and OSIs
cond = 1; % small or large pupil
shuffles = 1000;
directions = 0:30:330;
dirVectors = exp(directions./180.*pi .* 1i);
oriVectors = exp(directions./180.*2.*pi .* 1i);

meanResp = [];
shuffleResp = [];
for k = 1:length(responses)
    r = responses{k}(:,:,:,cond);
    meanResp = [meanResp; nanmean(r, 3)];
    r1 = reshape(r, size(r,1), []);
    shR = NaN(size(r,1),length(dirVectors),shuffles);
    for sh = 1:shuffles
        p = randperm(size(r1,2));
        r2 = reshape(r1(:,p), size(r));
        shR(:,:,sh) = nanmean(r2, 3);
    end
    shuffleResp = [shuffleResp; shR];
end

meanResp(isSuppr==1,:) = -meanResp(isSuppr==1,:);
shuffleResp(isSuppr==1,:,:) = -shuffleResp(isSuppr==1,:,:);
meanResp(meanResp<0) = 0;
shuffleResp(shuffleResp<0) = 0;
meanResp = meanResp ./ max(meanResp,[],2);
shuffleResp = shuffleResp ./ max(shuffleResp,[],2);

% Determine DSIs
vects = mean(dirVectors .* meanResp, 2);
shuffleVects = squeeze(mean(dirVectors .* shuffleResp, 2));
DSIs = abs(vects);
nullDSIs = abs(shuffleVects);
p_DSI = sum(nullDSIs > DSIs,2) ./ shuffles;
p_DSI(isnan(DSIs)) = NaN;

% Determine OSIs
vects = mean(oriVectors .* meanResp, 2);
shuffleVects = squeeze(mean(oriVectors .* shuffleResp, 2));
OSIs = abs(vects);
nullOSIs = abs(shuffleVects);
p_OSI = sum(nullOSIs > OSIs,2) ./ shuffles;
p_OSI(isnan(OSIs)) = NaN;

%% Plot population direction tuning curve
groups = [true(size(meanResp,1),1), ... % all
    p_DSI<.05|p_OSI<.05, ...            % exclude non-selective
    p_DSI<.05&p_OSI>=.05, ...           % only dir- but NOT ori-selective
    p_DSI>=.05&p_OSI<.05, ...           % only ori- but NOT dir-selective
    p_DSI<.05&p_OSI<.05];               % boht dir- AND ori-selective
groupNames = {'all','DS or OS','DS only','OS only','DS and OS'};
colors = [0 0 0; 0 0 0; 1 0 0; 0 0 1; .5 0 .5];

for g = 1:size(groups,2)
    figure
    hold on
    m = nanmean(meanResp(groups(:,g),[1:end,1]),1);
    sem = nanstd(meanResp(groups(:,g),[1:end,1]),0,1) ./ sqrt(sum(groups(:,g)));
    fill([directions 360 360 flip(directions)], [m+sem flip(m-sem)], ...
        'k', 'EdgeColor', 'none', 'FaceColor', colors(g,:), 'FaceAlpha', .2)
    plot([directions 360], m, 'Color', colors(g,:), 'LineWidth', 2)
    xlim([0 360])
    ylim([0.4 0.8])
    set(gca, 'box', 'off', 'XTick', 0:90:360)
    xlabel('Direction')
    ylabel('Mean population response')
    title(sprintf('%s, %s (n = %d)', label, groupNames{g}, sum(groups(:,g))))
end