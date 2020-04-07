% Analysis of 2p retinotopy, with responses as filters on stimulus onsets


addpath('\\ZSERVER\Code\Matteobox');


%% Load data, name variables, etc


load('AR087-2013-10-30_3_AR087');


calciumTrace    = outputStruct.calciumTrace([1 2 3], :); % show only cells 1 2 3
stimulusMatrix  = outputStruct.stimulusMatrix;
timeFrames      = outputStruct.timeFrames;

% name some variables, etc

Duration = range(timeFrames);
SamplingFrequency = length(timeFrames)/Duration; % sampling frequency

fs = SamplingFrequency;
tt = timeFrames'-timeFrames(1);

nT = length(tt);
nCells = size(calciumTrace,1);

rrr = calciumTrace';

% process the stimuli

nS = size(stimulusMatrix,1); % number of stimuli
% extract onsets
aa = [ max(0,diff(stimulusMatrix')); zeros(1,nS)];
% the last stimulus is a blank. Drop it
aa = aa(:,1:end-1); nS = nS - 1;
% find the times of the stimuli
[ii, ss] = find(aa);

% figure out the stimuli

load('output1-plotlabels.mat');
coords = zeros(nS,2);
for iS = 1:nS, coords(iS,:) = sscanf(condNames{iS},'%*d) x1 = %d\ny1 = %d'); end

%% find the filters

MyKernel = kernel('samplerate',fs);
MyKernel.downSampleRate = 4;
MyKernel.tSamples = 16;         % seconds, total duration of the kernel
MyKernel.tFutureSamples = 12;   % seconds, duration into the future

stim = aa;
[MyKernel,kernels, kstim, rrrDs] = MyKernel.buildKernel( stim, rrr);

pppDs = MyKernel.predictResponses(kstim); % takes a long time!
 
ShiftTimes = MyKernel.shiftValues;

% MyKernel.kernels is nS*nShiftsDs x nCells
fff = zeros(nS,MyKernel.numShifts,nCells);
for iCell = 1:nCells
    for iShift = 1:MyKernel.numShifts
        fff(:,iShift,iCell) = MyKernel.kernels((iShift-1)*nS+(1:nS),iCell);
    end
end

% predicted responses
ppp = imresize( pppDs, [nT, nCells] );



%% Separate time courses from tuning curves

% Baselines = zeros(nCells,nS);

DoGraphics = 0; % set it to 1 if you want to see what is happening

TimeCourses = cell(nCells,1);
Tunings = cell(nCells,1);
for iCell = 1:nCells
    % calculate the baseline and remove it
    % Baselines(iCell,:) = mean( fff(:, ShiftTimes < 0, iCell), 2 );
    % CorrectedFilter = fff(:,:,iCell) - repmat( Baselines(iCell,:)', [1 length(Shifts)]);
    CorrectedFilter = fff(:,:,iCell);
    
    % approximate with separable functions
    [TimeCourse,Tuning,BestScl] = MakeSeparable( CorrectedFilter,DoGraphics, 1);
    
    % make the TimeCourses peak at 1, so the amplitude is in the Tunings
    TimeCourses{iCell} = TimeCourse/max(TimeCourse);
    Tunings    {iCell} = Tuning*BestScl*max(TimeCourse);
    
end

%% try to predict responses with the separable filters

sss = zeros(size(fff)); % separable approximations to the filters
for iCell = 1:nCells
    for iS = 1:nS
       sss(iS,:,iCell) = TimeCourses{iCell}*Tunings{iCell}(iS);
    end
end

%% Graphics: compare predicted to actual responses

figure('Color','w'); clf; ax = gridplot(nCells+1,1);
set(ax,'NextPlot','add');

for iCell = 1:nCells
    rr = rrr(:,iCell);
    pp = ppp(:,iCell);
    
	plot(ax(iCell), tt, rr, '-', 'Color',0.5*[1 1 1] ); 
	plot(ax(iCell), tt, pp, 'r-' ); 
    ylabel(ax(iCell),'Response');
    sidetitle(ax(iCell),sprintf('Cell %d',iCell));
end
title(ax(1),'Gray: measured; Red: predicted');
xlabel(ax(end), 'Time (s)');

axes(ax(end));
plot( ii/fs, ss, 'k.' )
set(ax(1:nCells),'xcolor','w');
set(ax,'xlim',[min(tt) max(tt)],'ylim',[-inf inf]);
ylabel(ax(end),'Stimulus');

%% Graphics: look at the filters as displaced curves



figure('Color','w'); clf; ax = gridplot(1, nCells); CBs = zeros(nCells,1);
for iCell = 1:nCells
    axes(ax(iCell));
    offset = PlotDisplaced( ShiftTimes, fff(:,:,iCell),4,'std','k' ); hold on
    PlotDisplaced (         ShiftTimes, sss(:,:,iCell), offset, 'absolute', 'b' );
    title(sprintf('Cell %d',iCell));
end
xlabel(ax(1),'Time (s)');
set(ax, 'PlotBoxAspectRatio',[1 3 1]);
supertitle('Black: kernels; Blue: separable approximation');

%% compare "stimulus triggered responses" for data and model

DeltaT = 4; % seconds
nSteps = round(DeltaT*fs);

ccc = nan(2*nSteps+1,nS,nCells); % correlations present in the data
ddd = nan(2*nSteps+1,nS,nCells); % correlations predicted by the model

for iCell = 1:nCells
    for iS = 1:nS
        rr = rrr(:,iCell);
        pp = ppp(:,iCell);
        [ccc(:,iS,iCell),jj] = xcorr( rr, aa(:,iS), nSteps, 'unbiased' );
        [ddd(:,iS,iCell),~ ] = xcorr( pp, aa(:,iS), nSteps, 'unbiased' );
    end
end

figure('Color','w'); clf; ax = gridplot(1, nCells); CBs = zeros(nCells,1);
for iCell = 1:nCells
    axes(ax(iCell));
    offset = PlotDisplaced( jj/fs, ccc(:,:,iCell)',4,'std','k' ); hold on
    PlotDisplaced (         jj/fs, ddd(:,:,iCell)', offset, 'absolute', 'r' );
    title(sprintf('Cell %d',iCell));
end
xlabel(ax(1),'Time (s)');
set(ax, 'PlotBoxAspectRatio',[1 3 1]);
supertitle('Black: stimulus-triggered data; Red: predicted');

%% fit the spatial receptive fields

fitpars = zeros(4,2,nCells);

ii1 = (coords(:,2)==0);
ii2 = (coords(:,1)==0); 

for iCell = 1:nCells
    xx = coords(ii1,1); 
    zz = Tunings{iCell}(ii1);
    [~,fitpars(:,1,iCell)] = fitit('gaussian', zz, ....
        [min(xx), min(zz),         0,  30 ], [], ...
        [max(xx), max(zz), max(zz)/2, 300 ], [0], xx );
    
    xx = coords(ii2,2); 
    zz = Tunings{iCell}(ii2);
    [~,fitpars(:,2,iCell)] = fitit('gaussian', zz, ....
        [min(xx), min(zz),         0,  30 ], [], ...
        [max(xx), max(zz), max(zz)/2, 300 ], [0], xx );
end

%% plot the tuning curves

figure('color','w'); clf; ax = gridplot(nCells,4); set(ax,'NextPlot','add');

for iCell = 1:nCells
    xx = coords(ii1,1);
    many_x = min(xx):max(xx);
    zz = Tunings{iCell}(ii1);
    plot( ax(iCell,1), xx, zz, 'ko' ); % [xtop,ytop,y0,sigma]
    plot( ax(iCell,1), many_x,  gaussian(fitpars(:,1,iCell),many_x));
    yy = coords(ii2,2);
    many_y = min(yy):max(yy);
    zz = Tunings{iCell}(ii2);
    plot( ax(iCell,2), yy, zz, 'ko' ); % [xtop,ytop,y0,sigma]
    plot( ax(iCell,2), many_y,  gaussian(fitpars(:,2,iCell),many_y));
    
    zzz = gaussian(fitpars(:,2,iCell),many_y)'*gaussian(fitpars(:,1,iCell),many_x);
    axes(ax(iCell,3)); imagesc( many_y,many_x,zzz, [0 inf]); 
    set(gca,'plotboxaspectratio',[range(many_x) range(many_y) 1])
    
    plot( ax(iCell,4), ShiftTimes, TimeCourses{iCell});
    sidetitle( ax(iCell,4), sprintf('Cell %d',iCell) );
end
colormap hot
set(ax,'xlim',[-inf inf], 'ylim',[-inf inf],'plotboxaspectratio',[3 2 1]);
set(ax(:,1:2),'ylim',[0 inf])
for iCell = 1:nCells, matchy(ax(iCell,1:2)); end
set(ax(:,2),'ycolor','w');
set(ax(1:end-1,:),'xticklabel',[]);
xlabel(ax(end,1),'x position (deg)');
xlabel(ax(end,2),'y position (deg)');
xlabel(ax(end,4),'Time (s)');
set(ax(:,3),'xtick',[],'ytick',[],'ydir','reverse');
xlabel(ax(end,3),'x pos');
ylabel(ax(end,3),'y pos');




