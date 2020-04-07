function [EstAmps, EstPulse, baseline, prediction, numCosines, adjR2] = ...
    SumOfPulses(Sig, PulseOnsets, PulseLen, basisLen, doPlot, PulseLenEssential)
% function [EstAmps, EstPulse] = SumOfPulses(Sig, PulseOnsets, PulseLen, Pulse, RealAmps)
% splits the signal Sig into a sum of pulses of a single estimated shape EstPulse
% of user-specified length PulseLen, occuring at user-specified times PulseOnsets, 
% with estimated amplitudes EstAmps.
% StimID is just for plotting purposes, and the last two arguments are optional, only for if you test it on synthetic data

if nargin < 5 || isempty(doPlot)
    doPlot = 0;
end
if nargin < 6
    PulseLenEssential = PulseLen;
end

if mod(basisLen, 2) == 1
    basisLen = basisLen+1;
end

if doPlot == 1
    fig = figure('Position', [1921 1 1920 1123]);
end

EstAmps = [];
EstPulse = [];
baseline = [];

Sig = Sig(:);
SigLen = length(Sig);

PulseOnsets = PulseOnsets(:);
nPulses = length(PulseOnsets);

BreakTol = 0.01;
MaxIter = 500;

% The algorithm works by iterating two steps: computing the amplitudes
% then computing the pulse shape. Both of these are done by a backslash
% operation of a sparse matrix. To do this quickly, we first make up some 
% matrices that will help build it.

% for basis functions (cosines), which model the continuous baseline offset
% cosine = (1 - cos(linspace(0, 2*pi, basisLen))) / 2;
if basisLen > 0
    cosine = sin(linspace(0, pi, basisLen)).^2;
    numCosines = ceil(SigLen/basisLen*2)+2;
    xb = repmat(1:numCosines, basisLen, 1);
    yb = repmat((-1:numCosines-2)*basisLen/2+1, basisLen, 1) + ...
        repmat((0:basisLen-1)', 1, numCosines);
    cosines = repmat(cosine(:), numCosines, 1);
    bUseMe = find(yb>0 & yb<=SigLen);
    usedCols = ismember(sub2ind(size(yb),repmat(basisLen/2+1,1,size(yb,2)), ...
        1:size(yb,2)), bUseMe);
    matb = sparse(yb(bUseMe), xb(bUseMe), cosines(bUseMe), SigLen, max(xb(bUseMe)));

    if PulseLen < 1 % model the signal only with basis functions
        a = sparse(yb(bUseMe), xb(bUseMe), cosines(bUseMe), SigLen, max(xb(bUseMe)));
        b = a\Sig;
        baseline = a * b;
        prediction = baseline;
        adjR2 = 1 - (sum((Sig-baseline).^2) / (SigLen-numCosines)) / ...
            (sum((Sig-mean(Sig)).^2) / (SigLen-1));
        
        if doPlot == 1
            figure(fig)
            subplot(3,1,2), cla
            plot(yb(basisLen/2+1,usedCols), b(1:sum(usedCols)), 'b')
            xlabel('Samples')
            legend('baseline');
            title('Estimated amplitudes of baseline')
            ax1=gca;
            subplot(3,1,3); cla
            plot(Sig, 'k'); hold on
            plot(baseline, 'b');
            plot([PulseOnsets PulseOnsets]', repmat(ylim,length(PulseOnsets),1)', 'c:');
            legend('real', 'reconstructed');
            title(sprintf('Original and fitted signals (adj. R2: %.2f)',adjR2))
            ax2=gca;
            linkaxes([ax1 ax2],'x');
        end
        return
    end
else
    matb = [];
    numCosines = 0;
    
    if PulseLen < 1 % model the signal only with basis functions
        prediction = [];
        
        if doPlot == 1
            figure(fig)
            subplot(3,1,3)
            plot(Sig, 'k'); hold on
            plot([PulseOnsets PulseOnsets]', repmat(ylim,length(PulseOnsets),1)', 'c:');
            legend('real');
            title('Original signals')
            xlabel('Samples')
        end
        return
    end
end

% for the Shifted Pulse matrix, where Sig ~ ShiftedPulses*PulseAmps
xs = repmat(1:nPulses,PulseLenEssential,1);
ys = repmat((0:PulseLenEssential-1)',1,nPulses) + ...
    repmat(PulseOnsets(:)',PulseLenEssential,1);
sUseMe = find(ys<=SigLen); % this makes sure we don't write outside the matrix

% use with GPU
% Sig_gpu = gpuArray(Sig);
% Sig_gpu = full(Sig_gpu);

indNaN = isnan(Sig);
SigNoNaN = Sig(~indNaN);
ss = ones(nPulses, PulseLenEssential);
ShiftedPulses = [sparse(ys(sUseMe), xs(sUseMe), ss(sUseMe), ...
    SigLen, nPulses), matb];
ShiftedPulses(indNaN,:) = [];
indValid = sum(ShiftedPulses,1) >= .7*PulseLenEssential;
mPulses = sum(indValid);

% initial guess for the amplitudes: all ones.
EstAmps = ones(mPulses,1);

% for the SigLen by nPulses Toeplitz Matrix, where Sig ~ Toep*PulseShape
xt = repmat(1:PulseLen,mPulses,1);
yt = repmat(PulseOnsets(indValid),1,PulseLen) + repmat(0:(PulseLen-1),mPulses,1);
tUseMe = find(yt<=SigLen); % this makes sure we don't write outside the matrix

% for the Shifted Pulse matrix, where Sig ~ ShiftedPulses*PulseAmps
xs = repmat(1:mPulses,PulseLen,1);
ys = repmat((0:PulseLen-1)',1,mPulses) + repmat(PulseOnsets(indValid)',PulseLen,1);
sUseMe = find(ys<=SigLen); % this makes sure we don't write outside the matrix

for i=1:MaxIter
    % make the toeplitz matrix based on current estimated amplitudes
    st = repmat(EstAmps,1,PulseLen);
%     Toep = sparse([yt(tUseMe); yb(bUseMe)], [xt(tUseMe); xb(bUseMe)+PulseLen], ...
%         [st(tUseMe); cosines(bUseMe)], SigLen, PulseLen+max(xb(bUseMe)));
    Toep = [sparse(yt(tUseMe), xt(tUseMe), st(tUseMe), SigLen, PulseLen), ...
        matb];
%     Toep = gpuArray(Toep);
%     Toep = full(Toep);
    Toep(indNaN,:) = [];

    % simple linear regression
    b = Toep\SigNoNaN;
%     b = Toep\Sig_gpu;
%     b = gather(b);
    EstPulse = b(1:PulseLen);
    
    % normalize to have peak one
    EstPulse = EstPulse / max(abs(EstPulse));
    % make positive going
    if max(EstPulse) < abs(min(EstPulse))
        EstPulse = -EstPulse;
    end
    
    % construct a matrix of shifted pulses
    ss = repmat(EstPulse,mPulses,1);
%     ShiftedPulses = sparse([ys(sUseMe); yb(bUseMe)], [xs(sUseMe); xb(bUseMe)+nPulses], ...
%         [ss(sUseMe); cosines(bUseMe)], SigLen, nPulses+max(xb(bUseMe)));
    ShiftedPulses = [sparse(ys(sUseMe), xs(sUseMe), ss(sUseMe), ...
        SigLen, mPulses), matb];
%     ShiftedPulses = gpuArray(ShiftedPulses);
%     ShiftedPulses = full(ShiftedPulses);
    
    % now estimate the amplitude of each pulse
    OldEstAmps = EstAmps;
    
    % simple linear regression
    ShiftedPulses(indNaN,:) = [];
    b = ShiftedPulses\SigNoNaN;
%     b = ShiftedPulses\Sig_gpu;
    EstAmps = b(1:mPulses);
    baselines = b(mPulses+1:end);
    baseline = ShiftedPulses(:,mPulses+1:end) * baselines;
%     EstAmps = gather(EstAmps);
    
    % break if tolerance achieved
    if norm(EstAmps-OldEstAmps)<BreakTol; break; end
end
amps = NaN(nPulses,1);
amps(indValid) = EstAmps;
EstAmps = amps;
pred = ShiftedPulses * b;
prediction = NaN(size(Sig));
prediction(~indNaN) = pred;
% prediction = gather(prediction);
% baseline = gather(baseline);
adjR2 = 1 - (sum((SigNoNaN-pred).^2) / ...
    (SigLen-sum(indNaN)-PulseLen-mPulses-numCosines)) / ...
    (sum((SigNoNaN-mean(SigNoNaN)).^2) / (SigLen-sum(indNaN)-1));
if basisLen <= 0
    baseline = [];
else
    bl = NaN(size(Sig));
    bl(~indNaN) = baseline;
    baseline = bl;
end

if doPlot == 1
    figure(fig)
    % plots some stuff....
    subplot(3,2,1)
    plot(EstPulse)
    xlim([0 length(EstPulse)])
    title('response waveform');
    xlabel('Samples')
    subplot(3,1,2), cla
    stem(PulseOnsets, EstAmps, 'r', 'filled')
    label = {'stim resp'};
    if basisLen > 0
        hold on
        plot(yb(basisLen/2+1,usedCols), baselines(1:sum(usedCols)), 'b')
        %     plot(yb(basisLen/2+1,usedCols), b(1:sum(usedCols)), 'g')
        label = {'stim resp','baseline'};
    end
    legend(label);
    title('Estimated amplitudes')
    ax1=gca;
    subplot(3,1,3); cla
    plot(Sig, 'k'); hold on
    plot(prediction, 'r');
    label = {'original', 'reconstructed'};
    if basisLen > 0
        plot(baseline, 'b');
        label = {'real', 'reconstructed', 'baseline'};
    end
    plot([PulseOnsets PulseOnsets]', repmat(ylim,length(PulseOnsets),1)', 'c:');
    legend(label);
    xlabel('Samples')
    title(sprintf('Original and fitted signals (adj. R2: %.2f)',adjR2))
    ax2=gca;
    linkaxes([ax1 ax2],'x');
    xlim([0 length(Sig)])
end