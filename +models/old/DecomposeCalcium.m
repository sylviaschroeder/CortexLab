function [Templates, ScaledBy] = SumOfPulses(RawTrace, TrialOnset, StimID, TemplateLen)
% [EstAmps, EstPulse] = DecomposeCalcium(Sig, TrialOnsets, StimID, PulseLen, Pulse, RealAmps)
% decomposes the signal RawTrace into a sum of templates, triggered at each
% trial onset time, and scaled by a factor that varies from one trial to
% the next.
% 
% TrialOnset should give the start time of each trial (i.e. stimulus onset) in samples
% StimID should give the stimulus ID on that trial (an integer from 1 to nStims)
% TemplateLen is the length of the template it fits to each stimulus.
% 
% The outputs are Templates, a TemplateLen by nStims array giving the
% template found for each stimulus
% of user-specified length PulseLen, occuring at user-specified times PulseOnsets, 
% with estimated amplitudes EstAmps.

RawTrace = RawTrace(:);
TraceLen = length(RawTrace);

TrialOnset = TrialOnset(:);
StimID = StimID(:);
nPulses = length(TrialOnset);

nStims = max(StimID);

% the last two arguments are for when you test it on synthetic data
if nargin<5; Pulse = zeros(TemplateLen,1); end
if nargin<6; RealAmps = zeros(nPulses, 1); end

% The algorithm works in two steps. In the first, we estimate the pulse
% waveforms for each stimulus. To do this, we use a backslash
% to solve Sig ~ Toep*PulseShapes 
% where Toep is a SigLen by (PulseLen*nStims) Toeplitz-ish Matrix, 
% which consists of nStims SigLen*PulseLen submatrices. Each submatrix
% contains a diagonal stripe of ones for every stimulus of that StimID, 
% starting on the left at the stimulus onset time.

% First we make up some index matrices that will help build the Toeplitz matrix it.
% the (i,j)th element of these matrices gives the coordinates of the point 
% we add for time j into pulse i.
xt = repmat(1:TemplateLen,nPulses,1) + repmat((StimID-1)*TemplateLen, 1, TemplateLen);
yt = repmat(TrialOnset(:),1,TemplateLen) + repmat(0:(TemplateLen-1),nPulses,1);
tUseMe = find(yt<=TraceLen); % this makes sure we don't write outside the matrix

% make it (as a sparse matrix, for speed)
Toep = sparse(yt(tUseMe),xt(tUseMe),1,TraceLen,TemplateLen*nStims);

% compute pulse shapes and reshape into PulseLen by nStims matrix
AllPulses = Toep\RawTrace;
Templates = reshape(AllPulses, TemplateLen, nStims);

% In the second step, we find out how much these pulse waveforms need to be
% multiplied by on each individual stimulus presentation. To do this, we
% make a ShiftedPulse matrix, in which the jth column corresponds to the
% template waveform of the appropriate stimulus, shifted so it starts at
% stimulus onset. So Sig ~ ShiftedPulses*PulseAmps

% Again we first make up some index matrices. This time, the (i,j)th element 
% of these gives the coordinates of the point we add for time i into pulse j
xs = repmat(1:nPulses,TemplateLen,1);
ys = repmat((0:TemplateLen-1)',1,nPulses) + repmat(TrialOnset(:)',TemplateLen,1);
sUseMe = find(ys<=TraceLen); % this makes sure we don't write outside the matrix

% now we make the matrix of values to go in those coordinates. The (i,j)th
% element of this is the value of the template waveform at time i, for the 
% stimulus that was presented on trial j
ss = zeros(TemplateLen, nPulses);
for j=1:nPulses
    ss(:,j) = Templates(:,StimID(j));
end

ShiftedPulses = sparse(ys(sUseMe), xs(sUseMe), ss(sUseMe),TraceLen,nPulses);

% now we compute how much each pulse needs to be multipled by on a given
% trial:

ScaledBy = ShiftedPulses\RawTrace;

% now we see how well the reconstruction worked. To do this, we make a
% modified Toeplitz matrix, which instead of ones has EstMaps for each
% stimulus
st = repmat(ScaledBy,1,TemplateLen);
Toep2 = sparse(yt(tUseMe),xt(tUseMe),st(tUseMe),TraceLen,TemplateLen*nStims);

% estimated signal with and without scaling
EstSig = Toep*AllPulses;
EstSig2 = Toep2*AllPulses;

% Now plot diagnostics
figure(1); clf;
subplot(2,1,1);
imagesc(Templates');
title('Pulse waveforms');
xlabel('Time');
ylabel('Stimulus ID');
colorbar

subplot(2,1,2);
plot(TrialOnset, ScaledBy, '.');
xlabel('Time');
ylabel('Scalilng factor');


figure(2); clf;
plot([EstSig2, EstSig, RawTrace]); hold on
plot([TrialOnset TrialOnset]', repmat(ylim,length(TrialOnset),1)', 'k:');
textypos = min(RawTrace);
for i=1:nPulses
    text(TrialOnset(i), textypos, sprintf('stim %d\n mult %.1f', StimID(i), ScaledBy(i)));
    %num2str(StimID(i)));
end
legend('scaled', 'unscaled', 'original'); 
title('Raw signal');

return;
