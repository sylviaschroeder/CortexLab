function prediction = getStimulusResponsePrediction(stimMatrix, blankStim, ...
    alphas, kernel)

% get stimulus onsets (exclude blank)
nonBlanks = setdiff(1:size(stimMatrix,1), blankStim);
stimOnsetsAll = [0 diff(sum(stimMatrix,1))];
stimOnsetsAll = find(stimOnsetsAll == 1);
stimIDsAll = sum(bsxfun(@times, stimMatrix, (1:size(stimMatrix,1))'),1);
stimIDsAll = stimIDsAll(stimOnsetsAll);
stimOnsets = stimOnsetsAll;
stimOnsets(stimIDsAll == blankStim) = [];
stimIDs = stimIDsAll;
stimIDs(stimIDs==blankStim) = [];

% for the Shifted Pulse matrix, where Sig ~ ShiftedPulses*PulseAmps
nPulses = length(stimOnsets);
PulseLen = length(kernel);
xs = repmat(1:nPulses,PulseLen,1);
ys = repmat((0:PulseLen-1)',1,nPulses) + repmat(stimOnsets(:)',PulseLen,1);
sUseMe = find(ys<=size(stimMatrix,2)); % this makes sure we don't write outside the matrix

% construct a matrix of shifted pulses
ss = repmat(kernel,nPulses,1);
ShiftedPulses = sparse(ys(sUseMe), xs(sUseMe), ss(sUseMe), ...
    size(stimMatrix,2), nPulses);

% construct vector of response amplitudes
reps = repmat(1:(length(stimOnsets)/length(nonBlanks)),length(nonBlanks),1);
reps = reshape(reps,[],1);
inds = sub2ind(size(alphas),reps,stimIDs');
ampls = alphas(inds);

prediction = ShiftedPulses * ampls;