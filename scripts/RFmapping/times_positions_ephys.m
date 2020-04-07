subject = 'SS079';
date = '2017-09-01';
exp = 2;
TLexp = 1;
tag = 'K1';
 
protocolFolder = '\\ZSERVER.cortexlab.net\Data\trodes';
hardwareFolder = '\\ZSERVER.cortexlab.net\Data\expInfo';
subjectsFolder = '\\ZSERVER.cortexlab.net\Data\Subjects';
 
data = load(fullfile(protocolFolder, subject, date, num2str(exp), ...
    'Protocol.mat'));
pars = data.Protocol;
stimFile = str2func(strtok(pars.xfile, '.'));
% load myScreenInfo
load(fullfile(hardwareFolder, subject, date, num2str(exp), ...
    sprintf('%s_%d_%s_hardwareInfo.mat', date, exp, subject)));
myScreenInfo.windowPtr = NaN;
 
% call x-file to create stimuli
SS = stimFile(myScreenInfo, pars.pars);
stimFrames = cat(3, SS.ImageTextures{:});
 
framesPerImage = pars.pars(6,1);
frameTimes = (0 : size(stimFrames, 3)-1) * framesPerImage / myScreenInfo.FrameRate;
 
rows = size(stimFrames,1);
cols = size(stimFrames,2);
ind = find(stimFrames == 1);
t = ceil(ind / (rows * cols));
ind = mod(ind, (rows * cols));
ind(ind == 0) = rows * cols;
x = ceil(ind / rows);
y = mod(ind, rows);
y(y == 0) = rows;
time = frameTimes(t)';
 
alignDir = fullfile(subjectsFolder, subject, date, 'alignments');
bTLtoMaster = readNPY(fullfile(alignDir, ...
    sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, tag)));
stimOnTL = readNPY(fullfile(alignDir, ...
    sprintf('mpep_%d_onsets_in_timeline_%d.npy', exp, TLexp)));
stimOffTL = readNPY(fullfile(alignDir, ...
    sprintf('mpep_%d_offsets_in_timeline_%d.npy', exp, TLexp)));
 
stimOn = applyCorrection(stimOnTL, bTLtoMaster);
stimOff = applyCorrection(stimOffTL, bTLtoMaster);
    
% stimPosition = [pars.pars(2), pars.pars(3), pars.pars(4), pars.pars(5)] / 10;

x = repmat(x, 3, 1);
y = repmat(y, 3, 1);
times = reshape(bsxfun(@plus, stimOn', time), [], 1);
