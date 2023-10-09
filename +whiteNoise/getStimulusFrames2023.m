function [stimFrames, stimPosition] = getStimulusFrames2023(protocolPath, ...
    hardwarePath)

data = load(protocolPath);
pars = data.Protocol;
stimFile = str2func(strtok(pars.xfile, '.'));
% load myScreenInfo
data = load(hardwarePath);
myScreenInfo = data.myScreenInfo;
myScreenInfo.windowPtr = NaN;
% call x-file to create stimuli
SS = stimFile(myScreenInfo, pars.pars);
stimFrames = cat(3, SS.ImageTextures{:});

stimPosition = [pars.pars(2), pars.pars(3), pars.pars(4), pars.pars(5)] / 10;