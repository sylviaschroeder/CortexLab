function [psychoDat, allBlocks, allParams] = runPsychometric(mouseName, date, makeFig, varargin)

if ~isempty(varargin)
    blocksToInclude = varargin{1};
else
    blocksToInclude = 1:10000;
end

% psychoDatDir = ['\\zserver2\Data\2Photon\TL_RemoteRepository\' mouseName '\'];
% zsd = getZserverDir();
psychoDatDir = ['\\zserver\Data\expInfo\' mouseName '\'];
% psychoDatDir = fullfile(zsd, 'Data', 'expInfo', mouseName);
if makeFig 
    figDestination = ['\\zserver\Data\behavior\ChoiceWorld\' mouseName '\'];
%     figDestination = fullfile(zsd, 'Data', 'behavior', 'ChoiceWorld', mouseName);
else
    figDestination = 'none';
end

% d = dir([psychoDatDir date '*']);
% expRef = {d.name};
% 
% if ~iscell(expRef)
%     expRef = {expRef};
% end
% 
% for b = 1:length(expRef)    
%     load(fullfile(psychoDatDir, expRef{b}, [expRef{b} '_Block.mat']));
%     if b==1
%         allBlocks = block;
%     else
%         allBlocks(b) = block;
%     end
% end

blockDirs = dir(fullfile(psychoDatDir, date, '*'));
if isempty(blockDirs)
    return
end
b = 1;
for d = 1:length(blockDirs)
    if ~isempty(str2double(blockDirs(d).name)) && ismember(str2double(blockDirs(d).name), blocksToInclude)
        blockFile = dir(fullfile(psychoDatDir, date, blockDirs(d).name, '*Block.mat'));
        if ~isempty(blockFile)
            load(fullfile(psychoDatDir, date, blockDirs(d).name, blockFile.name));
            if b==1 && ~strcmp(block.endStatus, 'exception')
                allBlocks = block;
                b = b+1;
            elseif ~strcmp(block.endStatus, 'exception') && isfield(block, 'expType') && strcmp(block.expType, allBlocks(end).expType)
                allBlocks(b) = block;
                b = b+1;
            end
            
        end
        paramFile = dir(fullfile(psychoDatDir, date, blockDirs(d).name, '*parameters.mat'));
        if ~isempty(paramFile) && ~isempty(blockFile)
            load(fullfile(psychoDatDir, date, blockDirs(d).name, paramFile.name));
            if b==2 && ~strcmp(block.endStatus, 'exception')&& isfield(block, 'expType') && strcmp(block.expType, allBlocks(end).expType) % actually the first time through
                allParams = parameters;
            elseif ~strcmp(block.endStatus, 'exception') && isfield(block, 'expType') && strcmp(block.expType, allBlocks(end).expType)
                allParams(b-1) = parameters;
            end
            
        end
    end
end


% Determine what KIND of experiment this was: 2AFC, 2ADC, 4ADC, or
% 2-stim-2AFC


vcc = allParams(end).visCueContrast;
numRep = allParams(end).numRepeats;
paramsAlt = allParams(end).targetAltitude;
if size(paramsAlt,2)>1
    altsUsed = paramsAlt(numRep>0);
else
    altsUsed = paramsAlt;
end

if sum(min(vcc)>0 & numRep>0)>0
    % if there are any instances with stimuli on left *and* right, it was
    % 2-stim-2AFC
    psychoDat = generatePsychometric2Stim(allBlocks, figDestination);
elseif length(unique(altsUsed))>1
    psychoDat = generatePsychometric4Stim(allBlocks, figDestination);
else

    psychoDat = generatePsychometric(allBlocks, figDestination);
end