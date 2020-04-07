function [stimFrames, stimPosition] = getStimulusFrames(info)

SetDefaultDirs;

d = sscanf(info.expDate, '%d-%d-%d');
if d(1)>2017 || (d(1)==2017 && d(2)>=9)
    date = info.expDate;
else
    date = str2double(info.expDate(info.expDate ~= '-'));
end

pars = ProtocolLoad(info.subject, date, info.exp);
stimFile = str2func(strtok(pars.xfile, '.'));
% load myScreenInfo
if d(1)>2017 || (d(1)==2017 && d(2)>9)
    load(fullfile('\\zserver.cortexlab.net\Data\Subjects', info.subject, ...
        date, num2str(info.exp), strrep(info.basenameTL, 'Timeline', 'hardwareInfo')))
else
    load(fullfile(info.folderTL, strrep(info.basenameTL, 'Timeline', 'hardwareInfo')))
end
myScreenInfo.windowPtr = NaN;
% call x-file to create stimuli
SS = stimFile(myScreenInfo, pars.pars);
stimFrames = cat(3, SS.ImageTextures{:});

stimPosition = [pars.pars(2), pars.pars(3), pars.pars(4), pars.pars(5)] / 10;