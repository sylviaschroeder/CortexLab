function [pupilData, pupilTimes, rawTimes] = loadPupilData(info)

k = strfind(info.expDate, '-');
date = info.expDate;
date(k) = [];
date = str2num(date);

pupilData = [];
pupilTimes = [];

fileProcessed = fullfile('\\zserver.cortexlab.net\Data\EyeCamera', ...
    info.subject, info.expDate, num2str(info.exp), ...
    strrep(info.basename2p, '2P', 'eye_processed.mat'));
if ~exist(fileProcessed, 'file')
    disp('Pupil data not analysed.')
    return
end
data = load(fileProcessed, 'results');
pupilData = data.results;

fileTime = fullfile('\\zserver.cortexlab.net\Data\EyeCamera', ...
    info.subject, info.expDate, num2str(info.exp), ...
    strrep(info.basename2p, '2P', 'eye_TLtime.mat'));
if ~exist(fileTime, 'file')
    fileRaw = fullfile('\\zserver.cortexlab.net\Data\EyeCamera', ...
        info.subject, info.expDate, num2str(info.exp), ...
        strrep(info.basename2p, '2P', 'eye.mat'));
    if ~exist(fileRaw, 'file')
        disp('Pupil data does not exist.')
        return
    end
    pupilTimes = et.getFrameTimes(info.subject, date, info.exp);
    save(fileTime, 'pupilTimes')
else
    data = load(fileTime);
    pupilTimes = data.pupilTimes;
end

if nargout > 2
    warning off
    data = load(fileRaw, 'eyeLog');
    rawTimes = [data.eyeLog.TriggerData.Time];
end