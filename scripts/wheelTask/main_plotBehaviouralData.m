subject = 'SS091';
dates = 'last'; %{'last','2018-01-18','2018-01-19','2018-01-20','2018-01-21','2018-01-22','2018-01-23','2018-01-24','2018-01-25'}
[expRef, expDate] = dat.listExps(subject);
% convert the date(s) to vector of datenums
if ~isa(dates,'numeric')&&~isa(dates,'cell')
    switch lower(dates(1:3))
        case 'las' % last x sessions
            if numel(dates)>4; dates = expDate(end-str2double(dates(5:end))+1:end);
            else; dates = expDate(end);
            end
        case 'tod'; dates = floor(now); % today
        case 'yes'; dates = floor(now)-1; % yesterday
        case 'all'; dates = expDate; % all sessions
        otherwise, dates = sort(datenum(dates)); % specific date
    end
elseif isa(dates,'cell')
    dates = sort(datenum(dates));
end
% get date nums between specified range
if numel(dates)==2&&dateRange==1
    dates = sort(dates);
    dates = dates(1):dates(2);
end
if size(dates,2)>size(dates,1)
    dates = dates';
end
% find paths to existing block files
idx = cell2mat(arrayfun(@(x)find(expDate==x),dates,'UniformOutput',0));
if isempty(idx)
    fprintf('No experiments for this date (%s) recorded.\n', datestr(dates))
    return
end
filelist = mapToCell(@(r) dat.expFilePath(r, 'block', 'master'),expRef(idx));
existfiles = cell2mat(mapToCell(@(l) file.exists(l),filelist));
% load block files
allBlocks = cellfun(@load,filelist(existfiles));
allBlocks = {allBlocks.block};
psychoDat = beh.getBehData(allBlocks);
beh.plotChoiceMatrix(psychoDat);