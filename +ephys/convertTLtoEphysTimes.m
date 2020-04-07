

function whispTimes = convertTLtoEphysTimes(tlTimes, tlOffset)

if length(tlOffset) == 1
    whispTimes = tlTimes+tlOffset;
    
elseif length(tlOffset)==2 
    % this is a linear regression conversion
    
    if isrow(tlTimes); tlTimes = tlTimes'; end
    if isrow(tlOffset); tlOffset = tlOffset'; end;
    
    whispTimes = [tlTimes ones(size(tlTimes))]*tlOffset;
else
    
    whispTimes = tlTimes;
    disp(' warning: unrecognized tlOffset. No conversion applied.');
end

