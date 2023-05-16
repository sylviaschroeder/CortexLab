function frameTimes = getFrameTimes2023(timeLine)

nInputs = length(timeLine.hw.inputs);
for iInput = 1:nInputs
    if isequal(timeLine.hw.inputs(iInput).name, 'neuralFrames')
        ind = iInput;
        break;
    end
end

nTotalFrames = timeLine.rawDAQData(end, ind);
countStart = timeLine.rawDAQData(1,ind);
if countStart ~= 0
    warning('frame count in timeLine does NOT start with zero!')
    nTotalFrames = nTotalFrames - countStart;
end
TTLs = [0; diff(timeLine.rawDAQData(:, ind))];
idx = find(TTLs);
idx = idx(1:nTotalFrames);

frameTimes = reshape(timeLine.rawDAQTimestamps(idx), [], 1);