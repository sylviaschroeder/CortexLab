function matrix = makeLambdaMatrix(stimRows, stimColumns, numFrames)

matrix = eye(stimRows * stimColumns); % one column for each pixel; set that pixel to one, this will be the center pixel for the column
matrix = reshape(matrix, stimRows, stimColumns, []); % reshape matrix into image frames with one center pixel each
tmp = matrix;
matrix(1:end-1,:,:) = matrix(1:end-1,:,:) - tmp(2:end,:,:); % add -1 above each center pixel
matrix(2:end,:,:) = matrix(2:end,:,:) - tmp(1:end-1,:,:); % add -1 below each center pixel
matrix(:,1:end-1,:) = matrix(:,1:end-1,:) - tmp(:,2:end,:); % add -1 to left of each center pixel
matrix(:,2:end,:) = matrix(:,2:end,:) - tmp(:,1:end-1,:); % add -1 to right of each center pixel
matrix = matrix - tmp; % set center pixels to zero
matrix = reshape(matrix, [], size(matrix,3)); % reshape matrix to have one column per center pixel
% normalise so that all smoothing values per parameter sum to -1 (except on edges)
matrix = matrix ./ abs(median(sum(matrix,1)));
% set center pixel to 1
tmp = reshape(tmp, [], size(tmp,3));
matrix = matrix + tmp;
matrix = matrix'; % transpose matrix -> each row for one center pixel
tmp = matrix;
for n = 1:numFrames-1 % repeat matrix to match number of frames in RF
    matrix = blkdiag(matrix, tmp);
end