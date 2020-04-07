function [ k ] = smReg( predX, respY, intFlag, lambda, gpuFlag, pDim)
%smReg Linear regression with smoothing regularization
%   For intercept term set intFlag to true.
%   For GPU acceleration set gpuFlag to true.
%   pDim specifies dimensions of predictor (e.g. y by x for an image)

if nargin < 5
    gpuFlag = false;
    pDim = size(predX,2);
elseif nargin < 6
    pDim = size(predX,2);
end

if lambda > 0
    
    [predX,respY] = addLambdaSm(predX, respY, lambda, intFlag, pDim);
    
end

if gpuFlag 
    predX = gpuArray(double(predX));
    respY = gpuArray(double(respY));
    
    %     k = gather(predX' * predX\predX' * respY);
    k = gather(predX\respY);
else
    %     k = predX' * predX\predX' * respY;
    k = predX\respY;
end



end

