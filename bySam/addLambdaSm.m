function [XLambda,YLambda] = addLambdaSm(X, Y, lambda, intFlag, pDim)
%addLambdaRRSmAdds lambda values to design matrix for smoothing
%regularization
%   Adds lambda values to design matrix X and corresponding zeros to Y as
%   necessary for smoothing regularization. Intercept is added if
%   specified.

smMat = nan(size(X,2));

for p = 1:prod(pDim)
    
    zMat = zeros(pDim);
    
    [r,c] = ind2sub(pDim,p);
    
    ndim = pDim(1) > 1; % Check if predictor is 2d
    
    switch ndim
        
        case true
                    
            if r == 1 && c == 1                         % Corners
                zMat(r,c) = 2;
                zMat(r+1,c) = -1;
                zMat(r,c+1) = -1;
            elseif r == pDim(1) && c == pDim(2)
                zMat(r,c) = 2;
                zMat(r-1,c) = -1;
                zMat(r,c-1) = -1;
            elseif r == pDim(1) && c == 1
                zMat(r,c) = 2;
                zMat(r-1,c) = -1;
                zMat(r,c+1) = -1;
            elseif r == 1 && c == pDim(2)
                zMat(r,c) = 2;
                zMat(r+1,c) = -1;
                zMat(r,c-1) = -1;
            elseif r == 1 && c > 1 && c < pDim(2)       % Sides
                zMat(r,c) = 3;
                zMat(r,c+[-1 1]) = -1;
                zMat(r+1,c) = -1;
            elseif r == pDim(1) && c < pDim(2) && c > 1
                zMat(r,c) = 3;
                zMat(r,c+[-1 1]) = -1;
                zMat(r-1,c) = -1;
            elseif r > 1 && c == 1 && r < pDim(1)
                zMat(r,c) = 3;
                zMat(r+[-1 1],c) = -1;
                zMat(r,c+1) = -1;
            elseif r < pDim(1) && c == pDim(2) && r > 1
                zMat(r,c) = 3;
                zMat(r+[-1 1],c) = -1;
                zMat(r,c-1) =-1;
            else                                        % Middle
                zMat(r,c) = 4;
                zMat(r+[-1 1],c) = -1;
                zMat(r,c+[-1 1]) = -1;
            end
            
        case false
            
            if c == 1
                zMat(r,c) = 1;
                zMat(r,c+1) = -1;
            elseif c > 1 && c ~= pDim(2)
                zMat(r,c) = 2;
                zMat(r,c+[-1 1]) = -1;
            elseif c == pDim(2)
                zMat(r,c) = 1;
                zMat(r,c-1) = -1;
            end
            
    end
    
    smMat(p,:) = zMat(:)';
    
end

if intFlag
    
    X = addInt(X);
    smMat = [zeros(size(smMat,1),1) smMat];
    
end

XLambda = [X; smMat.*sqrt(lambda)];

YLambda = [Y; zeros(size(smMat,1),size(Y,2))];


end

