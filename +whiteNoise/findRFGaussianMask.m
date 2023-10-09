function [fitRFs, RFsigns, MSEs] =  findRFGaussianMask(rf)

% fit Gaussian mask to ON field only, OFF field only, or average of ON and
% OFF fields; determine which fit explains most of the RF structure

fitRFs = NaN(size(rf,1), size(rf,2), size(rf,3), 3);
RFsigns = NaN(1, 3);
MSEs = Inf(1, 3);

% start with separate ON and OFF fields
for sub = 1:3
    if sub < 3
        subRF = rf(:,:,sub);
    else
        subRF = (rf(:,:,1) .* RFsigns(1) + rf(:,:,2) .* RFsigns(2)) ./ 2;
    end
    for sign = [1 -1]
        if sub == 3 && sign == -1 % we already found the optimal signs for ON and OFF fields
            continue
        end
        [~, fitRF] = whiteNoise.fit2dGaussRF(subRF .* sign, false);
        mask = zeros(size(rf));
        if sub < 3
            mask(:,:,sub) = fitRF * sign;
        else
            mask(:,:,1) = fitRF * RFsigns(1);
            mask(:,:,2) = fitRF * RFsigns(2);
        end
        mse = sum((rf(:) - mask(:)) .^ 2) / numel(rf);
        if mse < MSEs(sub)
            fitRFs(:,:,:,sub) = mask;
            MSEs(sub) = mse;
            RFsigns(sub) = sign;
        end
    end
end