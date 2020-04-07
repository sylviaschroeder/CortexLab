function [R, filt] = get_responses3(F0, ton, ops, varargin)

filt = [];
switch ops.method
    case 'deconv'
        if ~isempty(varargin)
            [R, filt] = deconv(F0, ton, ops,varargin{1});
        else
            [R, filt] = deconv(F0, ton, ops);
        end
    case 'sta'
        R = get_sta(F0, ton, ops);
end
end

function [R, filt] = deconv(F0, ton, ops, varargin)
% initialize filters
if numel(varargin)<1
    filt = [0 0 exp(-[1:1:(ops.ntf-2)]/(ops.ntf/3))];
    filt = filt/sum(abs(filt(:)));
else
    filt = varargin{1};
    ops.niter = 1;
end
[NN, NT] = size(F0);

Nstims = numel(ton);
stims = zeros(Nstims,NT);
for i = 1:Nstims
    stims(i,ton{i}) = 1;
end
if ops.smoothORI
   stims(:, :)= my_conv_circ(stims(:, :)', ops.smoothORI)';    
end

pred2 = zeros(ops.ntf, NT);

F1 = sum(F0,1);

for n = 1:ops.niter
    
    pred = filter(filt', 1, stims')';
        
    X = [pred; ones(1,NT)];
    xtx = (X*X')/NT;
    f0x = F0*X'/NT;
    Lreg = eye(size(xtx,1));
    Lreg(Nstims+1:end, Nstims+1:end) = 0;
    B = f0x/(xtx + ops.lam * Lreg);
    %      vexp = mean(mean((F0 - B *X).^2)/VF;
    %     fprintf('error is %2.5f\n', vexp)
    
    if numel(varargin)<1
        for i = 1:ops.ntf
            for j = 1:Nstims
                pred2(i, ton{j}+i-1) = sum(B(:,j));
            end
        end
        
        X = [pred2; ones(1,NT);];
        xtx = (X*X')/NT;
        f0x = F1*X'/NT;
        B2 = f0x/xtx;
        filt = B2(1:ops.ntf);
        filt = filt/sum(abs(filt));
    end
    %     vexp = mean((F0 - B2 *X).^2)/VF;
    
%     fprintf('error is %2.5f\n', vexp)
end


R = B(1:NN,1:Nstims);

end

function R = get_sta(F0, ton)

dt = 0:1:30;
Nstims = length(ton);
[NN , ~] = size(F0);
R = zeros(NN, Nstims);

for i = 1:Nstims
    ton{i}(ton{i}>size(F0,2)-60) = [];
   tresp = repmat(ton{i}, numel(dt), 1) + repmat(dt', 1, numel(ton{i}));
   Fresp = reshape(F0(:, tresp(:)), NN, numel(dt),[]);
   R(:,i) = mean(mean(Fresp,2),3);
end

end