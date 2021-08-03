function colormap = getBlueWhiteRedMap(n)

colormap = ones(n,3);

if mod(n,2) == 0 % even length
    k = n / 2;
    gains = linspace(0, 1, 2*k)';
    gains(1:2:end) = [];
else % odd length
    k = (n-1) / 2;
    gains = linspace(0, 1, k+1)';
    gains(1) = [];
end

% blue half
colormap(1:k,:) = repmat([0 0 1],k,1) .* flip(gains) + ...
    ones(k,3) .* (1 - flip(gains));

% red half
colormap(end-k+1:end,:) = repmat([1 0 0],k,1) .* gains + ...
    ones(k,3) .* (1 - gains);