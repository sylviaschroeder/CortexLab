function out = takeRatio(a,b)
    ind = a .* b <= 0;
    a(ind) = NaN;
    b(ind) = NaN;
    out = log10((abs(b)+0.1) ./ (abs(a)+0.1));
end