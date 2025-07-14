function signal_new = interpolateWithNaNs(t, signal, t_new, method)

indNan = any(isnan(signal),2) | isnan(t);
t0 = t(~indNan);
s0 = signal(~indNan,:);
[t1, indUni] = unique(t0);
s1 = s0(indUni,:);
signal_new = interp1(t1, s1, t_new, method);
indNan_new = histcounts(t(indNan), [t_new; t_new(end)+1]);
signal_new(indNan_new>0,:) = NaN;