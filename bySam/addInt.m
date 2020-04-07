function [XInt] = addInt(X)
%addInt Adds column of 1's to design matrix
%   See summary

XInt = [ones(size(X,1),1) X];

end

