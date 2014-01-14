function[x]=cellconj(x)
%CELLCONJ  Complex conjugate of each element in a cell array.
%
%   XC=CELLCONJ(X) where X is a cell array of N arrays,
%
%       X{1}=X1, X{2}=X2,..., X{N}=XN
% 
%   returns the length N cell array XC with 
%  
%      XC(1)=CONJ(X1), XC(2)=CONJ(X2),..., XC(N)=CONJ(XN).
%
%   See also JCELL.
%
%   Usage: xc=cellconj(x);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009 J.M. Lilly --- type 'help jlab_license' for details

if ~iscell(x)
    error('X must be a cell array.')
end

for i=1:length(x)
    x{i}=conj(x{i});
end
