function[x]=cellangle(x)
%CELLANGLE  Complex argument (angle) of each element in a cell array.
%
%   XA=CELLANGLE(X) where X is a cell array of N arrays,
%
%       X{1}=X1, X{2}=X2,..., X{N}=XN
% 
%   returns the length N cell array of XA with 
%  
%      XA(1)=ANGLE(X1), XA(2)=ANGLE(X2),..., XA(N)=ANGLE(XN).
%
%   See also JCELL.
% 
%   Usage: xa=cellangle(x);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009 J.M. Lilly --- type 'help jlab_license' for details
 
if ~iscell(x)
    error('X must be a cell array.')
end
 
for i=1:length(x)
    x{i}=angle(x{i});
end
