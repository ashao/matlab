function[bool]=cellaresame(x,y,tol)
%CELLARESAME  Tests whether two cell arrays of arrays are the same.
%
%   BOOL=CELLARESAME(X,Y) where X and Y are cell arrays of N arrays,
%
%       X{1}=X1, X{2}=X2,..., X{N}=XN
%       Y{1}=Y1, Y{2}=Y2,..., Y{N}=YN 
%
%   returns the length N cell array BOOL such that
%  
%      BOOL(I)=ARESAME(XI,YI),
%
%   for I=1,2,...,N.  
%
%   BOOL=CELLREAL(X,Y,TOL) uses tolerance TOL.
%
%   See also ARESAME.
%
%   Usage: bool=cellaresame(x,y);
%          bool=cellaresame(x,y,tol);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2009 J.M. Lilly --- type 'help jlab_license' for details
 

if nargin==2
    tol=[];
end
bool=false(length(x),1);

for i=1:length(x)
    bool(i)=aresame(x{i},y{i},tol);
end
 
%reporttest('CELLARESAME',aresame())
