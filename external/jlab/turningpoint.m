function[bool]=turningpoint(x)
%TURNINGPOINT  True for turning points, i.e. local extrema, along rows.
%
%   BOOL=TURNINGPOINT(X) where X is a column vector returns BOOL, a
%   vector of the same size as X containing zeros and ones.  BOOL
%   is equal to one if the corresponding element of X is a local
%   maximum or minimum, and zero otherwise.
%
%   X may also be an array of any dimensionality, in which case 
%   local maxima and minima are found with respect to differentiation 
%   along the first dimension, i.e. along rows.
%
%   Note that if any two adjacent points are identical, TURNINGPOINT
%   adds a very small amount of numerical noise to one of the points
%   in order to make a choise about which is larger.
%
%   See also CROSSINGS.
%
%   'turningpoint --t' runs a test.
%
%   Usage: bool=turningpoint(x);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(x, '--t')
    turningpoint_test,return
end

%Add very small amout of noise to keep identical points from occurring
index=find(x==vshift(x,-1,1));
if ~isempty(index)
    x(index)=x(index)+randn(size(x(index)))*1e-10;
end

boolmin=~(x>vshift(x,1,1)|x>vshift(x,-1,1));
boolmax=x>vshift(x,1,1)&x>vshift(x,-1,1);
bool=boolmin|boolmax;
bool(1,:)=0;
bool(end,:)=0;

function[]=turningpoint_test

x=[ 4 5 6 7 6 5 6]';
b=[ 0 0 0 1 0 1 0]';

reporttest('TURNINGPOINT',aresame(b,turningpoint(x)))

x= [ 4 5 6 7 7 6 5 6]';
b1=[ 0 0 0 1 0 0 1 0]';
b2=[ 0 0 0 0 1 0 1 0]';

y1=turningpoint(x);
reporttest('TURNINGPOINT with repeated entry',aresame(y1,b1) || aresame(y1,b2))

