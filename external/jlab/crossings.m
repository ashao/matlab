function[varargout]=crossings(varargin)
%CROSSINGS  Find values of an array crossing a specified threshold.
%
%   INDEX=CROSSINGS(C,X) where X is a column vector finds the points
%   in X that cross the value of C either from above or below.
%
%   That is, moving down the rows of X, a "crossing point" is a point
%   which exceeds C while its immediate predecessor did not, or 
%   which does not exceed C while its immediate predecessor did.
%
%   INDEX is an index into the locations of these crossing points.
%
%   X may also be an ND array, in which case the crossing points
%   are still found for changes along rows. 
%
%   [I1,I2,...,IN]=CROSSING(C,X1,X2,... XN) finds the crossing
%   points of all input arrays.
%
%   See also TURNINGPOINTS, IND2SUB.
%  
%   'crossings --t' runs a test.
%
%   Usage: index=crossings(c,x);
%          [i1,i2,i3]=crossings(c,x1,x2,x3);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    crossings_test,return
end

c=varargin{1};
varargin=varargin(2:end);

for i=1:nargin-1
    varargout{i}=crossings1(c,varargin{i});
end


function[index]=crossings1(c,x)
x=x(:);

bool1=vshift(x,1,1)>c & x<c;
bool2=vshift(x,1,1)<c & x>c;
bool=bool1|bool2;
bool(1,:)=0;
bool(end,:)=0;

index=find(bool);

function[]=crossings_test
c=0.5;
x=[-1 -1 0 1 3 3 .4 0 -1]';
b=[ 0  0 1 0 0 1  0 0  0]';

reporttest('CROSSINGS',aresame(crossings(c,x),find(b)))

c=0.5;
x=[-1 -1 0 1 3 3 .4 0 -1]';
b=[ 0  0 1 0 0 1  0 0  0]';

[y1,y2]=crossings(c,x,x);
reporttest('CROSSINGS two input arguments',aresame(find(b),y1) && aresame(find(b),y2))

