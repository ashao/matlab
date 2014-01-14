function[y]=smartsmooth(x,f)
%SMARTSMOOTH  Fast light smoothing for large matrices
%
%   SMARTSMOOTH(X) for matrix X returns X smoothed by a 3x3 boxcar.  X
%   is zero-padded at the edges before smoothing.
%   
%   SMARTSMOOTH(X,F) uses 3x3 smoothing matrix F.
%
%   This is designed for fast filtering of very large matrices which
%   would take a painfully long time with CONV.  If you want to smooth
%   over more than a 3x3 box, call SMARTSMOOTH multiple times. For
%   example, SMARTSMOOTH(SMARTSMOOTH(X)) is equivalent to smoothing
%   over a 5x5 boxcar.
%
%   SMARTSMOOTH has the correct behavior for asymmetric F.     
%
%   See also VFILT  
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2003, 2004 J.M. Lilly --- type 'help jlab_license' for details      
  
if strcmp(x, '--t')
  smartsmooth_test; return
end

if nargin==1
  f=ones(3,3)/9;
end

y=zeros(size(x,1)+2,size(x,2)+2);
y(2:end-1,2:end-1)=x;

x=0*y;

for i=-1:1
  for j=-1:1
     x1=y;
     %Need to shift negatively to what you'd think to get the right relation
     x1=vshift(x1,-i,1); 
     x1=vshift(x1,-j,2);
     x=x+x1*f(i+2,j+2);
  end
end

y=x(2:end-1,2:end-1);

function[]=smartsmooth_test
x=zeros(5,5);
x(3,3)=2;
y=smartsmooth(x);
yans=[0 0 0 0 0; 0 1 1 1 0 ; 0 1 1 1 0 ; 0 1 1 1 0 ; 0 0 0 0 0];
yans=yans*2/9;
bool=allall(yans==y);
reporttest('SMARTSMOOTH delta-function centered in 5x5 matrix, boxcar',bool);
%reporttest('SMARTSMOOTH delta-function centered in 5x5 matrix, boxcar',bool);

f=[0 1 2; 3 4 5; 0 2 1];
x=zeros(5,5);
x(3,3)=1;
yans=[0 0 0 0 0; 0 f(1,:) 0; 0 f(2,:) 0; 0 f(3,:) 0; 0 0 0 0 0];
y=smartsmooth(x,f);
bool=allall(yans==y);
reporttest('SMARTSMOOTH delta-function centered in 5x5 matrix, asymmetric filter',bool);

f=[0 1 2; 3 4 5; 0 2 1];
x=zeros(5,5);
x(3,4)=1;
yans=[0 0 0 0 0; 0  0 f(1,:); 0 0 f(2,:); 0 0 f(3,:) ; 0 0 0 0 0];
y=smartsmooth(x,f);
bool=allall(yans==y);
reporttest('SMARTSMOOTH delta-function off-center in 5x5 matrix, asymmetric filter',bool);
