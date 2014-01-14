function[xg]=wavegrid(dk,N)
%WAVEGRID  Makes a complex-valued valued (x+iy) grid.  
%
%   WAVEGRID(DK,N) returns an NxN grid with complex-valued elements
%   (X+iY).  X and Y are both linearly spaced with spacing DK.  
%  
%   The maximum absolute values of X and Y will be DK*(N-1)/2.  For
%   DK==1 the grid values will occur at the integers for odd N and
%   halfway between integers for even N.  Thus the grid will include
%   (0,0) if N is odd.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004 J.M. Lilly --- type 'help jlab_license' for details  
  
%   
if strcmp(dk, '--t'),return,end

b=dk*(N-1)/2;  
x=linspace(-b,b,N)';
xg=osum(sqrt(-1)*x,x);  

function[xg]=wavegrid_test
%Test is not working right now
N=5;
dk=1;
K=wavegrid(dk,N);
[m,n]=k2sub(K,dk,N);
K2=sub2k(m,n,dk,N);
all(K(:)==K2(:))
all(m(2,:)==2) && all(n(2,:)==(1:N));

N=4;
dk=1;
K=wavegrid(dk,N);
[m,n]=k2sub(K,dk,N);
K2=sub2k(m,n,dk,N);
all(K(:)==K2(:));
all(m(2,:)==2) && all(n(2,:)==(1:N));

N=5;
dk=pi;
K=wavegrid(dk,N);
[m,n]=k2sub(K,dk,N);
K2=sub2k(m,n,dk,N);
all(K(:)==K2(:));
all(m(2,:)==2) && all(n(2,:)==(1:N));

N=4;
dk=pi;
K=wavegrid(dk,N);
[m,n]=k2sub(K,dk,N);
K2=sub2k(m,n,dk,N);
all(K(:)==K2(:));
all(m(2,:)==2) && all(n(2,:)==(1:N));
