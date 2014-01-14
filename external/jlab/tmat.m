function[x]=tmat()
%TMAT   2x2 complex grouping matrix.
%
%   T=TMAT returns the complex grouping matrix
%
%       T=[1  i; 
%          1 -i] / SQRT(2);
%
%   such that Y=T*X where X=[U V]^T leads Y=[U+iV; U-iV].  Note that
%   T*T'=EYE(2), which implies T'*Y=X.    
%  
%   T=TMAT(N) returns a 3-D array of N copies of the complex grouping
%   matrix, having dimension 2 x 2 x N.
%
%   See also VECTMULT, JMAT, IMAT, and KMAT.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004 J.M. Lilly --- type 'help jlab_license' for details        

if nargin==0
  N=1;
end

x=ones(2,2,N);
x(1,2,:)=sqrt(-1);
x(2,2,:)=-sqrt(-1);

x=frac(1,sqrt(2)).*x;
