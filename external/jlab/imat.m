function[x]=imat(N)
%IMAT  2x2 identify matrix.
%
%   I=IMAT returns the identify matrix
%
%       D=[1 0; 
%          0 1]
%
%   such that I*X equals X.
%
%   I=IMAT(N) returns a 3-D array of N copies of the identity matrix,
%   having dimension 2 x 2 x N.
%
%   See also VECTMULT, JMAT, KMAT, and TMAT.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004 J.M. Lilly --- type 'help jlab_license' for details        
  
if nargin==0
  N=1;
end

x=zeros(2,2,N);
x(1,1,:)=1;
x(2,2,:)=1;

