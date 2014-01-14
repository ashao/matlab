function[x]=dmat(N)
%DMAT  2-D unit reflection matrix.
%
%   D=DMAT returns the unit reflection matrix
%
%       D=[0 1; 
%          1 0]
%
%   such that D*X exchanges the first and second elements of the column
%   vector X, i.e. X is reflected across the diagonal line y=x.
%
%   D=DMAT(N) returns a 3-D array of N copies of the unit relection
%   matrix, having dimension 2 x 2 x N.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004 J.M. Lilly --- type 'help jlab_license' for details        
  
if nargin==0
  N=1;
end

x=zeros(2,2,N);
x(1,2,:)=1;
x(2,1,:)=1;

