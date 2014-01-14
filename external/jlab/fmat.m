function[f]=fmat(N,dt)
%FMAT  Unitary Fourier transform matrix
%
%   F=FMAT(N) returns a matrix an NxN matrix F such that the matrix
%   multiplication Y=F*X for X another NxN matrix returns a matrix Y
%   which is the unitary Fourier transform of X.  F' performs the
%   inverse Fourier transform.  
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004 J.M. Lilly --- type 'help jlab_license' for details        



if nargin==1
  dt=1;
end

T=N*dt;

f=zeros(N);

index=(1:length(f(:)))';
[ii,jj]=ind2sub(size(f),index);
f(index)=frac(1,sqrt(N)).*exp(-2*pi*sqrt(-1).*(dt./T).*(ii-1).*(jj-1));
 
