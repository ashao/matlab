function[h]=hmat(N,dt)
%HMAT Hilbert transform matrix
%
%  H=HMAT(N) returns a matrix an NxN matrix H such that the matrix
%  multiplication Y=H*X for X another NxN matrix returns a matrix Y,
%  which is the discrete Hilbert transform of X.
%
%  H=HMAT(N,DT) specifies the 'time' step, assumed to be one
%  otherwise, between the elements of the columns of X.  
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004 J.M. Lilly --- type 'help jlab_license' for details        
  
if nargin==1
  dt=1;
end

f=fmat(N,dt);

sgnom=ones(N,1);
sgnom(1)=0;
index=(1:N);
sgnom(index-1>N/2)=-1;

d=sqrt(-1)*diag(sgnom);

h=f'*d*f;


if 0
  h=zeros(N,N);

if iseven(N)
   t=(-(N-2)/2:N/2)';
elseif isodd(N)
   t=(-(N-1)/2:(N-1)/2)';
end

index=(1:length(h(:)))';
[ii,jj]=ind2sub(size(h),index);

h(index)=frac(dt,pi).*frac(1,t(ii)-t(jj));

h(~isfinite(h))=0;
end




