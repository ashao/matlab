function[p]=pmat(N,f,dt)
%PMAT Phase shift matrix
%
%   P=PMAT(N,F) returns a matrix an NxN matrix P such that the matrix
%   multiplication Y=H*X for X an Nx1 vector returns Y such that
%
%       Y(j)=X(j) .* exp(2*pi*f*(j - (N+1)/2)*dt);  
%
%   for j = 1, 2, ... N.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004 J.M. Lilly --- type 'help jlab_license' for details        
  
if nargin==2
  dt=1;
end

t=(1:N)-(N+1)/2;
t=t*dt;
p=diag(rot(2*pi*f*t));
