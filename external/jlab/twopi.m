function[y]=twopi(n)
% TWOPI   Raise (2 pi) to the specified power.
% 
%   TWOPI = 2 pi;  TWOPI(N) = (2 pi)^N;
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2004 J.M. Lilly --- type 'help jlab_license' for details
if(nargin==0)
   n=1;
end
y=(2*pi)^n;
