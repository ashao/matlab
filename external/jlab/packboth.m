function[h]=packboth(i1,i2,i3)
%PACKBOTH   Squeeze together rows and columns of the current figure.
%  
%   PACKBOTH(M,N) squeezes together all rows and all columns in the
%   current figure, which has M rows and N columns generated with
%   SUBPLOT. This is used when all subplots in a given column share a
%   x- and y-axes.
%  
%   X- and y-axes tickmarks and labels are adjusted as described in
%   PACKROWS and PACKCOLS.
%  
%   H=PACKBOTH(M,N) returns a vector of handles H into the subplots.
%
%   PACKBOTH(H,M,N) also works, where H is a vector of handles.
%
%   After calling PACKBOTH, do not call SUBPLOT again, as this will
%   destroy the subplots; instead, access the subplots through
%   AXES(H(I)) where I is an integer. 
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2007 J.M. Lilly --- type 'help jlab_license' for details   

if nargin==2
  m=i1;
  n=i2;
  h=subplots(m,n);
elseif nargin==3
  h=i1;
  m=i2;
  n=i3;
end

packrows(h,m,n);
packcols(h,m,n);

if nargout==0
  clear h
end

