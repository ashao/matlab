function [yi] = sinterp(y,n2)
%SINTERP  Spline-interpolates a column vector to a new length.
%
%   YI=SINTERP(Y,NI) spline-interpolates the column vector Y to the
%   new length NI.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information (C)
%   1993, 2004 J.M. Lilly --- type 'help jlab_license' for details

n1=length(y);
x=zeros(1,n1);
xi=zeros(1,n2);
x=(1:n1)';
xi=(1:(n1-1)/(n2-1):n1)';
yi=interp1(x,y,xi,'spline');  %This is much better than linear
