function[g,h]=edgeplot(n,p,t)
% EDGEPLOT  Draws limits of edge-effect region on wavelet transform.
%
%   EDGEPLOT(N,P,T) plots the limits of the edge effects of the
%   wavelet transform on a PCOLOR or CONTOUR plot having an x-axis of
%   T and a y-axis measured in central periods P.  N is an array of
%   length LENGTH(P) containing the approximate widths of the wavelets.
%
%   The edge effect limit is approximated by one-half wavelet width N
%   from the left- and right-hand boundaries.
%
%   [H1,H2]=EDGEPLOT(N,P,T) returns the handles to the left- and right-
%   hand side edge lines. 
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information (C)
%   1993, 2004 J.M. Lilly --- type 'help jlab_license' for details
 
n=row2col(floor(n/2));
p=row2col(p);

hold on

ledgex=t(n);
redgex=t(end)-(t(n)-t(1));
g=plot(ledgex+sqrt(-1)*p,'k--');
h=plot(redgex+sqrt(-1)*p,'k--');
set(g,'LineWidth',2.0)
set(h,'LineWidth',2.0)

