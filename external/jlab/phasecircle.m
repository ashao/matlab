function[y]=phasecircle(phi,x,r,ar)
%PHASECIRCLE Plots circles to indicate a phase angle.
%
%  PHASECIRCLE(PHI) plots a unit radius  "phase circle" with 
%  phase PHI centered at location (0,0).
%  
%  PHASECIRCLE(PHI,X) plots the circle at complex-valued location X.  
%  PHASECIRCLE(PHI,X,R) uses radius R.  
%  PHASECIRCLE(PHI,X,R,AR) rescales the circle for plot with
%  aspect ratio AR=Y/X.  
%
%  Multliple circles are plotted if PHI,X, and R are arrays of 
%  the same size. If PHI is an array but X is not, X is used
%  as a complex-valued offset in between circles, beginning at 0.
%
%  Y=PHASECIRCLE(...) returns a matrix without plotting.
%
%  Usage: phasecircle(phi)
%         phasecircle(phi,x)
%         phasecircle(phi,x,r)
%         phasecircle(phi,x,ar)
%
%   See also ELLIPSEPLOT, CIRC2ELL
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2003,2004 J.M. Lilly --- type 'help jlab_license' for details        

  
if nargin==3
  ar=1;
end

if nargin==2
  ar=1;
  r=1;
end

if nargin==1
  ar=1;
  r=1;
  x=0;
end

if nargin==0
  phi=0;
  ar=1;
  r=1;
  x=0;
end

%/********************************************************
%make things the right size
vcolon(x,phi,r);

N=max([length(r),length(phi),length(x)]);

if length(phi)==1
    phi=phi*ones(N,1);
end

if length(x)==1  &&  N>1
    x=exp(sqrt(-1)*angle(x)).*conj((0:abs(x):(N-1)*abs(x))');    
end

if length(r)==1
    r=r*ones(N,1);
end

%\********************************************************



y=phasecircle1;
x=ones(size(y))*conj(x');
phi=ones(size(y))*conj(phi');
y=(y*(r')).*exp(sqrt(-1).*phi);
y=real(y)+sqrt(-1)*ar*imag(y);   %stretch for aspect ratio
y=y+x;

if nargout==0
  plot(y,'k')
end

function[y]=phasecircle1
z=(0:.1:2*pi+.1)';
y=exp(sqrt(-1)*z);
y=[0+sqrt(-1)*0;y];
