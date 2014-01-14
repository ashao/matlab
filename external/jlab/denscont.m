function [hcont,hcont2,hlab,hfp]=denscont(depth,arg2,arg3,arg4,arg5)
%DENSCONT  Density contour overlay for oceanographic T/S plots.
%
%   DENSCONT draws contours of potential density referenced to the
%   sea surface ("sigma-theta") on a plot of temperature vs salinity.
%
%   DENSCONT(PRES) references the potential density to pressure
%   level PRES (in units of dbar), eg. DENSCONT(1500) uses sigma-1500.
%
%   DENSCONT(PRES,V1,V2,STY1,STY2), where V1 and V2 are vectors,
%   draws nonlabelled contours with style STY1 at levels specified
%   in V1, and labelled contours with STY2 at the V2 levels.
%
%   See JCONTOUR for more details about input and output arguments.
%
%   Note that DENSCONT requires the "SEAWATER" Matlab toolbox by
%   Phillip Morgan of CSIRO.  
%
%   DENSCONT also plots the freezing point of water as a heavy line.
%
%   Example:  axis([34 35 1 10]);denscont(0,0.1,0.5,'b','r');
%
%   Usage:  denscont(pres);
%           denscont(pres,v,sty);
%           denscont(pres,v1,v2,sty);
%           denscont(pres,v1,v2,sty1,sty2);
%           [h,h1,hlab,hfp]=denscont(pres,v1,v2,sty1,sty2);
%  
%   See also JCONTOUR, LINESTYLE
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 1999--2011 J.M. Lilly --- type 'help jlab_license' for details  
  
if strcmp(depth,'--t')
    return
end
if ~ishold
        hold on
end

if nargin==0
	depth=0;
end

ax=axis;tha=ax(3);thb=ax(4);sa=ax(1);sb=ax(2);
n=20;
y=linspace(tha,thb,n)'*ones(1,n);
x1=linspace(sa,sb,n);
x=ones(n,1)*x1;
px0=zeros(size(y));
px=depth*ones(size(y));
sigma0=sw_pden(x,sw_temp(x,y,px,px0),px,px);
z=sigma0-1000;

switch nargout
	case 0
		evalme='';
	case 1
		evalme='hcont=';
	case 2
		evalme='[hcont,hcont2]=';
	case 3
		evalme='[hcont,hcont2,hlab]=';
end

evalme=[evalme 'jcontour(x,y,z'];
for i=2:nargin
	evalme=[evalme ',arg' int2str(i)];
end
evalme=[evalme ');'];
eval(evalme);

fp = sw_fp(x1,depth+0*x1);
%figure,plot(x1,fp)
hfp=plot(x1,fp);linestyle -h hfp 2k




