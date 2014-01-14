function[v,psi]=rankineeddy(eta,ro,vo)
%RANKINEEDDY Velocity and streamfunction for a Rankine vortex.
%
%   RANKINEEDDY computes the velocity and streamfunction observed at
%   the origin due to a Rankine vortex located at a specified point,
%   as described in Lilly and Rhines (2002).
%
%   [v,psi]=rankineeddy(eta,ro,vo);
%
%     Input
%	eta: complex-valued position of eddy (km)
%	ro:  size of eddy (km)
%	vo:  velocity at edge, negative for anticyclone (cm/s) 
%
%     Output
%	v,psi:	azimuthal velocity and streamfunction at origin
%
%   'rankineeddy --f' generates a sample figure
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001, 2004 J.M. Lilly --- type 'help jlab_license' for details    

if strcmp(eta,'--f')
  rankineeddy_fig; return
end

eta=eta*1000;
ro=ro*1000;
vo=vo/100;

warning('off','MATLAB:divideByZero')

r=abs(eta);

if all(size(ro)==1)
	ro=ro.*ones(size(r));
end

vI=vo./ro.*r;
vO=vo.*ro./r;

iO=find(r>ro);
v=vI;
v(iO)=vO(iO);
v=-1*v*sqrt(-1).*exp(sqrt(-1)*angle(eta))*100;

if nargout>1
	psiI_T=(1/4).*(r.^2);
	psiO_T=(1/4).*(ro.^2).*(2*log(r./ro)+1);

	psi=psiI_T;
	psi(iO)=psiO_T(iO);
	psi=psi*vo;
end

function[]=rankineeddy_fig
ro=10;
vo=10;
eta=sqrt(-1)*5+(-50:.1:50)';
[v,psi]=rankineeddy(eta,ro,vo);

figure,
subplot(121),
uvplot(real(eta),v),
title('Cyclone sliced south of center, halfway to edge')
subplot(122)
polar(angle(v),abs(v))
title('Hodograph of Rankine eddy')
