function[v,psi]=gaussianeddy(eta,ro,vo)
%GAUSSIANEDDY Velocity and streamfunction for a Gaussian vortex.
%
%   GAUSSIANEDDY computes the velocity and streamfunction observed at
%   the origin due to an eddy having a Gaussian streamfunction profile
%   located at a specified point.
%
%   [v,psi]=gaussianeddy(eta,ro,vo);
%
%     Input
%	eta: complex-valued position of eddy (km)
%	ro:  size of eddy (km)
%	vo:  velocity at edge, negative for anticyclone (cm/s) 
%
%     Output
%	v,psi:	azimuthal velocity and streamfunction at origin
%
%   'gaussianeddy --f' generates a sample figure
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001, 2004 J.M. Lilly --- type 'help jlab_license' for details    
  
if strcmp(eta,'--f')
  gaussian_fig; return
end

eta=eta*1000;
ro=ro*1000;
vo=vo/100;

psi=-1*vo*ro*exp(1/2*(1-abs(eta).^2./ro.^2));
v=-1*eta*(vo./ro)*sqrt(-1).*exp(1/2*(1-abs(eta).^2./ro.^2))*100;


function[]=gaussian_fig
ro=10;
vo=10;
eta=sqrt(-1)*5+(-50:.1:50)';

[v,psi]=gaussianeddy(eta,ro,vo);

figure,
subplot(121),
uvplot(real(eta),v),
title('Cyclone sliced south of center, halfway to edge')
subplot(122)
polar(angle(v),abs(v))
title('Hodograph of Gaussian eddy')
