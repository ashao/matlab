function[v,psi,v2,psi2]=twolayereddy(eta,ro,vo,F1,F2)
%TWOLAYEREDDY Velocity and streamfunction for a 2-layer Rankine vortex.
%
%   TWOLAYEREDDY computes the velocity and streamfunction for a two-
%   layer Rankine vortex consisting of a circular potential vorticity 
%   anomaly in the upper layer only; see Lilly (2002) {Ph. D. thesis,
%   Appendix B.2} for solution details.
%  
%   [v,psi,v2,psi2]=two_layer(eta,ro,vo,F1,F2);
%
%   Input
%	eta: complex-valued position of eddy (km)
%	ro:  size of eddy (km)
%	vo:  upper layer velocity at edge; vo < 0 for anticyclone (cm/s) 
%	F1,F2:	constants characterizing upper and lower layers
%         
%   Output
%	v,psi:	   upper layer azimuthal velocity and streamfunction 
%	v2, psi2:  same, for lower layer
%
%   Note lambda^2 = 1./(F1+F2), where lambda is the Rossby radius.
%
%   'twolayereddy --f' generates a sample figure
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001, 2004 J.M. Lilly --- type 'help jlab_license' for details    
 
if strcmp(eta,'--t'),return,end

if strcmp(eta,'--f')
  twolayereddy_fig; return
end

r1=abs(eta);
ang=exp(sqrt(-1)*(angle(eta)-pi/2));  %JML 03-Mar-2004 note change
                                      %in convention
if nargout<=2
   [v,psi]=twolayereddy_r(r1,ro,F1,F2);
   v=v.*ang.*vo;
   psi=psi*vo;
else
   [v,psi,v2,psi2]=twolayereddy_r(r1,ro,F1,F2);
   v=v.*ang.*vo;
   v2=v2.*ang.*vo;
   psi=psi*vo;
   psi2=psi2*vo;
end

function[]=twolayereddy_fig
ro=10;
vo=10;
eta=sqrt(-1)*5+(-50:.1:50)';

lambda=ro*2;
F1=1/lambda^2;
F2=0;
[v,psi]=twolayereddy(eta,ro,vo,F1,F2);

figure,
subplot(121),
uvplot(real(eta),v),
title('Cyclone sliced south of center, halfway to edge')
subplot(122)
polar(angle(v),abs(v))
title('Hodograph of two-layer eddy, upper layer')



function[v,psi,v2,psi2,scale]=twolayereddy_r(r,R,F1,F2)
%TWOLAYEREDDY_R  Velocity magnitude and streamfunction for a
%2-layer Rankine vortex, as  function of distance from center.
%
%     Input 
%       r: distance to center of eddy (real-valued) 
%       R: size of eddy, may be a matrix 
%       F1,F2: constants characterizing upper and lower layers
%
%     Output
%	v,psi:	azimuthal velocity and streamfunction in upper layer
%	v2, psi2: same, for lower layer
%	
%       All output variables are normalized so max upper layer flow
%	(at r=R) is unity

warning('off','MATLAB:divideByZero')

F=F1+F2;
sqrtF=sqrt(F);

if all(size(R)==1)
	R=R.*ones(size(r));
end

ko=besselk(0,R.*sqrtF);
io=besseli(0,R.*sqrtF);
k1=besselk(1,R.*sqrtF);
i1=besseli(1,R.*sqrtF);
ko_norm=besselk(0,r*sqrtF)./ko;
io_norm=besseli(0,r*sqrtF)./io;

k1_norm=besselk(1,r*sqrtF)./ko;
i1_norm=besseli(1,r*sqrtF)./io;


phi=io.*k1./(i1.*ko);


%found an error 3/20/99 sqrtF should divide baroclinic part
vI_C=phi./(phi+1).*i1_norm./sqrtF;
vI_T=(F2./F).*(1/2).*r;
vO_C=1./(phi+1).*k1_norm./sqrtF;
vO_T=(F2./F).*(1/2).*(R.^2)./r;

iO=find(r>R);
v=vI_T + (F1./F).*vI_C;
v(iO)=vO_T(iO) + (F1./F).*vO_C(iO);

scale=(F2./F).*(1/2).*R + (F1./F).*phi./(phi+1).*(i1./io)./sqrtF;
v=v./scale;

if nargout>1
	psiI_C=-1*(1-phi./(phi+1).*io_norm)/F;
	psiI_T=(F2./F).*(1/4).*(r.^2);
	psiO_C=-1./(phi+1).*ko_norm/F;
	psiO_T=(F2./F).*(1/4).*(R.^2).*(2*log(r./R)+1);

	psi=psiI_T + (F1./F).*psiI_C;
	psi(iO)=psiO_T(iO) + (F1./F).*psiO_C(iO);
	psi=psi./scale;
end


if nargout>2
	v2=vI_T - (F2./F).*vI_C;
	v2(iO)=vO_T(iO) - (F2./F).*vO_C(iO);
	v2=v2./scale;

	psi2=psiI_T - (F2./F).*psiI_C;
	psi2(iO)=psiO_T(iO) - (F2./F).*psiO_C(iO);

	psi2=psi2./scale;
end


warning('on','MATLAB:divideByZero')
