function[upa,upb,upc,upd,upe]=ellband(varargin)
%ELLBAND  Bandwidth of modulated elliptical signals in two or three dimensions.
%
%   [A,B,C]=ELLBAND(KAPPA,LAMBDA,THETA,PHI) computes the instantaneous 
%   bandwidth of the elliptical signal characterized by RMS amplitude 
%   KAPPA, linearity LAMBDA, orientation THETA, and orbital phase PHI.
%
%   The three output arguments are
%
%          A  - Amplitude modulation bandwidth  
%          B  - Deformation bandwidth
%          C  - Precession bandwidth. 
%
%   and these satisfy UPSILON^2=A^2+B^2+C^2 where UPSILON is the joint 
%   instantaneous bandwidth of the bivariate signal.
%
%   The form of these terms is as follows:
%
%         A = 1/KAPPA d/dt KAPPA 
%         B = 1/2 * 1/SQRT(1-LAMBDA^2) *  d/dt LAMBDA 
%         C = LAMBDA d/dt THETA
%
%   For details see Lilly and Olhede (2010).
%   __________________________________________________________________
%
%   Three dimensions
%
%   ELLBAND can also compute the instantaneous bandwidth of modulated 
%   elliptical signals in three dimensions.
%
%   [A,B,C,D,E]=ELLBAND(KAPPA,LAMBDA,THETA,PHI,ALPHA,BETA) returns the
%   terms in the bandwidth from a modulated ellipical signal in a plane
%   with a normal vector having azimuth angle ALPHA and zenith angle BETA.
%   
%   The five output arguments are
%
%          A   - Amplitude modulation bandwidth, as in 2D 
%          B   - Deformation bandwidth, as in 2D
%          C   - Precession bandwidth, as in 2D
%          D   - Precession bandwidth with full 3D effects
%          E   - Bandwidth due to motion of the normal to the plane
%
%   and these, in principle, satisfy USPILON^2=A^2+B^2+D^2+|E|^2 where 
%   UPSILON is the joint instantaneous bandwidth of the trivariate signal.
%   See below for a caveat on this statement.
%
%   Terms A--C are just as in the bivariate case.  The new terms are:
%
%         D = LAMBDA [d/dt THETA + COS(BETA) * d/dt ALPHA]
%         E = N^T X_+ / |X_+^H X_+|    
%   
%   where N is the trivariate normal vector, X_+ is the trivariate analytic
%   signal vector, and ``T'' denotes the matrix transpose, and ``H'' the 
%   Hermitian transpose.  Note that term E may be complex-valued.
%
%   Note that term C does not contribute to the full bandwidth, but is 
%   output in order to compare the two-dimensional and three-dimensional
%   effects in the full precession bandwidth, term D.
%
%   An important point is that the trivariate ellipse parameters can be 
%   ill-defined for a nearly linear signal, and the elliptical bandwidth 
%   terms can give erroneously large values at isolated points.  To check
%   for this, compare with the bandwidth from JOINTFREQ.
%
%   For details see Lilly (2011).
%   __________________________________________________________________
%
%   ELLBAND(DT,...) sets the sample interval DT, which defaults to DT=1.
%   
%   See also ANATRANS, WAVETRANS, INSTFREQ, JOINTFREQ.
%
%   'ellband --t' runs a test.
%   'ellband --f' generates a figure from Lilly and Olhede (2010).
%
%   Usage:  [a,b,c]=ellband(kappa,lambda,theta,phi);
%           [a,b,c]=ellband(dt,kappa,lambda,theta,phi);  
%           [a,b,c,d,e]=ellband(kappa,lambda,theta,phi,alpha,beta); 
%           [a,b,c,d,e]=ellband(dt,kappa,lambda,theta,phi,alpha,beta);    
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006--2011 J.M. Lilly --- type 'help jlab_license' for details
 

if strcmp(varargin{1}, '--t')
    ellband_test,return
end

if strcmp(varargin{1}, '--f')
    ellband_figure,return
end

if nargin==5||nargin==7
    dt=varargin{1};
    varargin=varargin(2:end);
else
    dt=1;
end

kappa=varargin{1};
lambda=varargin{2};
theta=varargin{3};
phi=varargin{4};
if length(varargin)>4
    alpha=varargin{5};
    beta=varargin{6};
    dims=3;
else
    dims=2;
end

omtheta=frac(vdiff(unwrap(theta),1),dt);
zeta=sign(lambda).*sqrt(1-lambda.^2);

upa=frac((vdiff(log(kappa),1)),dt);
upb=frac(sqrt(frac(1,1-lambda.^2)).*(frac(1,2)*vdiff(lambda,1)),dt);

upb2=frac(sqrt(frac(1,1-zeta.^2)).*(frac(1,2)*vdiff(zeta,1)),dt);
upb(abs(lambda)>sqrt(1/2))=-upb2(abs(lambda)>sqrt(1/2));
%To correct numerical instability when lambda is close to unity

%This is the same
%[a,b]=kl2ab(kappa,lambda);
%upc=abs(frac(a.*b,a.^2+b.^2).*vdiff(log(abs(b./a)),1))./dt;
upc=(lambda.*omtheta);

if dims==3
     omalpha=frac(vdiff(unwrap(alpha),1),dt);
     ombeta=frac(vdiff(unwrap(beta),1),dt);
     
     upd=lambda.*(omtheta+omalpha.*cos(beta));
       
     [x,y,z]=ellsig(kappa,lambda,theta,phi,alpha,beta);
      
     %Now find x-tilde, 2-d vector
     [x,y,z]=vectmult(jmat3(-alpha,3),x,y,z);
     [x,y,z]=vectmult(jmat3(-beta,1),x,y,z);
     
     numer=-omalpha.*sin(beta).*x+ombeta.*y;     
     denom=sqrt(abs(x).^2+abs(y).^2);
     upe=frac(numer,denom); 
else 
    upd=[];
    upe=[];
end


function[]=ellband_test
ellband_test1;
ellband_test2;
ellband_test3;
ellband_test4;


function[]=ellband_test1

load npg2006
use npg2006

%Decide on frequencies
fs=2*pi./(logspace(log10(10),log10(100),50)');

%Compute wavelet transforms using generalized Morse wavelets
[wx,wy]=wavetrans(real(cx),imag(cx),{1,2,4,fs,'bandpass'},'mirror');
[ir,jr,wr,fr]=ridgewalk(dt,wx,wy,fs,{3,0,'amp'});
[kappa,lambda,theta,phi]=ellparams(wr(:,1),wr(:,2));
[ba,bd,bp]=ellband(kappa,lambda,theta,phi);   
[om,upbar]=instfreq(wr,2);


vindex(upbar,ba,bd,bp,2:length(upbar)-1,1);
err=vsum((upbar-sqrt(ba.^2+bd.^2+bp.^2)).^2,1)./vsum(upbar.^2,1);
reporttest('ELLBAND terms sum to bandwidth from INSTFREQ', err<1e-5); 

function[]=ellband_test2
t=(0:1:925)';
cxe=zeros(length(t),3);

kappa=3*exp(2*0.393*(t/1000-1));
lambda=0.4+0*t;
phi=(t/1000*5)*2*pi;
theta=pi/4+0*t;

om=vdiff(phi,1);  %Since theta is constant

[x,y]=ellsig(kappa,lambda,theta,phi);
[om,up]=instfreq([x y],2);
ups=sqrt(vsum(abs(up).^2,2));
[upsa,upsb,upsc]=ellband(kappa,lambda,theta,phi);
b1=aresame(vmean([ups upsa upsb upsc]./[om om om om],1),[1 1 0 0]*0.025,1e-4);

om=vdiff(phi,1);  %Since theta is constant

kappa=2.5+0*t;
lambda=zeros(size(t));
for i=2:length(t)
    lambda(i)=real(lambda(i-1)+2*sqrt(1-lambda(i-1).^2)*0.025.*om(i));
end
lambda(1)=nan;
lambda(lambda>1)=1;  

[x,y]=ellsig(kappa,lambda,theta,phi);
[om,ups]=instfreq([x y],2);
[upsa,upsb,upsc]=ellband(kappa,lambda,theta,phi);
b2=aresame(vmean([ups upsa upsb upsc]./[om om om om],1),[1 0 1 0]*0.025,1e-4);

[kappa,lambda]=ab2kl(3+zeros(size(t)),2+zeros(size(t)));

theta=phi/14.45;

[x,y]=ellsig(kappa,lambda,theta,phi);
[om,ups]=instfreq([x y],2);
[upsa,upsb,upsc]=ellband(kappa,lambda,theta,phi);
b3=aresame(vmean([ups upsa upsb upsc]./[om om om om],1),[1 0 0 1]*0.025,1e-4);

reporttest('ELLBAND three panels of figure each have bandwidth 0.025',b1&&b2&&b3)

function[]=ellband_test3

kappao=10;
lambdao=2/3;
omtheta=linspace(-2,2,100);
omphi=1;

dt=0.1;
t=[0:dt:100]';


zo=ellsig(kappao,lambdao,0,omphi.*t);
psi=sleptap(length(zo)); 

z=vzeros(length(zo),length(omtheta));
for i=1:length(omtheta)
    z(:,i)=rot(omtheta(i).*t).*zo;
end

%Taper for better computation
z=z.*vrep(psi(:,1),size(z,2),2);
[zp,zn]=anatrans(z);
    
[omp,upp]=instfreq(dt,zp);
[omn,upn]=instfreq(dt,zn);

kappa=sqrt(abs(zp).^2+abs(zn).^2);
om=frac(omp.*abs(zp).^2+omn.*abs(zn).^2,kappa.^2);
ombar=powermean(om,kappa,1);

sig=sqrt(frac((upp.^2+(omp-vrep(ombar,size(om,1),1)).^2).*abs(zp).^2+...
    (upn.^2+(omn-vrep(ombar,size(om,1),1)).^2).*abs(zn).^2,kappa.^2));
sigbar=sqrt(powermean(sig.^2,kappa,1));

[x,y]=vectmult(sqrt(2)*tmat',zp,zn);
[kappa,lambda,theta,phi]=ellparams(x,y);
[ba,bd,bp]=ellband(dt,kappa,lambda,theta,phi);
bom=om-vrep(ombar,size(om,1),1);


sig2=sqrt(ba.^2+bd.^2+bp.^2+bom.^2);
sigbar2=sqrt(powermean(sig2.^2,kappa,1));

%plot(sig,'b'),hold on,plot(sig2,'r')  %Visually identical
reporttest('ELLBAND rotary and elliptical bandwidths match for shifted ellipse, non-unit sample rate',aresame(sigbar,sigbar2,11e-3))

function[]=ellband_test4


load solomon 
use solomon

%Choose central portion where polarization is elliptical
vindex(x,y,z,420:580,1);
vfilt(x,y,z,10,'zeros');
x=anatrans(x,'mirror');
y=anatrans(y,'mirror');
z=anatrans(z,'mirror');


[kappa,lambda,theta,phi,alpha,beta]=ellparams(x,y,z);
[om,up]=instfreq([x,y,z]);
[ombar,upbar]=instfreq([x,y,z],2);
[a,b,c,d,e]=ellband(kappa,lambda,theta,phi,alpha,beta); 
upbar1=sqrt(a.^2+b.^2+d.^2+abs(e).^2);

err=abs(upbar-upbar1).^2./abs(upbar).^2;
err=flipud(sort(err));
err=err(20:end);
%figure, plot([upbar upbar1])
reporttest('ELLBAND trivariate case, sum of bandwith terms for Solomon Islands (removing worst outliers)',allall(err<0.05))


function[]=ellband_figure
t=(0:1:925)';
cxe=zeros(length(t),3);

kappa=3*exp(2*0.393*(t/1000-1));
lambda=0.4+0*t;
phi=(t/1000*5)*2*pi;
theta=pi/4+0*t;

cxe(:,1)=ellsig(kappa,lambda,theta,phi);
om=vdiff(phi,1);  %Since theta is constant

kappa=2.5+0*t;
lambda=zeros(size(t));
for i=2:length(t)
    lambda(i)=real(lambda(i-1)+2*sqrt(1-lambda(i-1).^2)*0.025.*om(i));
end
lambda(1)=nan;
lambda(lambda>1)=1;  

cxe(:,2)=ellsig(kappa,lambda,theta,phi);

[kappa,lambda]=ab2kl(3+zeros(size(t)),2+zeros(size(t)));

%theta=-phi/16.3;  Retrograde precession
theta=phi/14.45;
cxe(:,3)=ellsig(kappa,lambda,theta,phi);

titlestr{1}='Increasing Magnitude';
titlestr{2}='Increasing Eccentricity';
titlestr{3}='Precession';

figure
for i=1:3
    subplot(1,3,i)
    plot(cxe(:,i)),hold on
    plot(cxe(1:200,i)),linestyle 1k 2k
    plot(cxe(2,i),'ko','markersize',10), plot(cxe(end,i),'kx','markersize',10)
    title(titlestr{i})
    xlabel('Displacement East'),axis equal,axis square
    axis([-3.2 3.2 -3.2 3.2])
    ylabel('Displacement North')
    xtick(-3:3)   
    ytick(-3:3)
end
packcols(1,3)
letterlabels(1)

