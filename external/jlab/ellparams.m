function[k,l,theta,phi,alpha,beta]=ellparams(varargin)
%ELLPARAMS  Ellipse parameters of a modulated bivariate or trivariate oscillation.
%
%   [KAPPA,LAMBDA,THETA,PHI]=ELLPARAMS(X,Y) where X and Y are analytic 
%   signals, returns the parameters of the complex-valued signal 
%   Z=REAL(X)+i REAL(Y), expressed as a modulated ellipse.
%
%   Here KAPPA is the RMS ellipse amplitude, LAMBDA is the linearity, 
%   THETA is the orientation, and PHI is the instantaneous orbital phase.
%
%   See Lilly and Gascard (2006) and Lilly and Olhede (2010a) for details.
%
%   ELLPARAMS is inverted by ELLSIG, which returns the X and Y signals 
%   given the ellipse parameters.
%   __________________________________________________________________
%
%   Trivariate signals
%
%   ELLPARAMS also works for trivariate signals, which can be expressed as
%   a modulated ellipse in three dimensions.
%
%   [KAPPA,LAMBDA,THETA,PHI,ALPHA,BETA]=ELLPARAMS(X,Y,Z), where X, Y, and Z
%   are all analytic signals, also returns the zenith angle ALPHA and the 
%   azimuth angle BETA in addition to the other ellipse parameters.
%
%   See Lilly (2010) for details on the trivariate case.
%   __________________________________________________________________
%
%   Alternate form
%
%   [...]=ELLPARAMS(M) also works, where M is a two-column or three-column
%   matrix, i.e. M=[X Y] or M=[X Y Z], for bivariate or trivariate signals,
%   respectively. 
%   __________________________________________________________________
%
%   'ellparams --t' runs a test.
%
%   See also ELLSIG, ELLDIFF, ELLVEL, ELLRAD, KL2AB, AB2KL. 
%
%   Usage: [kappa,lambda,theta,phi]=ellparams(x,y);
%   [kappa,lambda,theta,phi,alpha,beta]=ellparams(x,y,z);
%   [kappa,lambda,theta,phi]=ellparams([x,y]);
%   [kappa,lambda,theta,phi,alpha,beta]=ellparams([x,y,z]);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009--2011 J.M. Lilly --- type 'help jlab_license' for details

if strcmp(varargin{1}, '--t')
    ellparams_test,return
end

[na,k,l,theta,phi,alpha,beta]=vempty;

if nargin==1
    if ~isempty(varargin{1})
        if size(varargin{1},2)==2
            na=2;
            x=varargin{1}(:,1);
            y=varargin{1}(:,2);
        elseif size(varargin{1},2)==3
            na=3;
            x=varargin{1}(:,1);
            y=varargin{1}(:,2);
            z=varargin{1}(:,3);
        else
            error('For ELLPARAMS(X) with one input argument, X must have either two or three columns.')
        end
    end
elseif nargin==2
    na=2;
    x=varargin{1};
    y=varargin{2};
elseif nargin==3
    na=3;
    x=varargin{1};
    y=varargin{2};
    z=varargin{3};
else
    error('ELLPARAMS must have one, two, or three input arguments.')
end    
            
if na==2
    [k,l,theta,phi]=ellconv_xy2kl(abs(x),abs(y),angle(x),angle(y));
    varargout{1}=k;
    varargout{2}=l;
    varargout{3}=theta;
    varargout{4}=phi;
elseif na==3
    [nx,ny,nz]=normvect(x,y,z);
    warning('off','MATLAB:log:logOfZero')
    alpha=imag(log(sqrt(-1)*nx-ny));
    warning('on','MATLAB:log:logOfZero')
    beta=imag(log(nz+sqrt(-1)*sqrt(nx.^2+ny.^2)));
    [x,y,z]=vectmult(jmat3(-alpha,3),x,y,z);
    [x,y,z]=vectmult(jmat3(-beta,1),x,y,z);
    [k,l,theta,phi]=ellconv_xy2kl(abs(x),abs(y),angle(x),angle(y));
end
      
function[kappa,lambda,theta,phi]=ellconv_xy2kl(X,Y,phix,phiy)
%phia=double(phix+phiy+pi/2)/2;
%phid=double(phix-phiy-pi/2)/2;
%P=double(frac(1,2)*sqrt(squared(X)+squared(Y)+2.*X.*Y.*cos(2*phid)));
%N=double(frac(1,2)*sqrt(squared(X)+squared(Y)-2.*X.*Y.*cos(2*phid)));

phia=(phix+phiy+pi/2)/2;
phid=(phix-phiy-pi/2)/2;

P=frac(1,2)*sqrt(squared(X)+squared(Y)+2.*X.*Y.*cos(2*phid));
N=frac(1,2)*sqrt(squared(X)+squared(Y)-2.*X.*Y.*cos(2*phid));

phip=unwrap(phia+imlog(X.*rot(phid)+Y.*rot(-phid)));
phin=unwrap(phia+imlog(X.*rot(phid)-Y.*rot(-phid)));

kappa=sqrt(P.^2+N.^2);
lambda=frac(2*P.*N.*sign(P-N),P.^2+N.^2);

%For vanishing linearity, put in very small number to have sign information 
lambda(lambda==0)=sign(P(lambda==0)-N(lambda==0))*(1e-10);

theta=phip/2-phin/2;
phi=  phip/2+phin/2;

theta=unwrap(theta);
phi=unwrap(phi);


function[]=ellparams_test
 
t=(0:1:925)';
kappa=3*exp(2*0.393*(t/1000-1));
lambda=0.5+0*kappa;
phi=(t/1000*5)*2*pi;
theta=pi/4+phi./14.45;

beta=pi/6+phi./14.45*lambda(1);
alpha=pi/6-phi./14.45*2*lambda(1)*sqrt(2);

[x,y,z]=ellsig(kappa,lambda,theta,phi,alpha,beta);
[kappa2,lambda2,theta2,phi2,alpha2,beta2]=ellparams(x,y,z);

x1=[kappa lambda theta phi alpha beta];
x2=[kappa2 lambda2 theta2 phi2 alpha2 beta2];
reporttest('ELLPARAMS rapidly changing trivariate ellipse',aresame(x1,x2,1e-8))


% 
% %/*******************************************************************
% %Flip ellipse parameters if unit normal is pointing radially inwards
% %Components of unit normal vector to surface of earth
% [nx,ny,nz]=normvect(xr,yr,zr);
% 
% %Replicate x, y, and z along columns
% xmat=vrep(x,size(nx,2),2)./radearth;
% ymat=vrep(y,size(ny,2),2)./radearth;
% zmat=vrep(z,size(nz,2),2)./radearth;
% 
% %Projection of normal vector to plane onto normal to sphere
% proj=xmat.*nx+ymat.*ny+zmat.*nz;
% 
% bool=(proj<0);
% theta(bool)=-theta(bool);
% lambda(bool)=-lambda(bool);
% beta(bool)=pi+beta(bool);
% nx(bool)=-nx(bool);
% ny(bool)=-ny(bool);
% nz(bool)=-nz(bool);
% 
% if length(find(bool))>0
%     [xr2,yr2,zr2]=ellsig(kappa,lambda,theta,phi,alpha,beta);
%     tol=1e-6;
%     reporttest('ELLIPSEXTRACT adjustment for sign of normal vector',aresame(xr2,xr,tol)&&aresame(yr2,yr,tol)&&aresame(zr2,zr,tol))
% end
% 
% %dev=1-abs(proj);
% %figure,plot(dev)
% [latn,lonn]=xyz2latlon(nx*radearth,ny*radearth,nz*radearth);
% %\*******************************************************************
