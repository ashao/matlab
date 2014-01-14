function[nx,ny,nz]=normvect(x,y,z)
%NORMVECT  Unit normal vector to the ellipse plane in three dimensions.
%
%   [NX,NY,NZ]=NORMVECT(X,Y,Z) returns the three components of the unit 
%   normal vector to the plane containing the ellipse specified by the 
%   three complex-valued arrays X, Y, and Z.
%
%   In vector notation the normal vector is defined as N=IMAG(X) x REAL(X),
%   where ``x'' is the vector cross product, and the unit normal is N/||N||.
%
%   The input arrays and output arrays are all the same size.
% 
%   See Lilly (2010) for details.
%
%   'normvect --t' runs a test.
%
%   Usage: [nx,ny,nz]=normvect(x,y,z);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2010 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(x, '--t')
    normvect_test,return
end
 

nx= (imag(y).*real(z)-imag(z).*real(y));
ny=-(imag(x).*real(z)-imag(z).*real(x));
nz= (imag(x).*real(y)-imag(y).*real(x));
  
denom=sqrt(nx.^2+ny.^2+nz.^2);
nx=frac(nx,denom);
ny=frac(ny,denom);
nz=frac(nz,denom);

function[]=normvect_test
 
load solomon 
use solomon

x=anatrans(x./1e4);
y=anatrans(y./1e4);
z=anatrans(z./1e4);

[nx,ny,nz]=normvect(x,y,z);

dot=nx.*x+ny.*y+nz.*z;

reporttest('NORMVECT parallel part of X_+ vanishes, Solomon Islands',allall(abs(dot)<1e-10))


%Choose central part where signal is elliptical
vindex(x,y,z,7000:10000,1);

[omx,upx]=instfreq(x);
[omy,upy]=instfreq(y);
[omz,upz]=instfreq(z);
     
dx=x.*(upx+sqrt(-1)*omx);
dy=y.*(upy+sqrt(-1)*omy);
dz=z.*(upz+sqrt(-1)*omz);

[nx,ny,nz]=normvect(x,y,z);

dot=nx.*dx+ny.*dy+nz.*dz; %Parallel part of derivative

[dnx,dny,dnz]=vdiff(nx,ny,nz,1);

dot2=-(dnx.*x+dny.*y+dnz.*z); 

err=abs(dot-dot2).^2./abs(dot).^2;
err=flipud(sort(err));
err=err(60:end);


reporttest('NORMVECT derivative of parallel part matches, Solomon Islands (removing worst outliers)',allall(err<0.05))

