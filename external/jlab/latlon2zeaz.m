function[el,az]=latlon2zeaz(varargin)
%LATLON2ZEAZ  Compute zenith and azimuth angles for satellite beam.
%
%   [ZE,AZ]=LATLON2ZEAZ(LAT,LON,LATO,LONO,RO) computes the zenith
%   angle ZE and azimuth angle AZ from a satellite to points on the 
%   surface of the earth, specified by [LAT,LON].
%
%   The satellite has nadir point [LATO, LONO] and is located RO 
%   kilometers about the surface of the earth.  
%
%   The radius of the earth is given by RADEARTH.
%
%   All angles are in degrees.
%
%   LAT and LON may be arrays of the same size, but LATO, LONO, and 
%   RO must be scalars.
%
%   ZE and AZ are set to be NANs for [LAT,LON] locations that lie
%   outside of the satellite's field of view.
%
%   LATLON2ZEAZ is inverted by ZEAZ2LATLON.
%
%   See also ZEAZ2LATLON, ZE2DIST, ZE2INC.
%
%   'latlon2zeaz --t' runs a test.
%
%   Usage: [ze,az]=latlon2zeaz(lat,lon,lato,lono,ro);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006--2008 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    latlon2zeaz_test,return
end
 
 
lat=varargin{1};
lon=varargin{2};
lato=varargin{3};
lono=varargin{4};
ro=varargin{5};

if ~(isscalar(lato) && isscalar(lono) && isscalar(ro))
   error('LATO, LONO, and RO should all be scalars.')
end

R=radearth;
[th,phi,tho,phio]=jdeg2rad(lon,lat,lono,lato);

[x,y,z]=latlon2xyz(lat,lon);
[x,y]=vectmult(jmat(-tho),x,y);
[x,z]=vectmult(jmat(-phio),x,z);

gamma=atan(frac(sqrt(y.^2+z.^2),R+ro-x));
rho=imlog(y+sqrt(-1)*z);

[az,el]=jrad2deg(rho,gamma);

gammamax=asin(frac(R,R+ro));  %Maximum zenith angle
rmax=sqrt((R+ro).^2-R.^2);    %Maximum distance to contact point
xmin=R+ro-rmax.*cos(gammamax);%Minimum value of x-coordinate  

index=find(x<xmin);
if ~isempty(index)
    az(index)=nan;
    el(index)=nan;
end


function[]=latlon2zeaz_test
xy2zeaz_test; 


function[]=xy2zeaz_test
N=100;
tol=1e-4;
x=rand(N,1)*1000-500;
y=rand(N,1)*1000-500;

ro=657;
lato=0;
lono=0;

[lat,lon]=xy2latlon(x,y,lato,lono);
[el,az]=latlon2zeaz(lat,lon,lato,lono,ro);
[lat2,lon2]=zeaz2latlon(el,az,lato,lono,ro);
[x2,y2]=latlon2xy(lat2,lon2,lato,lono);

b=aresame(x,x2,tol) && aresame(y,y2,tol);
reporttest('LATLON2ZEAZ inverts ZEAZ2LATLON for scalar [LATO,LONO]=[0,0]',b);


N=1000;
tol=1e-4;
x=(rand(N,1)-1/2)*10000;
y=(rand(N,1)-1/2)*10000;

ro=657;
lato=70;
lono=110;

[lat,lon]=xy2latlon(x,y,lato,lono);
[el,az]=latlon2zeaz(lat,lon,lato,lono,ro);
[lat2,lon2]=zeaz2latlon(el,az,lato,lono,ro);
[x2,y2]=latlon2xy(lat2,lon2,lato,lono);

vindex(x,y,x2,y2,find(~isnan(x2)),1);
b=aresame(x,x2,tol) && aresame(y,y2,tol);
reporttest('LATLON2ZEAZ inverts ZEAZ2LATLON for scalar [LATO,LONO]=[70,110]',b);


N=100;
tol=1e-4;
ro=657;

lono=deg180(rand(N,1)*360);
lato=rand(N,1)*180-90;

lon=deg180(rand(N,1)*360);
lat=rand(N,1)*180-90;

clear lat2 lon2 el az
for i=1:N
    [el(i),az(i)]=latlon2zeaz(lat(i),lon(i),lato(i),lono(i),ro);
    [lat2(i,1),lon2(i,1)]=zeaz2latlon(el(i),az(i),lato(i),lono(i),ro);
end

vindex(lat,lon,lat2,lon2,find(~isnan(lat2)),1);
b=aresame(lat,lat2,tol) && aresame(lon,lon2,tol);
reporttest('LATLON2ZEAZ inverts ZEAZ2LATLON for random [LATO,LONO] and [LAT,LON]',b);



