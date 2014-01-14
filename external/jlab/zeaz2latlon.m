function[lat,lon]=zeaz2latlon(varargin)
%ZEAZ2LATLON  Compute latitude and longitude viewed by satellite beam.
%
%   [LAT,LON]=ZEAZ2LATLON(ZE,AZ,LATO,LONO,RO) computes the latitude
%   and longitude point [LAT,LON] viewed by a satellite beam with
%   zenith angle ZE and azimuth angle AZ.
%
%   The satellite has nadir point [LATO, LONO] and is located RO 
%   kilometers about the surface of the earth.  
%
%   The radius of the earth is given by RADEARTH.
%  
%   All angles are in degrees.
%
%   ZE and AZ may be arrays of the same size, but LATO, LONO, and 
%   RO must be scalars.
%
%   ZE and AZ are set to be NANs for [LAT,LON] locations that lie
%   outside of the satellite's field of view.
%
%   ZEAZ2LATLON is inverted by LATLON2ZEAZ.
%
%   See also LATLON2ZEAZ, ZE2DIST, ZE2INC.
%
%   'zeaz2latlon --t' runs a test.
%
%   Usage: [lat,lon]=zeaz2latlon(ze,az,lato,lono,ro);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    zeaz2latlon_test,return
end
 
el=varargin{1};
az=varargin{2};
lato=varargin{3};
lono=varargin{4};
ro=varargin{5};

if ~(isscalar(lato) && isscalar(lono) && isscalar(ro))
   error('LATO, LONO, and RO should all be scalars.')
end

R=radearth;

[gamma,rho,tho,phio]=jdeg2rad(el,az,lono,lato);

%r1=(R+ro).*cos(gamma) - sqrt( squared((R+ro).*cos(gamma))-2*R.*ro-ro.^2);

r1=ze2dist(el,ro);

x=R+ro-r1.*cos(gamma);
y=r1.*sin(gamma).*cos(rho);
z=r1.*sin(gamma).*sin(rho);

gammamax=asin(frac(R,R+ro));  %Maximum zenith angle
rmax=sqrt((R+ro).^2-R.^2);    %Maximum distance to contact point
xmin=R+ro-rmax.*cos(gammamax);%Minimum value of x-coordinate  

index=find(x<xmin);
if ~isempty(index)
    lat(index)=nan;
    lon(index)=nan;
end

[x,z]=vectmult(jmat(phio),x,z);
[x,y]=vectmult(jmat(tho),x,y);
[lat,lon]=xyz2latlon(x,y,z);


function[]=zeaz2latlon_test
latlon2zeaz --t 

%reporttest('ZEAZ2LATLON ',aresame())

