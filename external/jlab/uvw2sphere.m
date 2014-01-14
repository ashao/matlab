function[dr,dth,dphi]=uvw2sphere(lat,lon,u,v,w)
%UVW2SPHERE  Converts a 3D Cartesian vector to a 3D spherical vector.
%
%   [V1,V2,V3]=UVW2SPHERE(LAT,LON,U,V,W) converts a Cartesian 3-vector
%   with components U, V, and W located at point (LAT, LON) into a 
%   vector in spherical coordinates with components V1, V2, and V3.
%
%   LAT and LON are in degrees.
%
%   The Cartesian vector [U,V,W] is in a reference frame with the
%   X-axis at zero degrees longitude and the Z-axis at the North Pole.  
%
%   The vector in the spherical coordinate system has radial component
%   V1, longitudinal component V2, and latitudinal component V3.
%
%   All input arguments should be arrays of the same size.
%
%   UVW2SPHERE is inverted by SPHERE2UVW.
%
%   See JSPHERE for related functions.
%
%   'uvw2sphere --t' runs a test.
%
%   Usage: [v1,v2,v3]=uvw2sphere(lat,lon,u,v,w);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(lat, '--t')
    uvw2sphere_test,return
end

[phi,theta]=jdeg2rad(lat,lon);

dr     =  cos(phi).*cos(theta).*u +  cos(phi).*sin(theta).*v  + sin(phi).*w;
dth    =           -sin(theta).*u +            cos(theta).*v;
dphi   = -sin(phi).*cos(theta).*u -  sin(phi).*sin(theta).*v  + cos(phi).*w;

function[]=uvw2sphere_test
 
lon=(1e-10:2:360)-180;
lat=(-90:2:90);
[lon,lat]=meshgrid(lon,lat);

u=randn(size(lat));
v=randn(size(lat));
w=randn(size(lat));
[dr,dth,dphi]=uvw2sphere(lat,lon,u,v,w);
[u2,v2,w2]=sphere2uvw(lat,lon,dr,dth,dphi);

tol=1e-6;
reporttest('UVW2SPHERE inverts SPHERE2UVW',aresame(u,u2,tol)&&aresame(v,v2,tol)&&aresame(w,w2,tol))
