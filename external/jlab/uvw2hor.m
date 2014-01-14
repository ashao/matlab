function[uh,vh]=uvw2hor(lat,lon,u,v,w)
%UVW2HOR  Projects a 3D Cartesian vector into a horizontal vector on a sphere.
%
%   [UH,VH]=UVW2HOR(LAT,LON,U,V,W) takes the 3D Cartesian vector with
%   components U, V, and W located at point (LAT,LON) and projects it to
%   find the local horizontal components of the vector UH and VH.  
%
%   LAT and LON are in degrees.  
%
%   U, V, and W are in a reference frame with the X-axis at zero 
%   degrees longitude and the Z-axis at the North Pole.  
%
%   All input arguments should be arrays of the same size.
%
%   UVW2HOR inverts HOR2UVW, but the reverse is not true, since some
%   information is lost during the projection.
%
%   See JSPHERE for related functions.
%  
%   'uvw2hor --t' runs a test.
%
%   Usage: [uh,vh]=uvw2hor(lat,lon,u,v,w);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(lat, '--t')
    uvw2hor_test,return
end
 
[phi,theta]=jdeg2rad(lat,lon);

uh=v.*cos(theta)-u.*sin(theta);
vh=w./cos(phi);

function[]=uvw2hor_test

lon=(1e-10:2:360)-180;
lat=(-90:2:90);
[lon,lat]=meshgrid(lon,lat);

uh=randn(size(lat));
vh=randn(size(lat));
[u,v,w]=hor2uvw(lat,lon,uh,vh);
[uh2,vh2]=uvw2hor(lat,lon,u,v,w);


tol=1e-6;
reporttest('UVW2HOR inverts HOR2UVW',aresame(uh,uh2,tol)&&aresame(vh,vh2,tol))
