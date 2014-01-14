function[u,v,w]=hor2uvw(lat,lon,uh,vh)
%HOR2UVW  Converts a horizontal vector on a sphere into a 3D Cartesian vector.
%
%   [U,V,W]=HOR2UVW(LAT,LON,UH,VH) converts the vector with local 
%   horizontal components UH and VH at point (LAT,LON) on a sphere 
%   into a 3D Cartesian vector having components U, V, and W.
%
%   LAT and LON are in degrees.  
%
%   U, V, and W are in a reference frame with the X-axis at zero 
%   degrees longitude and the Z-axis at the North Pole.  
%
%   All input arguments should be arrays of the same size.
%
%   HOR2UVW is inverted by UVW2HOR.
%
%   See JSPHERE for related functions.
%
%   'hor2uvw --t' runs a test.
%
%   Usage: [u,v,w]=hor2uvw(lat,lon,uh,vh);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2008 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(lat, '--t')
    hor2uvw_test,return
end

[phi,theta]=jdeg2rad(lat,lon);

u=-uh.*sin(theta)-vh.*cos(theta).*sin(phi);
v= uh.*cos(theta)-vh.*sin(theta).*sin(phi);
w= vh.*cos(phi);

function[]=hor2uvw_test
 
lat=[0  0  45         0        -90    45]';
lon=[0  90 0         90         0     0]';
uh= [1  0  1          1         1     0]';
vh= [0  1  0          0         0     1]';
u=  [0  0  0         -1         0     -sqrt(2)/2]'; 
v=  [1  0  1          0         1    0]';
w=  [0  1  0          0         0    sqrt(2)/2]';


[u2,v2,w2]=hor2uvw(lat,lon,uh,vh);

tol=1e-6;
reporttest('HOR2UVW example points',aresame(u,u2,tol)&&aresame(v,v2,tol)&&aresame(w,w2,tol))

lon=(0:2:358)-180;
lat=(-90:1:90);
[long,latg]=meshgrid(lon,lat);

uh=randn(size(latg));
vh=randn(size(latg));
[u,v,w]=hor2uvw(latg,long,uh,vh);
[v1,v2,v3]=uvw2sphere(latg,long,u,v,w);

tol=1e-10;
reporttest('HOR2UVW plus UVW2SPHERE for vanishing w',aresame(uh,v2,tol)&&aresame(vh,v3,tol)&&aresame(v1,0*v1,tol))

