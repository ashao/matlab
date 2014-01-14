%
% JSPHERE  Spherical geometry and derivatives.  
%
%   See also: 
%   CONTENTS, JARRAY, JMATH, JPOLY, JGRAPH, JSTRINGS, JSTATS, JSIGNAL,
%   JELLIPSE, JCELL, VTOOLS, JOCEANS, *JSPHERE*, JSATFUN, JTRIADS, JPAPERS
%   ______________________________________________________________________
%
% Distances
%   spheredist    - Computes great circle distances on a sphere.
%   spheresort    - Sorted great circle distances to nearby points on a sphere.
%   inregion      - Tests whether lat/lon points lie within a specified box.
%   lonshift      - Shifts longitude origin for plotting purposes.
%   [-- see also JOCEANS]
%
% Spherical geometry
%   xyz2latlon - Converts 3D Cartesian coordinates into latitude and longitude.
%   latlon2xyz - Converts latitude and longitude into 3D Cartesian coordinates.
%   uvw2sphere - Converts a 3D Cartesian vector to a 3D spherical vector.
%   sphere2uvw - Converts a 3D spherical vector to a 3D Cartesian vector.
%   uvw2hor    - Projects a 3D Cartesian vector into a horizontal vector on a sphere.
%   hor2uvw    - Converts a horizontal vector on a sphere into a 3D Cartesian vector.
%
% Div, curl, and grad
%   spherediv  - Divergence of a vector field on the surface of a sphere.
%   spheregrad - Gradient of a field on the surface of a sphere.
%   spherecurl - Curl of a vector field on the surface of a sphere.
%   spherelap  - Laplacian of a field on the surface of a sphere.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2011 J.M. Lilly --- type 'help jlab_license' for details   

help Jsphere

if 0          
           %Distances
             spheredist    %-      Computes great circle distances on a sphere.
             spheresort    %-      Sorted great circle distances to nearby points on a sphere.
             inregion      %-      Tests whether lat/lon points lie within a specified box.
             lonshift      %-      Shifts longitude origin for plotting purposes.
          
           %Spherical geometry
             xyz2latlon %-      Converts 3D Cartesian coordinates into latitude and longitude.
             latlon2xyz %-      Converts latitude and longitude into 3D Cartesian coordinates.
             uvw2sphere %-      Converts a 3D Cartesian vector to a 3D spherical vector.
             sphere2uvw %-      Converts a 3D spherical vector to a 3D Cartesian vector.
             uvw2hor    %-      Projects a 3D Cartesian vector into a horizontal vector on a sphere.
             hor2uvw    %-      Converts a horizontal vector on a sphere into a 3D Cartesian vector.
          
           %Div, curl, and grad
             spherediv  %-      Divergence of a vector field on the surface of a sphere.
             spheregrad %-      Gradient of a field on the surface of a sphere.
             spherecurl %-      Curl of a vector field on the surface of a sphere.
             spherelap  %-      Laplacian of a field on the surface of a sphere.
end

