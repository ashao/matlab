
function vel = geovel(ga,lat,lon)

% SW_GVEL    Geostrophic velocity
%===================================================================
% GEOVEL   $Revision: 1.3 $  $Date: 1994/10/10 05:04:30 $
%          Copyright (C) CSIRO, Phil Morgan 1992
%
% USAGE:  vel = geovel(ga,lat,lon)
%
% DESCRIPTION:
%    Calculates geostrophic velocity given the geopotential anomaly
%    and position of each station.
% 
% INPUT:
%    ga   = geopotential anomoly relative to the sea surface.
%           dim(mxnstations)
%    lat  = latitude  of each station (+ve = N, -ve = S) [ -90.. +90]
%    lon  = longitude of each station (+ve = E, -ve = W) [-180..+180]
%
% OUTPUT:
%    vel  = geostrophic velocity RELATIVE to the sea surface.
%           dim(m,nstations-1)
%
%    Note: negative vel means flow towards left of section,
%          positive vel means flow towards right of section!!!   SM
%
% AUTHOR:   Phil Morgan   1992/03/26  (morgan@ml.csiro.au)
%           
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.  
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCE: S. Pond & G.Pickard  2nd Edition 1986
%            Introductory Dynamical Oceanogrpahy
%            Pergamon Press Sydney.  ISBN 0-08-028728-X
%            Equation 8.9A p73  Pond & Pickard
%
% NOTE: This calls sw_dist.m.  You can replace the call to this
%       routine if you have a more appropraite distance routine.
% REPLACE sw_dist.m with dist.m 94-10-18 Gregory C. Johnson
% RENAMED dist.m to gj_dist.m 98-01-21 Sabine Mecking
% ADDED CHECK OF INPUT PARAMETERS 98-01-21 Sabine Mecking
%==================================================================

% CALLER:   general purpose
% CALLEE:   gj_dist.m
%
% ----------------------
% CHECK INPUT ARGUMENTS
% ----------------------
if nargin ~=3 
  error('sw_gvel.m: Must pass 3 inout arguments'); 
end

% DETERMINE GA, LAT, LON DIMENSIONS 
[mlat,nlat] = size(lat);
[mlon,nlon] = size(lon);
[m,n] = size(ga);

% MAKE ROW VECTOR OUT OF LAT/LON IF COLUMN VECTOR
if nlat == 1 
  lat = lat';
  [mlat,nlat] = size(lat);
end
if nlon == 1   % column vector
  lon = lon';
  [mlon,nlon] = size(lon);
end

% CHECK DIMENSIONS OF LAT/LON AGAINST GA
if (nlat ~= n) | (nlon~=n)
  error('sw_gvel.m: GA and LAT/LON must have same number of station entries');
end
  
% CONSTANTS
DEG2RAD = pi/180;
RAD2DEG = 180/pi;
OMEGA   = 7.292e-5;  % Angular velocity of Earth  [radians/sec]

%--------------------------------
% CALCULATE GEOSTROPHIC VELOCITY
%--------------------------------
% You may replace the call to sw_dist if you have
% a more appropriate distance routine.
distm = gj_dist(lat,lon);

f     = 2*OMEGA*sin( (lat(1:n-1)+lat(2:n))*DEG2RAD/2 );
lf    = f.*distm;
for im = 1:m
  LF(im,:) = lf;
end %for
vel   = -( ga(:,2:n)-ga(:,1:n-1) ) ./ LF;  

return
%--------------------------------------------------------------------

