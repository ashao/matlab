function[d]=ze2dist(el,ro,R)
%ZE2DIST Converts beam zenith angle into distance to surface.
%
%   D=ZE2DIST(ZE,RO) where RO is a satellite elevation above the surface
%   of the earth, and ZE is the beam zenith angle in degrees, returns 
%   the distance D along the beam to the surface of the Earth.  
% 
%   The zenith angle is defined as the angle between the beam and the
%   nadir line.  The Earth is approximated as a sphere of radius RADEARTH.
%
%   D=ZE2DIST(ZE,RO,R) optionally uses a sphere of radius R, in kilometers.
%
%   See also ZE2INC, ZEAZ2LATLON, LATLON2ZEAZ.
%   
%   Usage: d=ze2dist(ze,ro);
%          d=ze2dist(ze,ro,R);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006 J.M. Lilly --- type 'help jlab_license' for details    

if strcmp(el, '--t'),return,end

if nargin==2
   R=radearth;
end

if (~aresame(size(el),size(ro))) && (~isscalar(el) && ~isscalar(ro))
   error('ZE and RO must be the same size, or one must be a scalar.')
end
if ~isscalar(R)
   error('R must be a scalar.')
end

c=2*pi/360;
el=el*c;

[rp,d]=quadform(1,-2*cos(el).*(ro+R),2*R*ro+ro.^2);
%The larger root is the distance to the far side of the earth
index=find(d<0);
if ~isempty(index)
    d(index)=nan;
end


function[]=ze2dist_test

%reporttest('ze2dist',bool);

