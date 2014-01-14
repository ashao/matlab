function[h]=latratio(lat,h)
%LATRATIO  Set plot aspect ratio for latitude / longitude plot.
%
%   LATRATIO(LAT) sets the aspect ratio of the current axis  
%   correctly a Cartesian plot centered about the latitude LAT,
%   i.e. the x/y aspect ratio is set to COS(LAT).
%
%   Equal distances along the x- and y-axes then correspond to
%   the same physical distance.
%
%   LAT is measured in degrees.
%
%   LATRATIO(LAT,H) does the same for the axis with handle H.
%
%   LATRATIO with no input arguments sets the aspect ratio of 
%   the current axis using the midpoint of the y-axis limits.
%
%   Usage: latratio(lat,h);
%          latratio(lat);
%          latratio;
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007 J.M. Lilly --- type 'help jlab_license' for details
 
if nargin~=0    
    if strcmp(lat, '--t')
        latratio_test,return
    end
end

if nargin==0
   ax=axis;
   lat=frac(1,2)*(ax(3)+ax(4));
end

if nargin<2
    h=gca;
end

set(h,'dataaspectratio',[1./cos(jdeg2rad(lat)) 1 1])
if nargout==0
    clear h
end
function[]=latratio_test
 
%reporttest('LATRATIO',aresame())
