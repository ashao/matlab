function[bool]=inregion(region,lat,lon)
% INREGION  Tests whether lat/lon points lie within a specified box.
%
%   BOOL = INREGION(REGION,LAT,LON), where LAT and LON are arrays having
%   the same size, returns an array BOOL which is true (=1) for all
%   LAT/LON pairs which lie inside REGION, and false (=0) otherwise.
%
%   REGION has the format [WEST EAST SOUTH NORTH];
%
%   All input arrays are in degrees, but all longitudes may either be 
%   specified on the interval [-180, 180] or on the interval [0, 360].
%
%   The region may overlap the prime meridian (LON=0) or the dateline
%   (LON=180).  Region boundaries are interpreted to exclude the poles. 
%
%   Usage:  bool=inregion(region,lat,lon);
%
%   'inregion --t' runs a test 
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006--2009 J.M. Lilly --- type 'help jlab_license' for details        

if strcmp(region,'--t')
     inregion_test,return
end


west=region(1);
east=region(2);
south=region(3);
north=region(4);

west=deg180(west);
east=deg180(east);
lon=deg180(lon);

boollat=(lat>=south)  & (lat<=north);
if west<east
    %Dateline is not inside region
    boollon=(lon>=west)  & (lon<=east);
elseif east<west 
    boollon=(lon>=west) | (lon<=east);
end

bool=logical(boollat.*boollon);

function[]=inregion_test
lat=[58 76];
lon=[-52.5 -52.5];
region=[-63 -41 52 65];
ans1=inregion(region,lat,lon);
reporttest('INREGION Labrador Sea', aresame(ans1,[1 0]))

lat=[58 76];
lon=360+[-52.5 -52.5];
region=[-63 -41 52 65];
ans1=inregion(region,lat,lon);
reporttest('INREGION Labrador Sea, differing longitude conventions', aresame(ans1,[1 0]))

lat=[58 58];
lon=[170 150];
region=[160 -165 52 65];
ans1=inregion(region,lat,lon);
reporttest('INREGION enclosing dateline', aresame(ans1,[1 0]))

lat=[58 58];
lon=[170 150];
region=[160 195 52 65];
ans1=inregion(region,lat,lon);
reporttest('INREGION enclosing dateline', aresame(ans1,[1 0]))





