function [mat,lat,lon,bool]=read_sand(dirname,region)
%READ_SAND   Read topography data from the Smith and Sandwell Database. 
%
%   [TOPO,LAT,LON]=READ_SAND(DIRNAME,REGION) returns the global topography 
%   from Smith and Sandwell database in directory DIRNAME within REGION.
%
%   REGION is an array [WEST EAST SOUTH NORTH] with units of degrees.
%
%   TOPO is a matrix of topography at longitudes LON and latitudes LAT. 
%   LON is a row vector of longitudes and LAT a column vector of latitudes.
% 
%   TOPO is in units of kilometers and is positive for above sea level, 
%   and negative for below sea level.   
%
%   Note that the Smith and Sandwell database is defined for latitudes 
%   between -80.738 and 80.738.
%   __________________________________________________________________
%
%   Data and documentation
%
%   To use this function, you will need the file "topo_14.1.img" which is 
%   available via anonymouse ftp from 
%
%       ftp://topex.ucsd.edu/pub/global_topo_2min/topo_14.1.img 
%
%   and may be found together with other topography/bathymetry products
%   and documentation at
%
%       http://topex.ucsd.edu/WWW_html/mar_topo.html.
%
%   The reference for the Smith and Sandwell Database is
%
%      Smith, W. H. F., and D. T. Sandwell, Global seafloor topography 
%      from satellite altimetry and ship TOPO soundings, Science, 
%      v. 277, p. 1957-1962, 26 Sept., 1997.
%   __________________________________________________________________
%
%   License
%
%   READ_SAND has the same functionality as "extract_1m" by Catherine de 
%   Groot-Hedlin, but has been completely rewritten. Unlike the latter, 
%   you are free to use and re-distribute READ_SAND, under the terms of the 
%   JLAB license in "jlab_license.m".
%
%   Usage: [topo,lat,lon]=read_sand(dirname,[west east south north]);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006--2011 J.M. Lilly --- type 'help jlab_license' for details


%   __________________________________________________________________
%
%   Ship sounding flag
%
%   [TOPO,LAT,LON,BOOL]=READ_SAND(DIRNAME,REGION) also returns a boolean
%   matrix of the same size as TOPO, which has is true when the
%   corresponding element of TOPO is from a ship sounding and false when 
%   the corresponding element is an interpolated value.


%This is self-contained

if strcmp(dirname, '--t')
    read_sand_test,return
end


%dirname=jlab_settings('dirnames.sandwell');
%   See also READ_IBCAO, L_IBCAO.

minlat=-80.738;
maxlat=80.738;

west=deg180(region(1));
east=deg180(region(2));
south=max(region(3),minlat);
north=min(region(4),maxlat);

if east==0
    east=-1e-10;
end

ddeg=1/60;  
nbytes=21600*2;  % 2-byte integers
skip=17280;

if north<-80.738||south>80.738
  error('The Smith and Sandwell database is defined within +/- 80.738 degrees latitude.')
end

fid=fopen([dirname '/topo_14.1.img'], 'r','b');  
if (fid < 0)
     error(['Could not open database topo_14.1.img.']);
end

north=north(:);
south=south(:);
east=east(:);
west=west(:);

if (west<0)&&(east>=0)
    [mat1,lat,lon1]=read_sand1(fid,south,north,west,360-ddeg,ddeg,minlat,skip,nbytes);
    [mat2,lat,lon2]=read_sand1(fid,south,north,0,east,ddeg,minlat,skip,nbytes);
    mat=[mat1 mat2];
    lon=[lon1 lon2];
else
    [mat,lat,lon]=read_sand1(fid,south,north,west,east,ddeg,minlat,skip,nbytes);
end

fclose(fid);

if nargout==4
    bool=isodd(topo);
end
function[mat,lat,lon]=read_sand1(fid,south,north,west,east,ddeg,minlat,skip,nbytes);
[west,east]=deg360(west,east);

arg=log(tand(45+minlat./2)./tand(45+south./2));
isouth=floor(arg/jdeg2rad(ddeg)) + skip+1;

arg=log(tand(45+minlat./2)./tand(45+north./2));
inorth=floor(arg/jdeg2rad(ddeg)) + skip+1;

iwest=ceil(west./ddeg);
ieast= floor(east./ddeg)-1;

ilat=inorth:isouth;
ilon=iwest:ieast;

Nlats=length(ilat);
Nlons=length(ilon);

mat=zeros(Nlats,Nlons);

for i=1:Nlats
        offset=ilat(i)*nbytes+iwest*2;
        fseek(fid, offset, 'bof');
        mat(isouth-ilat(i)+1,:)=fread(fid,[1,Nlons],'int16');
end

lat=zeros(Nlats,1);
argo=tand(45+minlat/2);

%Can't use jdeg2rad here since it should be between 0 and 2pi
arg=exp(frac(2*pi,360).*(ddeg*(skip-(ilat-0.5)))).*argo;
lat(isouth-ilat+1)=2*atand(arg)-90;
lon=deg180(ddeg*(ilon+0.5));

function[]=read_sand_test
if exist('extract_1m')==2  
    dirname=jlab_settings('dirnames.sandwell');
    %dirname='/Users/lilly/Data/topography/sandwell';
    region=[-70 -35 50 70];
    [topo1,lat1,lon1]=read_sand(dirname,region);
    [topo2,lat2,lon2]=extract_1m([region(3:4) region(1:2)]);

    [long,latg]=meshgrid(lon1,lat1);
    bool=inregion(region,latg,long);
    reporttest('READSAND strictly within region',allall(bool==1))

    [long,latg]=meshgrid(lon2,lat2);
    bool=inregion(region,latg,long);
    topo2=topo2(bool);
    topo2=reshape(topo2,2492,2100);

    reporttest('READSAND matches EXTRACT_1M',aresame(topo1,topo2))
end




