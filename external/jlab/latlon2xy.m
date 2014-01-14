function[x,y,d]=latlon2xy(varargin)
% LATLON2XY  Converts latitude and longitude into local Cartesian coordinates.
%
%   [X,Y]=LATLON2XY(LAT,LON,LATO,LONO) converts (LAT,LON) with units
%   of degrees into displacements (X,Y) in a plane tangent to the 
%   earth at the point (LATO, LONO). X and Y have units of kilometers.
%
%   CX=LATLON2XY(LAT,LON,LATO,LONO) with one output argument returns 
%   the location as a complex-valued quantity X+SQRT(-1)*Y. NANs in
%   LAT or LON become NAN+SQRT(-1)*NAN.
%   
%   LAT and LON are arrays of the same size.  LATO and LONO are either
%   also arrays of this size, or else scalars.  X and Y have the same
%   size as the input arrays.  
%
%   X and Y are computed by projecting the tangent plane onto the 
%   sphere using full spherical geometry.  
%
%   The radius of the earth is given by the function RADEARTH.
%
%   Note that X and Y are set to NANs for points on the opposite side
%   of the earth from the tangent plane, that is, where the great 
%   circle distance would exceed RADEARTH * pi/2.  
%   ___________________________________________________________________
%
%   Great circle distance
%
%   [X,Y,D]=LATLON2XY(...) also returns the great circle distance D
%   between the two sets of points.  
%
%   For points on the same side on the earth from the tangent plane,
%   i.e. where the great circle distance is less than RADEARTH * pi/2,
%   LATLON2XY gives the same distance as SPHEREDIST.  
%
%   However, for points on the opposite side of the earth, LATLON2XY
%   returns NANs whereas SPHEREDIST returns the correct distance.
%
%   The great circle distance computed here is useful because it is 
%   a fast computation if X and Y are already known.
%   ___________________________________________________________________
%
%   Cell array input / output
%
%   LATLON2XY returns cell array output given cell array input.  
%
%   That is, if LAT, LON, LATO, and LONO are all cell arrays of length K, 
%   containing K different numerical arrays, then the output will also be 
%   cell arrays of length K.  
%   ___________________________________________________________________
%
%   Complex-valued input / output
%
%   If LON and LAT are complex-valued, the real and imaginary parts are 
%   converted separately to give complex-valued X and complex-valued Y.  
%   ___________________________________________________________________
%
%   LATLON2XY is inverted by XY2LATLON.
%
%   See also XY2LATLON, LATLON2UV, SPHEREDIST.
%
%   Usage:  [x,y]=latlon2xy(lat,lon,lato,lono);
%           cx=latlon2xy(lat,lon,lato,lono);
%
%   'latlon2xy --t' runs some tests.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2011 J.M. Lilly --- type 'help jlab_license' for details        
  
if strcmp(varargin{1}, '--t')
  latlon2xy_test,return
end

na=nargin;
if ischar(varargin{end})
    str=varargin{end};
    na=na-1;
else 
    str='sphere';
end


lat=varargin{1};
lon=varargin{2};
lato=varargin{3};
lono=varargin{4};

R=radearth;

if ~iscell(lat)
    if isreal(lat)
        [x,y,d]=latlon2xy_celloop(nargout,lat,lon,lato,lono,R,str);
    else
        [xr,yr,dr]=latlon2xy_celloop(nargout,real(lat),real(lon),lato,lono,R,str);
        [xi,yi,di]=latlon2xy_celloop(nargout,imag(lat),imag(lon),lato,lono,R,str);
        x=xr+sqrt(-1)*xi;
        y=yr+sqrt(-1)*yi;
        d=dr+sqrt(-1)*di;
    end
else
    for i=1:length(lat)
        if isreal(lat{i})
            [x{i},y{i},d{i}]=latlon2xy_celloop(nargout,lat{i},lon{i},lato{i},lono{i},R,str);
        else
           [xr,yr,dr]=latlon2xy_celloop(nargout,real(lat{i}),real(lon{i}),lato{i},lono{i},R,str);
           [xi,yi,di]=latlon2xy_celloop(nargout,imag(lat{i}),imag(lon{i}),lato{i},lono{i},R,str);
           x{i}=xr+sqrt(-1)*xi;
           y{i}=yr+sqrt(-1)*yi;
           d{i}=dr+sqrt(-1)*di;
        end  
    end
end

function[x,y,d]=latlon2xy_celloop(N,lat,lon,lato,lono,R,str)
d=[];
if numel(lato)==1
    lato=lato+0*lat;
end
if numel(lono)==1
    lono=lono+0*lat;
end

if strcmp(str(1:3),'sma')
     [x,y]=latlon2xy_cartesian(lat,lon-lono,lato,R);
elseif strcmp(str(1:3),'sph')
     if N<3
        [x,y]=latlon2xy_sphere(lat,lon-lono,lato,R);
     else
        [x,y,d]=latlon2xy_sphere(lat,lon-lono,lato,R);
     end
end

if N==1 
  x=real(x)+sqrt(-1)*real(y);
  index=find(isnan(real(x))|isnan(imag(x)));
  if ~isempty(index)
     x(index)=nan+sqrt(-1)*nan;
  end
end

    
function[x,y]=latlon2xy_cartesian(lat,lon,lato,R)

[lat,lon,lato]=jdeg2rad(lat,lon,lato);
r1=R.*cos(lato);
dlon=angle(rot(lon));
x=dlon.*r1;
y=(lat-lato).*R;

function[x,y,d]=latlon2xy_sphere(lat,lon,lato,R)

x=nan*lat;
y=nan*lon;

%  Check to see if I'm in the same hemisphere, or on the other side
%  Here I'm implementing latlon2xyz to find dot product of two vectors
%  and, I avoid recomputing angles

%if ~isreal(lato),lato,end
coslat=cosd(lat);
coslon=cosd(lon);
coslato=cosd(lato);
sinlato=sind(lato);
sinlat=sind(lat);

%  The next line is the dot product, i.e. x.*x2+z.*z2 with
%  x=cos(lat).*cos(lon);z=sin(lat);
%  x2=cos(lato);z2=sin(lato);

index=find(coslat.*coslon.*coslato+sinlato.*sinlat>0);

%  Now implement these equations
%    x=R*cos(lat).*sin(lon);
%    y=-R*cos(lat).*sin(lato).*cos(lon)+R.*cos(lato).*sin(lat);
%  but only for points in the same hemisphere.

%  To derive these equations, calculate the xyz position in space of both 
%  a lat/lon point, and a point on the tangent plane.  It's helpful to
%  define d(x,y)=perpendicular distance to Earth from tangent plane 
%  where  d(x,y)=R - sqrt(R.^2-x^2-y^2)

%  Also note, correctly, that y>0 at the north pole for delta lon=180 
%                         but y>0 at the south pole for delta lon=0;

if ~isempty(index)
    x(index)=R*coslat(index).*sind(lon(index));
    y(index)=-R*coslat(index).*sinlato(index).*coslon(index)+R.*coslato(index).*sinlat(index);
end

%d=2*R.*asin(frac(1,sqrt(2)).*sqrt(1-sqrt(1-frac(x.^2+y.^2,R.^2)))); 
if nargout==3
    if ~isempty(index)
        d=nan*lon;
        d(index)=2*R.*asin(frac(1,sqrt(2)).*sqrt(1-sqrt(1-frac(x(index).^2+y(index).^2,R.^2)))); 
    end
end





function[]=latlon2xy_test
latlon2xy_sphere_test1
latlon2xy_sphere_test2
latlon2xy_sphere_test3
latlon2xy_sphere_test4
latlon2xy_complex_test

function[]=latlon2xy_complex_test

load npg2006
use npg2006
lato=vmean(lat,1);
lono=vmean(lon,1);

lat=anatrans(lat);
lon=anatrans(lon);

lat=lat+lato*sqrt(-1);
lon=lon+lono*sqrt(-1);

[x,y]=latlon2xy(lat,lon,lato,lono);
[lat2,lon2]=xy2latlon(x,y,lato,lono);

tol=1e-8;
b=aresame(lat,lat2,tol) && aresame(lon,lon2,tol);
reporttest('XY2LATLON / LATLON2XY inversion for complex-valued signals',b);



function[]=latlon2xy_sphere_test1


N=100;
tol=1e-3;

lon=2*pi*rand(N,1)-pi;
lat=pi*rand(N,1)-pi/2;
[lat,lon]=jrad2deg(lat,lon);

lat=lat/1000;
lon=lon/1000;

[x,y]=latlon2xy(lat,lon,0,0,'sphere');
[x2,y2]=latlon2xy(lat,lon,0,0,'small');

b=aresame(x,x2,tol) && aresame(y,y2,tol);
reporttest('LATLON2XY Cartesian and spherical algorithms match for small LAT and LON about zero',b);

function[]=latlon2xy_sphere_test2

N=100;
tol1=1e-1;
tol2=1e-1;

lon=2*pi*rand(N,1)-pi;
lat=pi*rand(N,1)-pi/2;
[lat,lon]=jrad2deg(lat,lon);

lat=lat/1000;
lon=lon/1000;

lono=2*pi*rand(N,1)-pi;
lato=pi*rand(N,1)-pi/2;
[lato,lono]=jrad2deg(lato,lono);

lat=lat+lato;
lon=lon+lono;

clear x y lat2 lon2

for i=1:length(lato)
    [x(i,1),y(i,1)]=latlon2xy(lat(i),lon(i),lato(i),lono(i),'sphere');
    [lat2(i,1),lon2(i,1)]=xy2latlon(x(i),y(i),lato(i),lono(i),'sphere');
end

[x2,y2]=latlon2xy(lat,lon,lato,lono);

b=aresame(x,x2,tol1) && aresame(y,y2,tol1);
reporttest('LATLON2XY Cartesian and spherical algorithms match for small LAT and LON perturbations',b);

b=aresame(lat,lat2,tol2) && aresame(lon,lon2,tol2);
reporttest('XY2LATLON Cartesian and spherical algorithms match for small LAT and LON perturbations',b);



function[]=xy2latlon_test

latc=44;
lonc=0;
N=100;

x=randn(N,1)*5;
y=randn(N,1)*5;
[lat,lon]=xy2latlon(x,y,latc,lonc,'small');
[x2,y2]=latlon2xy(lat,lon,latc,lonc,'small');

tol=1e-10;
bool=aresame(x2,x,tol).*aresame(y2,y,tol);
reporttest('XY2LATLON / LATLON2XY conversion', bool)

latc=-44;
lonc=180;
N=100;

x=randn(N,1)*5;
y=randn(N,1)*5;
[lat,lon]=xy2latlon(x,y,latc,lonc,'small');
[x2,y2]=latlon2xy(lat,lon,latc,lonc,'small');

tol=1e-10;
bool(2)=aresame(x2,x,tol).*aresame(y2,y,tol);
reporttest('XY2LATLON / LATLON2XY conversion at 180', bool)

N=100;
tol=1e-6;

R=radearth;
x=frac(R,sqrt(2)).*(2.*rand(N,1)-1);
y=frac(R,sqrt(2)).*(2.*rand(N,1)-1);

lato=2*pi*rand(N,1)-pi;
lono=pi*rand(N,1)-pi/2;
[lato,lono]=jrad2deg(lato,lono);

clear lat lon x2 y2

for i=1:length(lato)
    [lat(i,1),lon(i,1)]=xy2latlon(x(i),y(i),lato(i),lono(i),'sphere');
    [x2(i,1),y2(i,1)]=latlon2xy(lat(i),lon(i),lato(i),lono(i),'sphere');
    [x3(i,1),y3(i,1)]=latlon2xy(lat(i),lon(i),lato(i),lono(i));
end


b=aresame(x,x2,tol) && aresame(y,y2,tol);
reporttest('LATLON2XY spherical algorithm inverts XY2LATLON spherical algorithm',b);

function[]=latlon2xy_sphere_test3
lon=(-180:5:175);
lat=(-90:5:90);
[latg,long]=meshgrid(lat,lon);

N=100;
lat1=180*rand(1,1,N)-90;    
lon1=360*rand(1,1,N);   

latg=vrep(latg,N,3);
long=vrep(long,N,3);

lat1=vrep(vrep(lat1,size(latg,1),1),size(latg,2),2);
lon1=vrep(vrep(lon1,size(latg,1),1),size(latg,2),2);

[x,y,d]=latlon2xy(lat1,lon1,latg,long);
d2=spheredist(lat1,lon1,latg,long);

b1=1;
b2=1;
index=find(d2<radearth*(pi/2));
if ~isempty(index)
    b1=aresame(d2(index),d(index),1e-4);
end
index=find(d2>radearth*(pi/2));
if ~isempty(index)
    b2=allall(isnan(d(index)));
end
        
reporttest('LATLON2XY matches SPHEREDIST for points in same hemisphere',b1);
reporttest('LATLON2XY gives NANs for points in other hemisphere',b2);


function[]=latlon2xy_sphere_test4


N=100;
tol=1e-3;

lon=2*pi*rand(N,1)-pi;
lat=pi*rand(N,1)-pi/2;
[lat,lon]=jrad2deg(lat,lon);

lat=lat/1000;
lon=lon/1000;

[x,y]=latlon2xy(lat,lon,0,0,'sphere');

latc{1}=lat;lonc{1}=lon;lonoc{1}=0;latoc{1}=0;
latc{2}=lat;lonc{2}=lon;lonoc{2}=0;latoc{2}=0;

[xc,yc]=latlon2xy(latc,lonc,latoc,lonoc,'sphere');

b=aresame(xc{1},x,tol) && aresame(yc{1},y,tol)&&aresame(xc{2},x,tol) && aresame(yc{2},y,tol);
reporttest('LATLON2XY cell input / output',b);

