function[g,w]=aquaplot(lato,lono,hdg)
%AQUAPLOT  Plot Aquarius satellite radiometer footprint.
%
%   AQUAPLOT(LATO,LONO,HDG) plots the footprint of the three Aquarius
%   satellite radiometers in a Mercator projection.
%
%   The satellite nadir point is [LATO,LONO], and HDG specifies the 
%   satellite heading in degrees measured counterclockwise from East.
%   LAT and LON are given in degrees.
%
%   A single contour is drawn around the half-power point.
%
%   The Aquarius color convention is used, that is
%
%           Inner Beam      Green
%           Middle Beam     Red
%           Outer Beam      Blue
%
%   and the nadir point is plotted with a cyan dot.
%
%   See also AQUAPRINT, AQUASAL.
%
%   'aquaplot --f1' makes a sample figure.  
%   'aquaplot --f2' makes a second sample figure, for which you will need 
%        optional matfile 'orbit.mat'. See 'jlab_matfiles' for details.
%
%   Usage:  aquaplot(lat,lon,hdg);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006--2008 J. M. Lilly --- type 'help jlab_license' for details    


if strcmp(lato, '--f1')
  aquaplot_fig1,return
elseif strcmp(lato, '--f2')
  aquaplot_fig2,return
elseif strcmp(lato, '--f')
  aquaplot_fig1,aquaplot_fig2,return
end
    
    

for i=1:length(lono)
    th=jdeg2rad(hdg(i));
    x=(250:4:720);   
    y=(-200:4:200);

    [x,y]=meshgrid(x,y);
    [x,y]=vectmult(jmat(th-pi/2),x,y);

    [lat,lon]=xy2latlon(x,y,lato(i),lono(i));
    [g,w]=aquaprint(lat,lon,lato(i),lono(i),hdg(i));

    ci=.5; 
    sty='grb';
    for j=1:3
        contour(lon,lat,g(:,:,j),ci,'w','linewidth',4),hold on
        contour(lon,lat,g(:,:,j),ci,sty(j),'linewidth',2);hold on
%        contour(lon,lat,g(:,:,j).*w./maxmax(w),ci,[sty(j) '--'],'linewidth',1),hold on
    end
    plot(lono(i),lato(i),'co','markerfacecolor','c','markersize',5)
end

function[]=aquaplot_fig1
figure
aquaplot(0,0,90);

grid on
title('Radiometer footprint with nadir at (0,0).')
orient landscape
fontsize 16 14 14 14
set(gca,'units','inches','position',[1 1  5+7/8 4+5/8])

disp('The first figure shows the approximate Aquarius radiometer beam patterns') 
disp('at the equator.  The axes are in degree.')

function[]=aquaplot_fig2
disp('Be patient, this one takes a while....')

[lon,lat]=vempty;

berror=0;
try 
    load orbit
    if ~isfield(orbit,'creator')
        berror=1;
    elseif ~aresame(orbit.creator,'Created by the script L_ORBIT by J. M. Lilly')
        berror=1;
    end
catch
    berror=1;
end
    
if berror
    error('You must have optional file "orbit.mat".  See "jlab_matfiles" for details.')
end
figure,
use orbit
iia=min(find(lon(:,4)>0));
vshift(num,lat,lon,hdg,iia,1);
iib=min(find(vdiff(lon(:,4),1)>90));
vindex(num,lat,lon,hdg,(1:iib-1)',1);

for i=1:3
    iia=min(find(lon(:,i)>0));
    lat(1:iia-1,i)=nan;
    lon(1:iia-1,i)=nan;
end

plot(lon,lat),hold on
linestyle g r b 2k

ii=1:25:length(lon);
ii=ii(deg180(lon(ii,4))<170);
aquaplot(lat(ii,4),lon(ii,4),hdg(ii));
plot(lon(ii,4),lat(ii,4),'co')
axis([-180 180 -90 90])
title('Aquarius beam patterns through one orbit')

disp('   ')
disp('This figure shows the Aquarius beam patterns through a single orbit,')
disp('correctly represented in a Mercator projection.')

if 0
orient landscape
fontsize 16 14 14 14
set(gca,'units','inches','position',[1 1  5+7/8 4+5/8])

print -dpng orbitpatterns
end
