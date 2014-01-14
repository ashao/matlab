function[g,w]=aquaprint(lat,lon,lato,lono,hdg)
%AQUAPRINT  Compute Aquarius satellite radiometer footprints.
%
%   G=AQUAPRINT(LAT,LON,LATO,LONO) returns the footprints of the three 
%   Aquarius satellite radiometers at points [LAT,LON] for a satellite
%   having a nadir point at [LATO,LONO].
%
%   All angles are given in degrees.
%
%   AQUAPRINT(..., HDG) optionally specifies satellite heading measured 
%   counterclockwise from East.  The default value is HDG=90, i.e. the
%   satellite heading is due North.
%
%   [G,W]=AQUAPRINT(...) also returns the weighting function W needed for
%   integrating the footprint over latitude and longitude coordinates on
%   the surface of the earth.  
%
%   W is a Jacobian, defined as the ratio of the differential area for 
%   the unit sphere around the satellite to the differential area of the 
%   surface of the earth.
%
%   See also AQUAPLOT, AQUASAL.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006--2007 J.M. Lilly --- type 'help jlab_license' for details    


if nargin ==4
    hdg=90;
end
[ro,elo,azo,omo]=aquaparam;
omo=jdeg2rad(omo);
% Leave other angles as degrees

[el,az]=latlon2zeaz(lat,lon,lato,lono,ro);

azo=azo+(hdg-90);


[dom,g]=vzeros(size(el,1),size(el,2),3);
for i=1:3
    dom(:,:,i)=spheredist(90-el,az,90-elo(i),azo(i),1);
    sig=omo(i)/(sqrt(2*log(2)));
    g(:,:,i)=sig.*sqrt(2*pi).*simplepdf(dom(:,:,i),0,sig,'gaussian');
end
 
r1=ze2dist(el,ro);

inc=ze2inc(el,ro);
inc=jdeg2rad(inc);
w=frac(1,r1.^2.*cos(inc));

function[ro,elo,azo,omo]=aquaparam
%AQUAPARAM  Parameters for the Aquarius satellite.
%
%   [RO,ELO,AZO,OMO]=AQUAPARAM returns parameters of the Aquarius satellite. 
%   RO is the nominal radial distance, in kilometers,  from the surface of  
%   the earth to the  Aquarius satellite.  ELO and AZO are the elevation 
%   angle (look angle) and azimuth angle of each of the three radiometer 
%   beams, in degrees.  OMO is the angular distance, in degrees, to the 
%   half-power point for each of the three beams.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006 J.M. Lilly --- type 'help jlab_license' for details    
          


ro=657;
elo=[25.8 33.8 40.3]';
azo=[9.8 -15.3 6.5]';
omo=[6.1 6.3 6.6]'/2;




