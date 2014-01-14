function[fc]=corfreq(lat)
%CORFREQ  Coriolis frequency in cycles per hour.
%
%   FC=CORFREQ(LAT) returns the Coliolis frequency at latitude LAT
%   in cycles per hour.  FC is defined to be non-negative.
%
%   'corfreq --t' runs a test.
%
%   Usage: fc=corfreq(lat);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(lat, '--t')
    corfreq_test,return
end
 

omega=7.292e-5;
fc=frac(1,2*pi)*2*abs(sind(lat)).*omega.*(3600);

function[]=corfreq_test
 
reporttest('CORFREQ at 30 degrees is about 1/24',aresame(1./corfreq(30),24,1e-1))
