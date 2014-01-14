function[eke]=ssh2eke(lat,atd,ssh,N)
%SSH2EKE  Converts alongtrack sea surface height to eddy kinetic energy.
%
%   EKE=SSH2EKE(LAT,ATD,SSH,N) computes the eddy kinetic energy from 
%   the altimeter sea surface height anomaly SSH, smoothed alongtrack 
%   with an N-point Hanning filter.
%
%   LAT is the latitude, and ATD is the along-track distance in 
%   kilometers.  LAT and ATD the same size as SSH(:,:,1).
%
%   You should make sure that SSH has the mean sea surface height 
%   removed in order to avoid a contribution from aliasing the geoid.
%
%   SSH2EKE assumes that the contribution to EKE from the cross-track
%   and along-track velocity anomalies are equal, that is, that 
%   the velocity anomaly field is isotropic.
%
%   Lilly et. al (2003) use a smoothing of N=5 in the Labrador Sea.
%
%   Usage: eke=ssh2eke(lat,atd,ssh,N);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2008 J.M. Lilly --- type 'help jlab_license' for details
 

%disp('Beware! change')
dx=1000*abs(vdiff(atd,1));

%dx=1000*abs(atd);

c=9.81./corfreq(lat)./dx/3600;
%c=9.81./sw_f(lat)./dx;
c=vrep(c,size(ssh,3),3);

%ssh=ssh-vrep(mssh,size(ssh,3),3);
%The absence of a factor of 1/2 is because we assume isotropy---
%i.e. we only have one of two orthogonal directions, so we have 
%to count it twice.

if N>1
    eke=squared(c.*vdiff(vfilt(ssh,N),1));
else
    eke=squared(c.*vdiff(ssh,1));
end
