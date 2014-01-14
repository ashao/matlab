function[g,T,km,sm]=gc_params
%GC_PARAMS  Parameters for gravity-capillary waves.
%
%   [G,T,K/M,SM]=GC_PARAMS returns gravity (G), surface tension
%   divided by density (T), and the wavenumber (KM) and frequency (SM) 
%   at which the phase speed is minimized, in cgs units.  
%
%   Note that KM=SQRT(G/T) and SM=(4G^3/T)^(1/4) assume the standard
%   dispersion relation S^2=GK+TK^3.
%
%   The units are: 
%        G [cm s^-2], T [cm^-1 s^-2], KM [rad cm^-1], SM [rad s^-1]
%
%   See also OM, KMIN.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001--2006 J.M. Lilly --- type 'help jlab_license' for details    
  

rho=1000/(100^3);
g=981;
T=0.07288/rho; %Gives T/rho=72.88 as in McGoldrick

km=sqrt(g/T);
sm=(4*g^3/T).^(1/4);
