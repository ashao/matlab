function[s]=kmin
%KMIN  Wavenumber of minimum phase speed for gravity-capillary waves.
%
%   K=KMIN with no arguments returns the wave number K of the 
%   minimum phase speed, equal to SQRT(G/T), with units [cm s^-1].
%
%   See also GC_PARAMS, OM.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001--2006 J.M. Lilly --- type 'help jlab_license' for details    

[g,T,km,sm]=gc_params;  %cgs units
s=sqrt(g./T);

