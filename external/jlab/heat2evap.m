function[E]=heat2evap(Q)
%HEAT2EVAP  Transform latent heat loss into units of evaporation.
%
%   E=HEAT2EVAP(Q) where Q is the latent heat loss from the ocean,
%   in units of Watts per square meter, returns the associated 
%   evaporation in units of millimeters per day.
%
%   Negative values of Q, corresponding to condensation, are 
%   returned as NaNs.
%
%   Usage: E=heat2evap(Q);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(Q, '--t')
    heat2evap_test,return
end

%      s^2/m^2     s / day   m^3 /kg    mm/m
c=frac(1,2.5*10^6)*3600*24*frac(1,1000)*1000;
E=Q*c;

function[]=heat2evap_test
 %   'heat2evap --t' runs a test.
 
%reporttest('HEAT2EVAP',aresame())
