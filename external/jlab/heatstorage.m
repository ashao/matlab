function[Q,H,Qstar,sstar,thstar,sbar,thbar,ds,dth]=heatstorage(s,th,p,gamma)
%HEATSTORAGE Water column heat storage from 1-D mixing.
%
%   HEATSTORAGE implements the one-dimensional mixing model described in
%
%       Lilly et. al (2003), "Observations of the Labrador Sea eddy field"
%       Progress in Oceanography
%       Appendix B: One-dimensional mixing and oceanic heat storage
%
%   It is used for determining the mixed layer depth and water mass
%   properties created by surface heat loss and subsequent mixing.
%
%   Note that HEATSTORAGE does not take into account the possibility
%   of freezing, which may arrest the deepening of the mixed layer.
%
%   [Q,H,QSTAR,SSTAR,THSTAR,SBAR,THBAR]=HEATSTORAGE(S,TH,P,GAMMA)
%
%   Input variables 
%     S:      Salinity column vector or matrix
%     TH:     Potential temperatute column vector or matrix
%     P:      Pressure column vector or matrix  
%     GAMMA:  GAMMA=dS/dTH is the ratio of water column salinity gain
%             (dS) to heat gain (dTH) due to sea surface fluxes.
%
%		             GAMMA=rho_o * c_p * S_o * (P-E)/Q 
%
%             where rho_o is a reference density, c_p is the specific
%             heat of water, S_o is a reference salinity, P-E is net
%             precipitation, and Q is total heat loss assuming
%             radiation and precipitation heat fluxes are negligable.
%
%   Output variables
%     Q:      Heat content in 10^9 Joules / m^2
%     H:      Available heat content in 10^9 Joules / m^2; H_i is the 
%             surface heat loss it takes to convect to depth #i
%     QSTAR:  "Unavailable" Heat content in 10^9 Joules / m^2; QSTAR_i 
%             is the heat content of the convected mixed layer whose 
%             base is at depth #i
%     SSTAR:  SSTAR_i is the salinity of the mixed layer whose base is 
%             at depth #i
%     THSTAR: Same as SSTAR but for potential temperature
%     SBAR:   SBAR_i is the salinity of original water column averaged 
%             from the surface to depth #i
%     THBAR:  Same as SBAR but for potential temperature
%
%   Note that HEATSTORGE requires the "SEAWATER" Matlab toolbox by
%   Phillip Morgan of CSIRO. 
%  
%   Usage: [q,h,qstar,sstar,thstar,sbar,thbar]=heatstorage(s,th,p,gamma)
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2000, 2004 J.M. Lilly --- type 'help jlab_license' for details  
  
%Specific heat of water defined as 4000

%output: in giga (1e9) Joules
%note 500 W/m2 for 90 days ~ 4 GJ
 
  
t=sw_temp(s,th,p,0);
rho=sw_dens(s,t,p);

dp=diff([0*p(1,:);p]);

%Set values of zero to a small number
index=find(p==0);
if ~isempty(index)
  p(index)=0.1;
end

sbar=cumsum(s.*dp)./p;
thbar=cumsum(th.*dp)./p;

tbar=sw_temp(sbar,thbar,p,0);
rhobar=sw_dens(sbar,tbar,p);

alpha=sw_alpha(sbar,tbar,p);
beta=sw_beta(sbar,tbar,p);

alpha=cumsum(alpha.*dp)./p;
beta=cumsum(beta.*dp)./p;

%alpha=sw_alpha(34.85,3,1000);
%beta=sw_beta(34.85,3,1000);

dth=(rhobar-rho)./(1000.*alpha.*(1-beta./alpha.*gamma));
ds=gamma*dth;

thstar=thbar+dth;
sstar=sbar+ds;

H=-1000*4e3*dth.*p/1e9;
Q=1000*4e3*p.*thbar/1e9;
Qstar=1000*4e3*p.*thstar/1e9;













