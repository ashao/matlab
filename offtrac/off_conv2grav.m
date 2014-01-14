function [ conc_kg ] = off_conv2grav( conc_L, T, S, depth ) 
% OFF_CONV2GRAV Converts volumetric output (mol/L) from offtrac to
% gravimetric (mol/kg)
% -------------------------------------------------------------------------
% Dependencies:
%   sw_dens.m from seawater package
% -------------------------------------------------------------------------
% Input:
%   conc_L: concentration in (m,p,f)mol/L
%   T: Temperature in degrees C
%   S: Salinity in PSU
% Output:
%   conc_kg: concentration in (m,p,f)mol/kg

if nargin<4
    depth=0;
    if ndims(T)>2
        T=squeeze(T(1,:,:));
    end
    
    if ndims(S)>2
        S=squeeze(S(1,:,:));
    end
end

% Convert to mol/m^3
conc_m3=conc_L*1000;
% Calculate density of seawater packet
rho=sw_dens(T,S,depth);
% Convert to gravimetric
conc_kg=conc_m3./rho;


end
