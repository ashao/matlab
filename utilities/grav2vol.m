function [ con_vol ] = grav2vol ( con_grav, rho )
% GRAV2VOL Converts gravimetric concentrations to volumetric units
% Input:
%	con_grav: Concentration (#mol/kg)
%	rho: density of seawater (kg/m^3)
% Output:
%	con_vol: Concentration (#mol/L)

rho=rho/1000; %Convert from kg/m^3 to kg/L
con_vol=con_grav.*rho;

end

