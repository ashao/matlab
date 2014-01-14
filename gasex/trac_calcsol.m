function [ trac_sol ] = trac_calcsol(TempK,S,sol_coeffs)
% CALCSOL: Calculate solubility of a gas
% INPUTS:
%   T: Temperature in Kelvin
%   S: Salinity in PSU
%   sol_coeffs: Coefficients for empirical formula

trac_sol = sol_coeffs(1) + sol_coeffs(2)*(100 ./ TempK)+ ...
    sol_coeffs(3) * log(TempK ./ 100) + sol_coeffs(4) * (TempK ./ 100) ...
    + S .* (sol_coeffs(5) + sol_coeffs(6) * (TempK ./ 100) ...
    + sol_coeffs(7) * (TempK ./ 100).^2);
trac_sol = exp(trac_sol);

end