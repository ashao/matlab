function tracprop = tracer_properties ( )
%% Make structure with solubilities and Schmidt number based on Warner and E

% CFC Schmidt number coeffs from Zheng et al. 1985
tracprop.cfc11.sc_coeffs=[3501.8 -210.31 6.1851 -0.07513];
tracprop.cfc12.sc_coeffs=[3845.4 -228.95 6.1908 -0.067430];

% SF6 Schmidt number (Wanninkhof 1992)
tracprop.sf6.sc_coeffs=[3531.6 231.40 7.2168 0.090558];

% CFC Solubilities Warner and Weiss 1985
tracprop.cfc11.sol_coeffs.vol=...
    [-229.9261 319.6552 119.4471 -1.39165 -0.142382 0.091459 -0.0157274];
tracprop.cfc11.sol_coeffs.grav=...
    [-232.0411 322.5546 120.4956 -1.39165 -0.146531 0.093621 -0.0160693];
tracprop.cfc12.sol_coeffs.vol=...
    [-218.0971 298.9702 113.8049 -1.39165 -0.143566 0.09105 -0.0153924];
tracprop.cfc12.sol_coeffs.grav=...
    [-220.2120 301.8695 114.8533 -1.39165 -0.147718 0.093175 -0.0157340];

% SF6 Solubility Bullister 2002 NOTE: EXTRA ZERO ADDED BECAUSE SOLUBILITY
% FUNCTION IS DOES NOT USE A (T/100)^2 term
tracprop.sf6.sol_coeffs.vol=...
    [-80.0343 117.232 29.5817 0 0.0335183 -0.0373942 0.00774862 ]; 
tracprop.sf6.sol_coeffs.grav=...
    [-82.1639 120.152 30.6372 0 0.0293201 -0.0351974 0.00740056];
end