function [ ar n2 o2 ] = ar_n2_o2_properties( )
%% Sets a bunc hof gas properties for Ar O2 N2

ar.F_sol_coeffs.vol = [-173.5146 245.4510 141.8222 -21.8020 -0.034474 ...
    0.014934 -0.0017729];
n2.F_sol_coeffs.vol = [-172.4965 248.4262 143.0738 -21.7120 -0.049781 ...
    0.025018 -0.0034861];
o2.F_sol_coeffs.vol = [-173.4292 249.6339 143.3483 -21.8492 -0.033096 ...
    0.014256 -0.0017000];
ar.atmconc = 0.00934;
n2.atmconc = 0.780840;
o2.atmconc = 0.209476;

end