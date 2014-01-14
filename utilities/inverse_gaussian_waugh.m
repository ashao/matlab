function [ G ] = inverse_gaussian_waugh(parameters,t)
% Inverse Gaussian with parameters as defined in Waugh et al. 2003, JGR
    gamma = parameters(1);
    delta = parameters(2);
    
    a1 = sqrt ( gamma^3 ./ (4*pi*delta^2.*t.^3) );
    a2 = -gamma*(t-gamma).^2./(4*delta^2*t);
    G = a1 .* exp( a2 );
    G( t==0 ) = 0;
    
end