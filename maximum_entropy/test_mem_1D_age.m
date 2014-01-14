function [ G, prior, actual_ttd ] = test_mem_1D_age( prior_type )
%% test_mem_1D_age(prior_type)
%Demonstrates the maximum entropy method to deconvolve Green's function G
% whose exact solution is a 1D Inverse Gaussian TTD with gamma=delta=240
% To get interior concentration, this is convolved with a simple
% age tracer, Csurf(t)=t.
% REFERENCE:
%   Holzer et al. 2010, JGR
%       Where and how long ago was water in the western North Atlantic
%       ventilated? Maximum entropy inversions of bottle data from WOCE
%       line A20
% INPUT:
%   prior_type: Two defined prior distributions
%                   1) 'uniform' Uniform Distribution
%                   2) 'inverse_gaussian' Inverse Gaussian PDF with
%                      gamma=delta=120
% OUTPUT:
%   G:  Maximum entropy method derived TTD
%   prior: Value of the prior used
%   actual_ttd: The actual inverse_gaussian TTD used to derive interior age

if nargin == 0
    prior_type='uniform';
end

mean_age = 240; % Define actual mean age
time = 1:10000; % Define time over which to convolve

% Actual Green's function for a hypothetical 1D Transport Equation
actual_ttd=inverse_gaussian([mean_age,mean_age],1:10000);

% Define a surface boundary condition at each timestep
surf_bc = 1:10000;

% Determine what kind of prior to use
switch lower(prior_type)
    case 'uniform'
        prior = repmat(1/trapz(time),[1 length(time)]); %% Uniform Prior
    case 'inverse_gaussian'
        % Inverse Gaussian with delta=gamma=120
        prior = inverse_gaussian([120 120],1:10000);
end

% Solve MEM Problem using method of Lagrange Multipliers (Last Paragraph
% Section 2.1
options=optimset('Display','None');
lambda_opt = fsolve(@constraints,rand(1),options);

% Calculate MEM-derived Green's Function (Equation 8)
G=prior/calcZ(lambda_opt).*exp(-lambda_opt.*surf_bc);

% Calculate first moment of TTD
mean_age=trapz(G.*time);

fprintf('Lagrange Multiplier: %f\n',lambda_opt)
fprintf('Mean Age of MEM Greens Function: %f\n',mean_age)

%% NESTED SUBFUNCTIONS

    function [ cost ] = constraints( lambda )
        % The constraint equation substituting Eq. 8 into Eq. 2
        partition = calcZ( lambda ); 
        cost = trapz( prior./partition.*exp(-lambda*surf_bc).*surf_bc ) ...
            -mean_age;        
    end

    function [ Z ] = calcZ( lambda )
        %% Calculate partition function (Equation 9)
        Z=trapz(prior.*exp(-lambda.*surf_bc));
    end

end

%% INDEPENDENT SUBFUNCTIONS
function [ dist ] = inverse_gaussian( coeffs, time )
gamma=coeffs(1);
delta=coeffs(2);

dist=sqrt(gamma^3./(4*pi*delta^2*time.^3)).*...
    exp(-gamma*(time-gamma).^2./(4*delta^2*time));
end