function theta_mean = mean_angle( theta, range )
%% theta_mean = mean_angle ( theta, range )
% Compute the mean from a set of angles
% INPUT:
%   theta (required): Vector of thetas from which to calcualte the mean
%   range (optional) string:
%       '2pi': Output goes from [0 2Pi]
%       'pi': Output from [-Pi Pi] (default)
% OUTPUT:
%   theta_mean: Mean of the input thetas output in the specified range
% NOTE: theta in radians
% Algorithm:
%   1) Assume the angles are phases of some complex number
%   2) Calculate mean complex number
%   3) Calculate phase of this mean complex number
% April 10, 2012

if nargin < 2
    range='pi';
end
range=lower(range);

theta_mean=angle( nanmean ( exp(1i*theta )) );

switch range
    case '2pi'
        theta_mean(theta_mean<0)=theta_mean+2*pi;
end


end