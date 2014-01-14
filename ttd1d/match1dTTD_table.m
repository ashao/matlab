function [ gamma delta fval flag ] = match1dTTD( pCFC, sat, trac, meastime, varargin )
%% match1dTTD: Attempts to find an inverse Gaussian that when convolved
%  with the atmospheric source functon yields the bottle concentation
%   Input:
%       pCFC (ntrac): partial pressure of the tracers
%       sat: (ntrac,ntime): saturation level over time
%           OR sat( 1 ): specifies one level of saturation (e.g. 0.9 for
%           90%)
%       trac (struct):
%           atmconc (ntrac, ntime): Atmospheric source function
%           time (ntime): The time stamp of each value in atmconc
%       time (1): The year the measurement was made
%   Output:
%       gamma: The mean age of the fitted TTD
%       delta: The width parameter of the fitted TTD
%       fval: The fval in percent of each tracer (summed in quadrature)
%       flag: 0 if the routine is within the maximum fval allowed, 1 if
%           converged
%   Routine:
%       This routine attempts to match the analytic solution to the 1D linear
%       transport equation, an inverse Gaussian, with the observed
%       concentration and an atmospheric source scaled by the estimated
%       saturation level. This is done by creating sets of increasingly
%       finer resolution tables of Gamma, Delta values
%% Parse various options
i=1;
pCFC=double(pCFC);
while i <= length(varargin)
    
    if isstr(varargin{i})
        argin = upper(varargin{i});
        if strcmp(argin,'METHOD') % SIMULATEDANNEALING or PATTERNSEARCH
            method = varargin{i+1};
            i = i+2;
        elseif strcmp(argin,'MINPEC')
            minpecletsqr = sqrt(varargin{i+1});
            i = i+2;
        elseif strcmp(argin,'MAXPEC')
            maxpecletsqr = sqrt(varargin{i+1});
            i = i+2;p
        elseif strcmp(argin,'LOWERBOUND')
            lb = varargin{i+1};
            i = i+2;
        elseif strcmp(argin,'UPPERBOUND')
            ub = varargin{i+1};
            i = i+2;
        elseif strcmp(argin,'X0')
            x0 = varargin{i+1};
            i = i+2;
        elseif strcmp(argin,'MAXERR')
            maxerr = varargin{i+1};
            i = i+2;
        elseif strcmp(argin,'CONSTRAINT') % PECLET or FIXED
            confun = varargin{i+1};
            i = i+2;
        elseif strcmp(argin,'PECVAL') % For use with FIXED, defines a preset peclet ration
            pe = sqrt(varargin{i+1});
            i=i+2;
        elseif strcmp(argin,'MAXITER') % Maximum number of grid refinements
            maxiter = varargin{i+1};
            i=i+2;
        elseif strcmp(argin,'RES') % Number of points to use within a grid
            res = varargin{i+1};
            i=i+2;
        else
            error_str = ['Unknown identifier ''' varargin{i} '''!'];
            error(error_str);
        end
    else
        error('arg must be a string')
    end
end


%% Define some basic parameters of this routine, these are tunable

if ~exist('maxerr','var')
    maxerr = 1e-6; %Default 1e-6, set this pretty low, otherwise false convergence
end
if ~exist('minpecletsqr','var')
    minpecletsqr = 0.01; % Min ratio of delta/gamma,
end
if ~exist('maxpecletsqr','var')
    maxpecletsqr = 100.0; % Max ratio of delta/gamma
end
if ~exist('x0','var')
    x0=[100 100];
end
if ~exist('lb','var')
    lb = [0.1 0.1]; % Lower bound for mean age and width
end
if ~exist('ub','var')
    ub = [2000 2000]; % Upper bound for mean age and width
end
if ~exist('method','var')
    method = 'patternsearch';
end
if ~exist('confun','var')
    confun='none';
end
if ~exist('maxiter','var')
    maxiter=1000;
end
if ~exist('res','var')
    res=100;
end
if strcmpi(confun,'fixed') & ~exist('pe','var')
    pe = 1.0;
end

%% Get sizes of input variables
ntrac = length(pCFC);
% fprintf('Number of Tracers: %d\n',ntrac);

%% Scale the source function by estimated saturation
if length(sat)==1 % Do this if a constant value of saturation assumed
    sat=ones(size(trac.atmconc))*sat;
end
intime=trac.time<=meastime;

% Note we're flipping here to make sure the convolution with TTD is
% performed correctly

source_scaled = fliplr(sat(:,intime).*trac.atmconc(:,intime));
source_scaled(isnan(source_scaled))=0;
time_scaled = meastime - trac.time(intime); % Calculate as time since measurement
time_scaled = flipud(time_scaled)';

iter = 0;
tracerr = Inf;
range = 1;
while tracerr > maxerr & iter < maxiter
    
    range = range*0.9;
    iter = iter+1;
    
    if iter > 2
        [gammagrid deltagrid] = makegrid(x0, res,range);
    else
        [gammagrid deltagrid] = meshgrid(linspace(0.1, 1000, 1000));
    end
        
    
    if strcmpi(confun,'fixed')
        deltagrid = sqrt(gammagrid.^2./pe);
        gammagrid = [diag(gammagrid) ; x0(1)];
        deltagrid = [diag(deltagrid) ; x0(2)];
        
        
    end
    Gconv = convTTDsource(gammagrid(:), deltagrid(:),time_scaled,source_scaled);    
    [tracerr bestidx] = min(abs(Gconv-pCFC));
    gamma = gammagrid(bestidx);
    delta = deltagrid(bestidx);
    x0 = [gamma delta];
    fprintf('Iter: %d Gamma:%f Delta:%f Error: %f\n',iter,gamma,delta,tracerr)
end

gamma = gammagrid(bestidx);
delta = deltagrid(bestidx);
fval = tracerr;
if tracerr > maxerr
   flag = false;
else
    flag = true;
end

end
%% INDEPENDENT FUNCTIONS
function [ gammagrid deltagrid ] = makegrid(x0,res,range)

bounds = [1-range 1+range];

gammarange = x0(1)*bounds;
deltarange = x0(2)*bounds;
gamma = linspace(gammarange(1),gammarange(2),res);
delta = linspace(deltarange(1),deltarange(2),res);
[gammagrid deltagrid] = meshgrid(gamma,delta);

end

function [ Gconv ] = convTTDsource(gammagrid,deltagrid,time,source)

timeconv = linspace(0,max(time),1000);
sourceconv = interp1(time,source,timeconv);
npts = length(gammagrid);
Gconv = zeros(npts,1);

for pt = 1:npts
    G = inverse_gaussian([gammagrid(pt) deltagrid(pt)],timeconv);
    Gconv(pt) = trapz(timeconv,G.*sourceconv);
end

end

function [ G ] = inverse_gaussian(parameters,t)
% Inverse Gaussian with parameters as defined in Waugh et al. 2003, JGR
gamma = parameters(1);
delta = parameters(2);

a1 = sqrt ( gamma.^3 ./ (4*pi*delta.^2.*t.^3) );
a2 = -gamma*(t-gamma).^2./(4*delta.^2*t);
G = a1 .* exp( a2 );
G( t==0 ) = 0;

end