function [ gamma delta fval flag Gconv ] = match1dTTD( pCFC, sat, trac, meastime, varargin )
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
%       saturation level. This is done by using simulated annealing to
%       minimize the cost function ( (conv(TTD,tracsource)-pCFC)./pCFC ).^2
%% Parse various options

i=1;
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
            i = i+2;
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
    x0=randi([10,500],2,1); % Initial guess, doesn't seem to have much effect
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
if strcmpi(confun,'fixed') & ~exist('pe','var')
    pe = 1.0;
end
%% Get sizes of input variables
ntrac = length(pCFC);
Gconv = zeros(1,ntrac);
% fprintf('Number of Tracers: %d\n',ntrac);

%% Scale the source function by estimated saturation
if length(sat)==1 % Do this if a constant value of saturation assumed
    sat=ones(size(trac.atmconc))*sat;
end
intime=trac.time<meastime;

% Note we're flipping here to make sure the convolution with TTD is
% performed correctly

source_scaled = fliplr(sat(:,intime).*trac.atmconc(:,intime));
% source_scaled(isnan(source_scaled))=0;
time_scaled = meastime - trac.time(intime); % Calculate as time since measurement
time_scaled = flipud(time_scaled)';

% pCFC
% pCFC
% source_scaled
% % whos
% pause
%
%% Perform the minimization using two different techniques
switch lower(method)
    case 'simulatedannealing'
        options=saoptimset('TolFun',maxerr,'Display','none', ...
            'HybridFcn',@fmincon);
        [optpar fval] = simulannealbnd(@evalTTD,x0,lb,ub,options);
        
    case 'patternsearch'
        options=psoptimset('TolFun',sqrt(maxerr),'Display','none', ...
            'UseParallel','never','InitialMeshSize',.1,...
            'CompletePoll','off','MaxFunEvals',inf,'MaxIter',10000, ...
            'CompleteSearch','off');
        
        if strcmpi(confun,'peclet') % Constraints on Peclet Number            
            [optpar fval] = patternsearch(@evalTTD_peclet, ...
                x0,[],[],[],[],lb,ub,@pecconstraints,options);
            gamma = optpar(1);
            delta = optpar(2);
        elseif strcmpi(confun,'fixed') % Fixed Peclet Number
            [optpar fval] = patternsearch(@evalTTD_fixed, ...
                x0,[],[],[],[],lb,ub,options);
            gamma=optpar;
            delta=gamma/pe;
        elseif strcmpi(confun,'none') % No constraints on Peclet
            [optpar fval] = patternsearch(@evalTTD_peclet, ...
                x0,[],[],[],[],lb,ub,options);
            gamma = optpar(1);
            delta = optpar(2);
        end
        
        
        fval=sqrt(fval);
        if fval < maxerr
            flag=true;
        else
            flag=false;
        end
end
%% NESTED FUNCTIONS
% Next two functions, takes the parameters from the minimization routine
% concolves TTD with source and calculates the fval associated with it

% For use with CONSTRAINT=PECLET or NONE, allows gamma/delta to vary
    function [ fval ] = evalTTD_peclet( parameters )
        
        G=inverse_gaussian(parameters,time_scaled);
        
        for tracidx=1:ntrac            
            Gconv(tracidx) = trapz(time_scaled,source_scaled(tracidx,:).*G);
        end
        
        fval = sum( (1-Gconv./pCFC).^2);
    end

% For use with CONSTRAINT=FIXED, assumes delta = gamma/Pe
    function [ fval ] = evalTTD_fixed( gamma )
        
        G=inverse_gaussian([gamma gamma/pe],time_scaled);
        for tracidx=1:ntrac
            Gconv(tracidx) = trapz(time_scaled,source_scaled(tracidx,:).*G);
        end
        
        fval = sum( (1-Gconv./pCFC).^2);
        
    end

    function [c,ceq,gradc,gradceq] = pecconstraints( parameters )
        x1 = parameters(1);
        x2 = parameters(2);
        pecletsqr = x1/x2;
        
        c(1) = pecletsqr - maxpecletsqr;
        c(2) = minpecletsqr-pecletsqr;
        ceq = [];
        if nargout > 2
            gradc = [1/x2 -x1/x2^2; -1/x2 x1/x2^2];
            gradceq = 0;
        end
        c
    end

end

%% INDEPENDENT FUNCTIONS
function [ G ] = inverse_gaussian(parameters,t)
% Inverse Gaussian with parameters as defined in Waugh et al. 2003, JGR
gamma = parameters(1);
delta = parameters(2);

a1 = sqrt ( gamma^3 ./ (4*pi*delta^2.*t.^3) );
a2 = -gamma*(t-gamma).^2./(4*delta^2*t);
G = a1 .* exp( a2 );
G( t==0 ) = 0;

end
