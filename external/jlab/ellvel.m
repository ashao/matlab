function[vm,vbar,vphi,vi]=ellvel(varargin)
% ELLVEL  Average and instantaneous ellipse velocities. 
%
%   [VM,VBAR,VPHI,VI]=ELLVEL(KAPPA,LAMBDA,THETA,PHI) where KAPPA and LAMBDA
%   are the amplitude and linearity of a time-varying ellise, THETA is its
%   time-varying orientation, and PHI is its time-varying phase, returns
%   quantities related to the ellipse "velocity":
%
%       VM    Geometric mean velocity 
%       VBAR  Period-averaged speed 
%       VPHI  Instantaneous azimuthal velocity 
%       VI    Instantaneous speed
%
%   Note that all these quantities are defined to be positive when the
%   ellipse is orbited in the mathematically positive (counterclockwise) 
%   sense, and negative when the ellipse is orbited in the mathematically 
%   negative sense.
%   
%   ELLVEL(DT,...,) optionally uses DT as the data sample rate, with a 
%   default value of DT=1.  DT is a scalar.
%
%   ELLVEL(DT,...,FACT) optionally converts the physical units of velocity
%   through a multiplication by FACT, with a default value of FACT=1.  For 
%   example, FACT=1e5 converts kilometers into centimeters.
%
%   See Lilly and Gascard (2006) for details.
%   ____________________________________________________________________
%   
%   Cell array input/output
%
%   If ELLVEL is given cell array input, it returns cell array output.
%
%   Thus KAPPA, LAMBDA, THETA, and PHI may each be cell arrays of the 
%   same size, where each element in the cell array is a numerical array.
%
%   Then VM, VBAR, VPHI, and VI will be also cell arrays of this size.
%
%   In this case ELLVEL(DT,...) also works with DT a scalar or an 
%   array whose length is the number of elements in the cell arrays.
%   ____________________________________________________________________
%
%   See also ELLRAD, ELLPARAMS, ELLDIFF.
%
%   Usage:  vm=ellvel(kappa,lambda,theta,phi);
%           [vm,vbar,vphi,vi]=ellvel(kappa,lambda,theta,phi);
%           [vm,vbar,vphi,vi]=ellvel(dt,kappa,lambda,theta,phi);
%
%   'ellvel --f' generates a sample figure
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2005--2009 J.M. Lilly --- type 'help jlab_license' for details    



if strcmp(varargin{1}, '--f')
  ellvel_fig,return
end

if length(varargin{1})==1
    dt=varargin{1};
    varargin=varargin(2:end);
else 
    dt=1;
end
if length(varargin{end})==1
    fact=varargin{end};
    varargin=varargin(1:end-1);
else 
    fact=1;
end
k=varargin{1};
l=varargin{2};

theta=varargin{3};
phi=varargin{4};

if ~iscell(k)
    if nargout==1
        vm=ellvel_one(dt,k,l,theta,phi,fact);
    else
        [vm,vi,vphi,vbar]=ellvel_one(dt,k,l,theta,phi,fact);
    end
else
    for i=1:length(k)
        if isscalar(dt)
            dti=dt;
        else
            dti=dt(i);
        end
        if isscalar(fact)
            facti=fact;
        else 
            facti=fact{i};
        end
        if nargout==1
             vm{i}=ellvel_one(dti,k{i},l{i},theta{i},phi{i},facti);
        else
            [vm{i},vi{i},vphi{i},vbar{i}]=ellvel_one(dti,k{i},l{i},theta{i},phi{i},facti);
        end
        
    end
end


function[vm,vi,vphi,vbar]=ellvel_one(dt,k,l,theta,phi,fact)


str='endpoint';

omphi=vdiff(phi,1,str).*frac(1,dt).*fact;
[k2,l2,theta2,phi2]=elldiff(dt,k,l,theta,phi,fact);

if nargout ==1
    vm=sign(l).*ellrad(k2,l2,phi2);
else
    [vm,vbar,vi]=ellrad(k2,l2,phi2);

    vi=vi.*sign(l);
    vm=vm.*sign(l);
    vbar=vbar.*sign(l);
    
    z=ellsig(k,l,theta,phi);

    omphi=frac(1,dt).*vdiff(unwrap(angle(z)),1,str);
    [vm2,vbar2,vi2]=ellrad(k,l,phi);
    vphi=fact.*vi2.*omphi;
end

%These are equivalent methods
%zprime=rot(theta2).*(a2.*cos(phi2)+sqrt(-1).*b2.*sin(phi2));
%vphi=imag(rot(-angle(z)).*zprime);


function[]=ellvel_fig
lambda=(0:.001:1)';
kappa=1+0*lambda;
phi=(1:length(lambda))'/10;
theta=0*lambda;

[vm,vbar,vphi,vi]=ellvel(kappa,lambda,theta,phi);

figure
plot(lambda,[vbar,vm]./maxmax([vbar(:);vm(:)]));
linestyle k 2k 
e=[.2486 .967];
axis([0 1 0 1]),axis square 
vlines(e.^2./(2-e.^2),'k-')
axis([0 1 0.5 1]),axis square 
title('Velocity measures')
xlabel('Ellipse parameter \lambda')
ylabel('Mean velocity measures')
set(gcf,'paperposition', [2 2 3.5 3.5])
xtick(.1),ytick(.1),fixlabels(-1)
text(0.75,0.93,'V_{Bar}')
text(0.55,0.83,'V_M') 
fontsize 14 14 14 14
%fontsize jpofigure
%cd_figures
%print -depsc ellipsemeans.eps


