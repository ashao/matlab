function[rm,rbar,ri]=ellrad(k,l,phi)
%ELLRAD  Average and instantaneous ellipse radius. 
%
%   [RM,RA,RI]=ELLRAD(KAPPA,LAMBDA,PHI) where KAPPA and LAMBDA are the 
%   amplitude and linearity of a time-varying ellise, and PHI is its time-
%   varying phase, returns quantities related to the ellipse 'radius',
%   i.e. the distance from the ellipse curve to the origin:
%
%       RM    Geometric-mean radius
%       RBAR  Period-averaged distance from the origin
%       RI    Instantaneous distance from the origin
%  
%   See Lilly and Gascard (2006) for details.
%   ____________________________________________________________________
%   
%   Cell array input/output
%
%   If ELLRAD is given cell array input, it returns cell array output.
%
%   Thus KAPPA, LAMBDA, and PHI may each be cell arrays of the same size, 
%   where each element in the cell array is a numerical array.
%
%   Then RM, RA, and RI will be also cell arrays of this size.
%   ____________________________________________________________________
%
%   See also ELLVEL, ELLPARAMS, ELLDIFF.
%
%   Usage:  rm=ellrad(kappa,lambda);
%           [rm,ri,rbar]=ellrad(kappa,lambda,phi);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2005--2009 J.M. Lilly --- type 'help jlab_license' for details    
  
if ~iscell(k)
    [a,b]=kl2ab(k,l);
    rm=sqrt(a.*abs(b));
    if nargout>1
        [ri,rbar]=ellrad_one(a,b,phi);
    end
else
    for i=1:length(k)
        [a,b]=kl2ab(k{i},l{i});
        rm{i}=sqrt(a.*abs(b));
        if nargout>1
            [ri{i},rbar{i}]=ellrad_one(a,b,phi{i});
        end
    end
end

function[ri,rbar]=ellrad_one(a,b,phi)

[kappa,lambda]=ab2kl(a,b);
ecc=ecconv(lambda,'lin2ecc');
[K,E]=ellipke(ecc.^2);

rbar=frac(2*kappa,pi).*sqrt(1+abs(lambda)).*E;
ri=kappa.*sqrt(1+abs(lambda).*cos(2*phi));


%   ra=kappa.*(1-lambda.^2/16);
%   RA is obtained from a power series expansion in terms of the
%   eccentricity parameter lambda. 