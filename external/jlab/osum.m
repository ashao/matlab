function[z]=osum(x,y)
%   OSUM  "Outer" sum of two column vectors.  
%        OSUM(X,Y) <==> X*(1+0*Y')+ (1+0*X)*conj(Y)'
%
%   See also OPROD, ODOT
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2003--2011 J.M. Lilly --- type 'help jlab_license' for details        
  
if ~(jiscol(x) || isscalar(x)) || ~(jiscol(y) || isscalar(y))
  error('X and Y must both be column vectors or scalars.')
else
  z=x*(1+0*y')+(1+0*x)*conj(y)';  
end

