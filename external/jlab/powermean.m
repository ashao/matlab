function[varargout]=powermean(varargin)
%POWERMEAN  Power-weighted mean along a specified dimension.
%
%   FM=POWERMEAN(F,X,DIM) takes the mean of all finite elements of F along
%   dimension DIM, weighted by the squared magnitude of X.
%   
%   The power-weighted mean of F is defined as 
%
%        FM = SUM (ABS(X)^2.*F, DIM) / SUM(ABS(X)^2, DIM)
%
%   where X and F are arrays of the same size.  
%                                                                         
%   [FM1,FM2,...FMN]=POWERMEAN(F1,F2,...FN,X,DIM) also works.
%
%   OWERMEAN(F1,F2,...FN,X,DIM);  with no output arguments overwrites the 
%   original input variables.
%
%   POWERMEAN with X a set of analytic signals is used to construct the
%   joint instantaneous moments.  For details on these quantities, see
%
%       Lilly and Olhede (2010).  Bivariate instantaneous frequency and
%           bandwidth.  IEEE Trans. Sig. Proc., in press.
%
%   See also VMEAN.
%
%   Usage: fm=powermean(f,x,dim);
%          [fm1,fm2]=powermean(f1,f2,x,dim);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009--2011 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    powermean_test,return
end

dim=varargin{end};
x=varargin{end-1};
power=vsum(abs(x).^2,dim);

vswap(power,0,nan);
for i=1:length(varargin)-2
  varargout{i}=vsum(abs(x).^2.*varargin{i},dim)./power;
end

eval(to_overwrite(nargin-2))
 
function[]=powermean_test

reporttest('POWERMEAN',aresame(powermean([1 2],[2 3],2),frac(4+2*9,13)))
