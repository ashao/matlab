function[varargout]=mom2dev(varargin)
%MOM2DEV  Converts instantaneous moments into deviation vectors.
%
%   D1=MOM2DEV(X,OMEGA,UPSILON,DIM) returns the first deviation vector X1 
%   given a multivariate signal X, the components of which have 
%   instantaneous frequencies OMEGA and bandwidths UPSILON.  
%
%   X, UPSILON, and D1 are all arrays of the same size with the first 
%   dimension being "time".  DIM is the dimesional along which these arrays
%   contain multiple signals.
%
%   [D1,D2]=MOM2DEV(X,OMEGA,UPSILON,XI,DIM) where XI is the instantaneous 
%   curvature of the signals also returns D2, the second deviation vector.
%
%   Note that OMEGA, UPSILON, and XI are the moments of the individual
%   signal component and not the joint moments of the multivariate signal.
%
%   See also INSTFREQ, JOINTFREQ.
%
%   Usage: d1=mom2dev(x,omega,upsilon,dim);
%          [d1,d2]=mom2dev(x,omega,upsilon,xi,dim);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2011 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    mom2dev_test,return
end
 
dim=varargin{end};
x=varargin{1};
varargin=varargin(2:end-1);

om=varargin{1};
if length(varargin)>1
    upsilon=varargin{2};
end
if length(varargin)>2
   xi=varargin{3};
end

ombar=vrep(powermean(om,x,dim),size(x,dim),dim);
   
varargout{1}=x.*(upsilon+sqrt(-1)*(om-ombar));
if nargout>1
    varargout{2}=x.*(xi+2*sqrt(-1)*upsilon.*(om-ombar)-(om-ombar).^2);
end
if nargout>2
     error('Sorry, MOM2DEV only outputs the first two deviation vectors for joint moments.')
end


function[]=mom2dev_test
 
%reporttest('MOM2DEV',aresame())
