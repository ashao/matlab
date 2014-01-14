function[varargout]=jointfreq(varargin)
%JOINTFREQ  Joint instantaneous frequency, bandwidth, and curvature.
%
%   OMEGAX=JOINTFREQ(X,OMEGA,DIM) where X is an array of multiple analytic
%   signals and OMEGA is an array of their instantaneous frequencies, gives
%   the joint instaneous frequency OMEGAX averaged over dimension DIM.
%
%   X is presumed to have time oriented in rows.  The output matrix OMEGAX
%   are the same size as X except along dimension DIM, where the output
%   matrices will have length one.  X and OMEGA are the same size.
%   
%   [OMEGAX,UPSILONX]=JOINTFREQ(X,UPSILON,DIM) also works, where UPSILON 
%   are the individual bandwidths, and UPSILONX is the joint quantity. 
%   UPSILON is the same size as X and OMEGA.
%
%   Finally [OMEGAX,UPSILONX,XIX]=JOINTFREQ(X,UPSILON,XI,DIM) also returns 
%   the joint instantaneous curvature XIX given individual curvatures XI.
%
%   For details, see 
%  
%       Lilly & Olhede (2010), "Bivariate instantaneous frequency and 
%            bandwidth", IEEE Trans. Sig. Proc., 58 (2), 591--603.
% 
%   See also INSTFREQ, MOM2DEV.
%
%   Usage: om=jointfreq(x,om,dim);
%          [om,up]=jointfreq(x,om,up,dim);
%          [om,up,xi]=jointfreq(x,om,xi,dim);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009--2011 J.M. Lilly --- type 'help jlab_license' for details
 


if strcmp(varargin{1},'--t')
  jointfreq_test;return
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

varargout{1}=powermean(om,x,dim);
ombar=vrep(varargout{1},size(x,dim),dim);
   
if nargout>1
    varargout{2}=sqrt(powermean(abs(upsilon+sqrt(-1)*(om-ombar)).^2,x,dim));
end
if nargout>2
    varargout{3}=sqrt(powermean(abs(xi+2*sqrt(-1)*upsilon.*(om-ombar)-(om-ombar).^2).^2,x,dim));
end
if nargout>3
     error('Sorry, INSTFREQ only outputs the first two deviation vectors for joint moments.')
end

function[]=jointfreq_test
load solomon
use solomon

x=anatrans(x);
y=anatrans(y);

[omx,upx]=instfreq(x);
[omy,upy]=instfreq(y);

ombar=frac(abs(x).^2.*omx+abs(y).^2.*omy,abs(x).^2+abs(y).^2);
upbar=sqrt(frac(abs(x).^2.*(upx.^2+(omx-ombar).^2)+abs(y).^2.*(upy.^2+(omy-ombar).^2),abs(x).^2+abs(y).^2));

[om,up]=jointfreq([x y],[omx omy],[upx upy],2);


reporttest('JOINTFREQ using Solomon Islands frequency',aresame(om,ombar,1e-8))
reporttest('JOINTFREQ using Solomon Islands bandwidth',aresame(up,upbar,1e-8))
