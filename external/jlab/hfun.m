function[H]=hfun(varargin)
% HFUN  Sea surface height interaction functions for wave triads.
%
%   HFUN(I,K1,K2) returns the i-th H-function (H_i) operating on
%   wavenumeber arrays K1 and K2, defined in Lilly and Lettin (2006)
%   and based on Elfouhaily et al. (2000).
%  
%   HFUN(S1,S2,K1,K2) also works, where S1 and S2 are each either 
%   positive or negative one.  See SS2I for the relation between these
%   two conventions.  
%
%   The H-functions come in two variations.  HFUN(...,'tilde') returns
%   the interaction functions resulting from expanding the potential 
%   about the sea surface.  This is the default if no string is input.
%   HFUN(...,'notilde') returns the interaction functions resulting 
%   from expanding the potential about zero.
%   
%   Note that several typographic errors from Elfouhaily et al. (2000) 
%   have been corrected, as discussed in Lilly and Lettin (2006).
%
%   The wavenumber arrays are in complex form, that is, K = Kx + i*Ky,
%   and have units of [rad cm^-1].  Values for gravity and surface 
%   tension are specified in GC_PARAMS, and the dispersion relation
%   is given in OM.  
%
%   See also: DFUN, HFUN, RESCOEFF. 
%
%   Usage: H=hfun(i,k1,k2);
%          H=hfun(s1,s2,k1,k2);
%  
%   'hfun --t' runs some tests.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001--2006 J.M. Lilly --- type 'help jlab_license' for details    


if strcmp(varargin{1},'--t')
  hfun_test;return
end
  
str='tilde'; 
str2='elf';
na=nargin;  
if ischar(varargin{end})
  na=na-1; 
  str=varargin{end};
  varargin=varargin(1:end-1);
end
if ischar(varargin{end})
  na=na-1; 
  str2=str;
  str=varargin{end};
  varargin=varargin(1:end-1);
end

if na==3
  jj=varargin{1};
  k1=varargin{2};
  k2=varargin{3};
  
  [s1,s2]=i2ss(jj);
elseif na==4
  s1=varargin{1};
  s2=varargin{2};
  k1=varargin{3};  
  k2=varargin{4};  
end

[g,T]=gc_params;
  
if strcmp(str2(1:3),'elf')
   D=-sqrt(-1)*dfun(s1,s2,k1,k2,str);
   Ha=2*frac((s1.*om(k1)+s2.*om(k2)).*D,om(k1+k2).^2-(s1.*om(k1)+s2.*om(k2)).^2);
   Hb=-cdot(k1,k2);
   Hfact=frac(1,g+T.*abs(k1+k2).^2).*frac(s1.*om(k1).*s2.*om(k2),abs(k1).*abs(k2));
   if  strcmp(str,'notilde')
     Hc=(1+frac(s1.*om(k1),s2.*om(k2))+frac(s2.*om(k2),s1.*om(k1))).*abs(k1).*abs(k2);
   elseif  strcmp(str,'tilde')
     Hc=-abs(k1).*abs(k2); 
   end
   H=frac(1,2)*Hfact.*(Ha+Hb+Hc); 
   HP=frac(1,2)*Hfact.*Ha;
   HB=frac(1,2)*Hfact.*(Hb+Hc);
elseif strcmp(str2(1:3),'hol')
   Ha=-frac(kfun(s1,s2,s1.*k1,s2.*k2,str),om(k1+k2)-s1.*om(k1)-s2.*om(k2));
   Hb=-frac(kfun(-s1,-s2,s1.*k1,s2.*k2,str),om(k1+k2)+s1.*om(k1)+s2.*om(k2));
   H=Ha+Hb; 
end

function[]=hfun_test
  
kg=wavegrid(.126,32)*kmin;
k=real(kg(1,:))';

k1=2*kmin+sqrt(-1)*1e-10;
  
for i=1:4
 H1(:,:,i)=hfun(i,k1,kg,'tilde','elf');
 H2(:,:,i)=hfun(i,k1,kg,'tilde','hol');
 H3(:,:,i)=hfun(i,k1,kg,'notilde','elf');
 H4(:,:,i)=hfun(i,k1,kg,'notilde','hol');
end

tol=1e-6;
reporttest('HFUN Holliday and Elfouhaily tilde versions match at 2*kmin',aresame(H1,H2,tol))
reporttest('HFUN Holliday and Elfouhaily no-tilde versions match at 2*kmin',aresame(H3,H4,tol))

k1=4*kmin+sqrt(-1)*1e-10;
  
for i=1:4
 H1(:,:,i)=hfun(i,k1,kg,'tilde','elf');
 H2(:,:,i)=hfun(i,k1,kg,'tilde','hol');
 H3(:,:,i)=hfun(i,k1,kg,'notilde','elf');
 H4(:,:,i)=hfun(i,k1,kg,'notilde','hol');
end

tol=1e-6;
reporttest('HFUN Holliday and Elfouhaily tilde versions match at 4*kmin',aresame(H1,H2,tol))
reporttest('HFUN Holliday and Elfouhaily no-tilde versions match at 4*kmin',aresame(H3,H4,tol))
