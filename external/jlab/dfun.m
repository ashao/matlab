function[D]=dfun(varargin)
%DFUN	Velocity potential interaction functions for wave triads.
%
%   DFUN(I,K1,K2) returns the i-th D-function (D_i) operating on
%   wavenumber arrays K1 and K2, defined in Lilly and Lettin (2006) 
%   and based on Holliday (1977) and Elfouhaily et. al (2000).  
%  
%   DFUN(S1,S2,K1,K2) also works, where S1 and S2 are each either 
%   positive or negative one.  See SS2I for the relation between these
%   two conventions.  
%
%   The D-functions come in two variations.  DFUN(...,'tilde') returns
%   the interaction functions resulting from expanding the potential 
%   about the sea surface.  This is the default if no string is input.
%   DFUN(...,'notilde') returns the interaction functions resulting 
%   from expanding the potential about zero.
%
%   Additionally, DFUN(...,'limit') returns the limiting version of 
%   either the 'tilde' or 'notilde'  D-function for the case of 
%   negligible surface tension.  
%
%   Note that several typographic errors from Elfouhaily et al. (2000) 
%   have been corrected, as discussed in Lilly and Lettin (2006).
%
%   The wavenumber arrays are in complex form, that is, K = Kx + i*Ky,
%   and have units of [rad cm^-1].  Values for gravity and surface 
%   tension are specified in GC_PARAMS, and the dispersion relation
%   is given in OM.  
%
%   See also: KFUN, HFUN, RESCOEFF. 
%
%   Usage: D=dfun(i,k1,k2);
%          D=dfun(s1,s2,k1,k2);
%  
%   'dfun --t' runs some tests.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001--2006 J.M. Lilly --- type 'help jlab_license' for details    

if strcmp(varargin{1},'--t')
  dfun_test;return
end

str='notilde'; 
str2='nolimit';
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

a1=abs(k1);
a2=abs(k2);
a0=abs(k1+k2);%Holliday's Appendix B is incorrect about this

%Written in Elfouhaily's notation
%om1=s1.*om(k1);
%om2=s2.*om(k2);
% if strcmp(str,'notilde')  
%   Da=(om1+om2).*(a1.*a2-cdot(k1,k2));
%   Db=om1.*om2.*(om1+om2).*(frac(a2,g+T.*a1.^2)+frac(a1,g+T.*a2.^2));
%   Dc=-(g+T.*a0.^2).*(om2.*frac(a1.^2+cdot(k1,k2),g+T.*a2.^2) ...
%   		  +om1.*frac(a2.^2+cdot(k1,k2),g+T.*a1.^2));
%   D=frac(sqrt(-1),2).*(Da+Db+Dc);
% elseif strcmp(str,'tilde')   
%   Da=(om1+om2).*(a1.*a2-cdot(k1,k2));
%   Db=-2.*(om1+om2).*a1.*a2;  %Note correction here from Elfouhaily (2000) 
%   Dc=-(g+T.*a0.^2).*(om2.*frac(a1.^2-a0.*a1+cdot(k1,k2),g+T.*a2.^2) ...
%   		  +om1.*frac(a2.^2-a0.*a2+cdot(k1,k2),g+T.*a1.^2));
%   D=frac(sqrt(-1),2).*(Da+Db+Dc);
% end

%Written in my notation

if strcmp(str2,'nolimit')
  if strcmp(str,'notilde')  
    Da=(s1.*om(k1)+s2.*om(k2)).*(a1.*a2-cdot(k1,k2));
    Db=s1.*s2.*a1.*a2.*(s1.*om(k1)+s2.*om(k2)).*frac(om(k1).*om(k1)+om(k2).*om(k2),om(k1).*om(k2));
    Dc=-frac(om(k1+k2).*om(k1+k2),a0).*(frac(a2,s2.*om(k2)).*(a1.^2+cdot(k1,k2))+frac(a1,s1.*om(k1)).*(a2.^2+cdot(k1,k2)));
    D=frac(sqrt(-1),2).*(Da+Db+Dc);
  elseif strcmp(str,'tilde')   
    Da=(s1.*om(k1)+s2.*om(k2)).*(a1.*a2-cdot(k1,k2));
    Db=-2.*(s1.*om(k1)+s2.*om(k2)).*a1.*a2;  %Note correction here from Elfouhaily (2000) 
    Dc=-frac(om(k1+k2).*om(k1+k2),a0).*(frac(a2,s2.*om(k2)).*(a1.^2-a0.*a1+cdot(k1,k2))+frac(a1,s1.*om(k1)).*(a2.^2-a0.*a2+cdot(k1,k2)));
    D=frac(sqrt(-1),2).*(Da+Db+Dc);
  end
elseif strcmp(str2,'limit')
  if strcmp(str,'notilde')  
    D=sqrt(-1)*(s1.*om(k1)+s2.*om(k2)).*(a1.*a2-cdot(k1,k2));
  elseif strcmp(str,'tilde')   
    Da=(s1.*om(k1)+s2.*om(k2)).*(a1.*a2-cdot(k1,k2));
    Db=s1.*om(k1).*(a2.^2+3*a1.*a2-a0.*a2)+s2.*om(k2).*(a1.^2+3*a1.*a2-a0.*a1);
    %Note correction here from Elfouhaily's (2000) incorrect version: 
    %Db=s1.*om(k1).*(a2.^2+a1.*a2-a0.*a2)+s2.*om(k2).*(a1.^2+a1.*a2-a0.*a1);
    D=sqrt(-1)*(Da-frac(1,2).*Db);
  end
end


function[]=dfun_test

k1=1/10000;
kg=wavegrid(.126,32)*k1/2;
k=real(kg(1,:))';

for i=1:4
  D1(:,:,i)=dfun(i,k1,kg,'notilde');
  D2(:,:,i)=dfun(i,k1,kg,'notilde','limit');
  D3(:,:,i)=dfun(i,k1,kg,'tilde');
  D4(:,:,i)=dfun(i,k1,kg,'tilde','limit');
end


tol=1e-9;
reporttest('DFUN notilde limit and nolimit versions match at k=1/10000',aresame(D1,D2,tol))
reporttest('DFUN tilde limit and nolimit versions match at k=1/10000',aresame(D3,D4,tol))

if 0
kg=wavegrid(.126,32)*kmin;
k=real(kg(1,:))';

k1=2*kmin+sqrt(-1)*1e-10;

for i=1:4
  D1(:,:,i)=dfun(i,k1,kg,'notilde');
  D2(:,:,i)=dfun(i,k1,kg,'tilde');
end

subplot(2,2,1),pcolor(k./kmin,k./kmin,log10(abs(squeeze(D1(:,:,1))))),shading interp
subplot(2,2,2),pcolor(k./kmin,k./kmin,log10(abs(squeeze(D1(:,:,2))))),shading interp
subplot(2,2,3),pcolor(k./kmin,k./kmin,log10(abs(squeeze(D1(:,:,3))))),shading interp
subplot(2,2,4),pcolor(k./kmin,k./kmin,log10(abs(squeeze(D1(:,:,4))))),shading interp

subplot(2,2,1),pcolor(k./kmin,k./kmin,log10(abs(squeeze(D1(:,:,1))))),shading interp
subplot(2,2,2),pcolor(k./kmin,k./kmin,log10(abs(squeeze(D1(:,:,3))))),shading interp
subplot(2,2,3),pcolor(k./kmin,k./kmin,log10(abs(squeeze(D2(:,:,1))))),shading interp
subplot(2,2,4),pcolor(k./kmin,k./kmin,log10(abs(squeeze(D2(:,:,3))))),shading interp

subplot(2,2,1),pcolor(k./kmin,k./kmin,log10(abs(squeeze(D1(:,:,1))))),shading interp
subplot(2,2,2),pcolor(k./kmin,k./kmin,log10(abs(squeeze(D1(:,:,3))))),shading interp
subplot(2,2,3),pcolor(k./kmin,k./kmin,log10(abs(squeeze(D1(:,:,1)-D2(:,:,1))))),shading interp
subplot(2,2,4),pcolor(k./kmin,k./kmin,log10(abs(squeeze(D1(:,:,3)-D2(:,:,3))))),shading interp
plot(k./kmin,squeeze(D1(:,:,1)-D2(:,:,1)))
plot(k./kmin,squeeze(D1(:,:,3)-D2(:,:,3)))

subplot(2,2,1),pcolor(k./kmin,k./kmin,log10(abs(squeeze(D3(:,:,1))))),shading interp
subplot(2,2,2),pcolor(k./kmin,k./kmin,log10(abs(squeeze(D4(:,:,1))))),shading interp
subplot(2,2,3),pcolor(k./kmin,k./kmin,log10(abs(squeeze(D3(:,:,3))))),shading interp
subplot(2,2,4),pcolor(k./kmin,k./kmin,log10(abs(squeeze(D4(:,:,3))))),shading interp
end
