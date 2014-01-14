function[K]=kfun(varargin)
% KFUN  Combined interaction functions for wave triads.
%
%   KFUN(I,K1,K2) returns the i-th K-function operating on the 
%   wavenumber arrays K1 and K2 as defined Lilly and Lettvin (2006)
%   and based on Holliday (1977). 
%
%   KFUN(S1,S2,K1,K2) also works, where S1 and S2 are each either 
%   positive or negative one. See SS2I for the relation between these
%   two conventions.  
%
%   The K-functions have been made to be symmetric, and also differ 
%   from Holliday's original versions due to a change of variables
%   in the summation of the forcing equation, as discussed in Lilly 
%   and Lettvin (2006).
%
%   The K-functions come in two variations.  KFUN(...,'tilde') returns
%   the interaction functions resulting from expanding the potential 
%   about the sea surface.  This is the default if no string is input.
%   KFUN(...,'notilde') returns the interaction functions resulting 
%   from expanding the potential about zero.
%  
%   The wavenumber arrays are in complex form, that is, K = Kx + i*Ky,
%   and have units of [rad cm^-1].  Values for gravity and surface 
%   tension are specified in GC_PARAMS, and the dispersion relation
%   is given in OM.  
%
%   See also: DFUN, HFUN, RESCOEFF. 
%
%   Usage: K=kfun(i,k1,k2);
%          K=kfun(s1,s2,k1,k2);
%  
%   'KFUN  --t' runs some tests.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001--2006 J.M. Lilly --- type 'help jlab_license' for details    


if strcmp(varargin{1},'--t')
  kfun_test;return
end

str='tilde'; 
str2='sim';
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

if strcmp(str2(1:3),'hol')
  K=kfun_hol(s1,s2,k1,k2,str);
elseif strcmp(str2(1:3),'elf')
  K=kfun_elf(s1,s2,k1,k2,str);
elseif strcmp(str2(1:3),'sim')
  K=kfun_sim(s1,s2,k1,k2,str);
elseif strcmp(str2(1:3),'lil')
  K=kfun_lil(s1,s2,k1,k2,str);
end

function[K]=kfun_hol(s1,s2,k1,k2,str) 
%Implements symmetric versions Holliday's (1977) K-functions
K=frac(1,2)*kfun_hol_asym(s1,s2,s1.*k1,s2.*k2,str)+frac(1,2)*kfun_hol_asym(s2,s1,s2.*k2,s1.*k1,str); 


function[K]=kfun_hol_asym(s1,s2,k1,k2,str) 
%Implements asymmetric versions Holliday's (1977) K-functions
 
warning('off','MATLAB:divideByZero')
K = frac(s1,2).*cfun(1,k1,k2,str)-frac(1,2).*cfun(2,k1,k2,str)+frac(s1.*s2,2).*cfun(3,k1,k2,str);
warning('on','MATLAB:divideByZero')

function[x]=cfun(jj,k1,k2,str) 
%Implements Holliday's (1977) C- and C~-functions

warning('off','MATLAB:divideByZero')
if strcmp(str,'tilde')    
  switch jj
         case 1
           x=frac(om(k1),abs(k1)).*...
  		(cdot(k1,k2)+abs(k1).*(abs(k1)-abs(k1+k2)));
         case 2
           x=0;
         case 3
           x=frac(abs(k1+k2),om(k1+k2)).*...
  		frac(om(k1).*om(k2),abs(k1).*abs(k2)).*...
  		1/2.*(cdot(k1,k2)+abs(k1).*abs(k2));
  end
elseif  strcmp(str,'notilde')   
  switch jj
         case 1
           x=frac(om(k1),abs(k1)).*(cdot(k1,k2)+abs(k1).^2);
         case 2
           x=frac(abs(k1+k2),om(k1+k2)).*om(k1).*om(k1);
         case 3
           x=frac(abs(k1+k2),om(k1+k2)).*...
  		frac(om(k1).*om(k2),abs(k1).*abs(k2)).*...
  		1/2.*(cdot(k1,k2)-abs(k1).*abs(k2));
  end
end
warning('on','MATLAB:divideByZero')

x=vswap(x,nan,0);

warning('off','MATLAB:divideByZero')


function[K]=kfun_lil(s1,s2,k1,k2,str)
%K-function versions with Lilly simplification of Holliday
k3=s1.*k1+s2.*k2;

om1=om(k1);
om2=om(k2);
om3=om(k3);

Ka=frac(1,4)*frac(om1,abs(k1)).*(s1.*abs(k1).*(abs(k1)-abs(k3))+s2.*cdot(k1,k2));
Kb=frac(1,4)*frac(om2,abs(k2)).*(s2.*abs(k2).*(abs(k2)-abs(k3))+s1.*cdot(k1,k2));
Kc=frac(1,4)*frac(abs(k3),om3).*frac(om1.*om2,abs(k1).*abs(k2)).*(s1.*s2.*abs(k1).*abs(k2)+cdot(k1,k2));
K=Ka+Kb+Kc;

if strcmp(str,'notilde')   
    K=K+frac(abs(k3),4).*frac(s1.*om1+s2.*om2,om3).*(om3-s1.*om1-s2.*om2);
end


function[K]=kfun_sim(s1,s2,k1,k2,str)
%K-function versions based on Simmons
k3=s1.*k1+s2.*k2;

k1hat=k1./abs(k1);
k2hat=k2./abs(k2);
k3hat=k3./abs(k3);

om1=om(k1);
om2=om(k2);
om3=om(k3);

Ka=frac(1,4)*frac(abs(k3),om3).*om1.*om2.*(s1.*s2+cdot(k1hat,k2hat));
Kb=frac(1,4)*frac(abs(k3),om3).*om1.*om3.*(s1    -cdot(k1hat,k3hat));
Kc=frac(1,4)*frac(abs(k3),om3).*om2.*om3.*(s2    -cdot(k2hat,k3hat));
K=Ka-Kb-Kc;

if strcmp(str,'notilde')   
    K=K+frac(abs(k3),4).*frac(s1.*om1+s2.*om2,om3).*(om3-s1.*om1-s2.*om2);
end


function[K]=kfun_elf(s1,s2,k1,k2,str)
%K-function versions based on Elfouhaily  
Ka=-(om(k1+k2)-s1.*om(k1)-s2.*om(k2)).*gfun(s1,s2,k1,k2,str,'elf');
Kb=-frac(abs(k1+k2),om(k1+k2)).*gfun(s1,s2,k1,k2,str,'elf')./(om(k1+k2)+s1.*om(k1)+s2.*om(k2));
K=frac(1,2).*(Ka+Kb);      

function[]=kfun_test
kfun_test_simmons
kfun_test_lilly
kfun_test_tilde_vs_notilde
%kfun_test_elfouhaily

function[]=kfun_test_tilde_vs_notilde

k1=2*kmin+sqrt(-1)*1e-10;
kg=wavegrid(.126,32)*kmin;
k=real(kg(1,:))';

om1=om(k1);
om2=om(kg);
  
clear K1 K2 K3 K4
for i=1:4
 [s1,s2]=i2ss(i); 
 k3=s1.*k1+s2.*kg;
 om3=om(k3);
 
 K_tilde=kfun(i,k1,kg,'tilde','sim');
 K1(:,:,i)=K_tilde +  frac(abs(k3),4).*frac(s1.*om1+s2.*om2,om3).*(om3-s1.*om1-s2.*om2);
 K2(:,:,i)=kfun(i,k1,kg,'notilde','sim');
end

tol=1e-6;
reporttest('KFUN Hol tilde minus notilde matches prediction at 2*kmin off-resonance',aresame(K1,K2,tol))


k1=4*kmin+sqrt(-1)*1e-10;
kg=wavegrid(.126,32)*kmin;
k=real(kg(1,:))';

om1=om(k1);
om2=om(kg);
  
clear K1 K2 K3 K4
for i=1:4
 [s1,s2]=i2ss(i); 
 k3=s1.*k1+s2.*kg;
 om3=om(k3);
 
 K_tilde=kfun(i,k1,kg,'tilde','sim');
 K1(:,:,i)=K_tilde +  frac(abs(k3),4).*frac(s1.*om1+s2.*om2,om3).*(om3-s1.*om1-s2.*om2);
 K2(:,:,i)=kfun(i,k1,kg,'notilde','sim');
end

tol=1e-6;
reporttest('KFUN Hol tilde minus notilde matches prediction at 4*kmin off-resonance',aresame(K1,K2,tol))

function[]=kfun_test_simmons
  
kg=wavegrid(.126,32)*kmin;
k=real(kg(1,:))';

k1=2*kmin+sqrt(-1)*1e-10;
  
clear K1 K2 K3 K4
for i=1:4
 [s1,s2]=i2ss(i);
  
 K1(:,:,i)=kfun(i,k1,kg,'tilde','sim');
 K2(:,:,i)=kfun(i,k1,kg,'tilde','hol');
 K3(:,:,i)=kfun(i,k1,kg,'notilde','sim');
 K4(:,:,i)=kfun(i,k1,kg,'notilde','hol');
end

tol=1e-6;
reporttest('KFUN Holiday and Simmons tilde versions match at 2*kmin off-resonance',aresame(K1,K2,tol))
reporttest('KFUN Holiday and Simmons no-tilde versions match at 2*kmin off-resonance',aresame(K3,K4,tol))

k1=4*kmin+sqrt(-1)*1e-10;
  
for i=1:4
 [s1,s2]=i2ss(i);
 K1(:,:,i)=kfun(i,k1,kg,'tilde','sim');
 K2(:,:,i)=kfun(i,k1,kg,'tilde','hol');
 K3(:,:,i)=kfun(i,k1,kg,'notilde','sim');
 K4(:,:,i)=kfun(i,k1,kg,'notilde','hol');
end

tol=1e-5;
reporttest('KFUN Holiday and Simmons tilde versions match at 4*kmin off-resonance',aresame(K1,K2,tol))
reporttest('KFUN Holiday and Simmons no-tilde versions match at 4*kmin off-resonance',aresame(K3,K4,tol))

function[]=kfun_test_lilly
  
kg=wavegrid(.126,32)*kmin;
k=real(kg(1,:))';

k1=2*kmin+sqrt(-1)*1e-10;
  
clear K1 K2 K3 K4
for i=1:4
 [s1,s2]=i2ss(i);
  
 K1(:,:,i)=kfun(i,k1,kg,'tilde','lil'); 
 K2(:,:,i)=kfun(i,k1,kg,'tilde','hol');
 K3(:,:,i)=kfun(i,k1,kg,'notilde','lil'); 
 K4(:,:,i)=kfun(i,k1,kg,'notilde','hol');
end

tol=1e-6;
reporttest('KFUN Holiday and Lilly tilde versions match at 2*kmin off-resonance',aresame(K1,K2,tol))
reporttest('KFUN Holiday and Lilly no-tilde versions match at 2*kmin off-resonance',aresame(K3,K4,tol))

k1=4*kmin+sqrt(-1)*1e-10;
  
for i=1:4
 [s1,s2]=i2ss(i);
 K1(:,:,i)=kfun(i,k1,kg,'tilde','lil');
 K2(:,:,i)=kfun(i,k1,kg,'tilde','hol');
 K3(:,:,i)=kfun(i,k1,kg,'notilde','lil');
 K4(:,:,i)=kfun(i,k1,kg,'notilde','hol');
end

tol=1e-5;
reporttest('KFUN Holiday and Lilly tilde versions match at 4*kmin off-resonance',aresame(K1,K2,tol))
reporttest('KFUN Holiday and Lilly no-tilde versions match at 4*kmin off-resonance',aresame(K3,K4,tol))

function[]=kfun_test_elfouhaily
  
kg=wavegrid(.126,32)*kmin;
k=real(kg(1,:))';

k1=2*kmin+sqrt(-1)*1e-10;
  
clear K1 K2 K3 K4
for i=1:4
 K1(:,:,i)=kfun(i,k1,kg,'tilde','elf');
 K2(:,:,i)=kfun(i,k1,kg,'tilde','hol');
 K3(:,:,i)=kfun(i,k1,kg,'notilde','elf');
 K4(:,:,i)=kfun(i,k1,kg,'notilde','hol');
end

tol=1e-6;
reporttest('KFUN Hol and Elf tilde versions match at 2*kmin off-resonance',aresame(K1,K2,tol))
reporttest('KFUN Hol and Elf no-tilde versions match at 2*kmin off-resonance',aresame(K3,K4,tol))

k1=4*kmin+sqrt(-1)*1e-10;
  
for i=1:4
 K1(:,:,i)=kfun(i,k1,kg,'tilde','elf');
 K2(:,:,i)=kfun(i,k1,kg,'tilde','hol');
 K3(:,:,i)=kfun(i,k1,kg,'notilde','elf');
 K4(:,:,i)=kfun(i,k1,kg,'notilde','hol');
end

tol=1e-5;
reporttest('KFUN Hol and Elf tilde versions match at 4*kmin off-resonance',aresame(K1,K2,tol))
reporttest('KFUN Hol and Elf no-tilde versions match at 4*kmin off-resonance',aresame(K3,K4,tol))
