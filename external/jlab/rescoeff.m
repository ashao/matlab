function[g3,g2,g1]=rescoeff(k1,k2,k3,str)
%RESCOEFF  Interaction coefficients for a resonant wave triad.
%
%   [ALPHA1,ALPHA2,ALPHA3]=RESCOEFF(K1,K2,K3) returns the interaction 
%   coefficients for a resonant wave triad, as presented in Lilly and 
%   Lettvin (2006) and based on McGoldrick (1965) and Simmons (1969).
%   The input wavenumbers should be a resonant triad with K1+K2=K3.
%  
%   RESCOEFF(...,STR) specifies using one of four formulations, which 
%   turn out to be numerically identical.  STR='mcg' [the default] uses 
%   McGoldrick (1965), STR='sim' uses Simmons (1969), STR='hol' uses   
%   Holliday (1977), and STR='elf' uses a corrected version of 
%   Elfouhaily et al. (2000).
%
%   The wavenumber arrays are in complex form, that is, K = Kx + i*Ky,
%   and have units of [rad cm^-1].  Values for gravity and surface 
%   tension are specified in GC_PARAMS, and the dispersion relation is 
%   given in OM.  
%
%   See also: KFUN, DFUN, HFUN, TRIADEVOLVE. 
%
%   Usage: [alpha1,alpha2,alpha3]=rescoeff(k1,k2,k3);
%
%   'rescoeff --t' some tests.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001--2006 J.M. Lilly --- type 'help jlab_license' for details    


if strcmp(k1,'--t')
  rescoeff_test;return
end

if nargin==3
  str='mcg';
end

if strcmp(str(1:3),'mcg')
  [B12,B13,B23]=rescoeff_mcg(k1,k2,k3);
elseif strcmp(str(1:3),'sim')
  [B12,B13,B23]=rescoeff_sim(k1,k2,k3);
elseif strcmp(str(1:3),'hol')
  [B12,B13,B23]=rescoeff_hol(k1,k2,k3);
elseif strcmp(str(1:3),'elf')
  [B12,B13,B23]=rescoeff_elf(k1,k2,k3);
else
  error(['STR = ' str ' is not supported.']) 
end

vswap(B12,nan,0);
vswap(B23,nan,0);
vswap(B13,nan,0);

%All get multiplied by a factor of two since I write 2*cos  
g3=B12*2;
g2=B13*2;
g1=B23*2;

function[B12,B13,B23]=rescoeff_mcg(k1,k2,k3)
  
[g,T]=gc_params;

theta13=angle(k3)-angle(k1);
theta12=angle(k2)-angle(k1);
theta23=angle(k2)-angle(k3);

k1=abs(k1);
k2=abs(k2);
k3=abs(k3);

om1=om(k1);
om2=om(k2);
om3=om(k3);


B12=om1.*om2./(4.*om3.^2).*k3.*...
	((k1.*om3.^2)./(k3.*om2)+(k2.*om3.^2)./(k3.*om1)-om1.^2./om2-om2.^2./om1...
	-4.*om3.*(sin(theta12./2).^2)-T.*cos(theta12).* ...
	(k1.^3./om1+k2.^3./om2-(k3.^2.*k1)./om1-(k3.^2.*k2)./om2));

B13=om1.*om3./(4.*om2.^2).*k2.*...
	(om1.^2./om3-om3.^2./om1+(k3.*om2.^2)./(k2.*om1)-(k1.*om2.^2)./(k2.*om3)...
	+4.*om2.*(cos(theta13./2).^2)-T.*cos(theta13).* ...
	(k3.^3./om3-k1.^3./om1-(k2.^2.*k3)./om3+(k2.^2.*k1)./om1));

B23=om2.*om3./(4.*om1.^2).*k1.*...
	(om2.^2./om3-om3.^2./om2+(k3.*om1.^2)./(k1.*om2)-(k2.*om1.^2)./(k1.*om3)...
	+4.*om1.*(cos(theta23./2).^2)-T.*cos(theta23).* ...
	(k3.^3./om3-k2.^3./om2-(k1.^2.*k3)./om3+(k1.^2.*k2)./om2));

%note the typo in McGoldrick 3.13, B23.  4sigma_3 omhould be 4sigma_1

function[B12,B13,B23]=rescoeff_sim(k1,k2,k3)
  
om1=om(k1);
om2=om(k2);
om3=om(k3);

n1=k1./abs(k1);
n2=k2./abs(k2);
n3=k3./abs(k3);

Ja=om1.*om2.*(1+cdot(n1,n2));
Jb=om2.*om3.*(1-cdot(n2,n3));
Jc=om1.*om3.*(1-cdot(n1,n3));

J=-frac(1,4)*(Ja-Jb-Jc);

B12=-frac(abs(k3),om3).*J;
B13=-frac(abs(k2),om2).*J;
B23=-frac(abs(k1),om1).*J;

function[B12,B13,B23]=rescoeff_hol(k1,k2,k3)
  
B12=kfun(1,1,k1,k2);
B13=kfun(1,-1,k3,k1);
B23=kfun(1,-1,k3,k2);

function[B12,B13,B23]=rescoeff_elf(k1,k2,k3)

%These all differ from the verions of Elfouhaily (2000) by a minus sign
fact1=-frac(abs(k3).*om(k1).*om(k2),2*om(k3).*(om(k1)+om(k2)).*abs(k1).*abs(k2));
fact2=frac(abs(k2).*om(k1).*om(k3),2*om(k2).*(om(k3)-om(k1)).*abs(k3).*abs(k1));
fact3=frac(abs(k1).*om(k2).*om(k3),2*om(k1).*(om(k3)-om(k2)).*abs(k3).*abs(k2));

B12=-sqrt(-1)*fact1.*dfun(1,1,k1,k2);
B13=-sqrt(-1)*fact2.*dfun(1,-1,k3,-k1);
B23=-sqrt(-1)*fact3.*dfun(1,-1,k3,-k2);

%These seem to have trouble at resonace, though they have the right limits 
%B12=kfun(1,1,k1,k2,'notilde','elf');
%B13=kfun(1,-1,k3,-k1,'notilde','elf');
%B23=kfun(1,-1,k3,-k2,'notilde','elf');

function[]=rescoeff_test
  
k3=[16  8   4     2  1.5   1     1/2     1/4 ].*kmin; 
%k2i=(.00000001:.01:10]'*k3;
k2i=(.01:.1:10)'*k3;

clear k1 k2
for i=1:length(k3)
    [k1(:,:,i),k2(:,:,i)]=vtriadres(k3(i),k2i(:,i));
end

k1=squeeze(k1(:,1,:));
k2=squeeze(k2(:,1,:));
k3=osum(zeros(size(k1(:,1))),k3(:))+sqrt(-1)*1e-10;

[B12m,B13m,B23m]=rescoeff(k1,k2,k3,'mcg');
[B12s,B13s,B23s]=rescoeff(k1,k2,k3,'sim');
[B12h,B13h,B23h]=rescoeff(k1,k2,k3,'hol');
[B12e,B13e,B23e]=rescoeff(k1,k2,k3,'elf');

tol=1e-5;
bool(1)=aresame(B12m,B12s,tol);
bool(2)=aresame(B23m,B23s,tol);
bool(3)=aresame(B13m,B13s,tol);
bool(4)=aresame(B12h,B12s,tol);
bool(5)=aresame(B23h,B23s,tol);
bool(6)=aresame(B13h,B13s,tol);
bool(7)=aresame(B12e,B12s,tol);
bool(8)=aresame(B23e,B23s,tol);
bool(9)=aresame(B13e,B13s,tol);

reporttest('RESCOEFF McGoldrick, Simmons, Holliday, and Elouhaily are the same',all(bool));


%S_THSTHD  Script to solve for resonance contitions in |k_1|, |k_2| space
%/********************************************************
%first, I make a wavenumber magnitude grid and solve for angle difference
N=60;
k1=logspace(-3,3,N)'*kmin;

clear th
for i=1:length(k1)
  ths(:,i)=triadres(k1,k1(i))';
end

%test that these satisfy resonance
kg1s=osum(k1,0*k1);
kg2s=osum(0*k1,k1).*rot(ths);
bool=isres(1,1,kg1s,kg2s);
reporttest('RESCOEFF k1s and k2s satisfy (+,+) resonance condition',allall(bool(~isnan(kg2s))))
%looks great


%now for the opposite, difference 
clear thd
for i=1:length(k1)
  temp1=vtriadres(k1(i),k1);
  thd(:,i)=angle(temp1(:,1));
end
thd=thd(1:length(k1),:);
thd=thd+diag(nan*k1);
%something funny on the main diagonal... what is this?

%test that these satisfy resonance
kg1d=osum(k1,0*k1);
kg2d=osum(0*k1,k1).*rot(thd);
index=find(~isnan(kg2d));
bool=isres(-1,1,kg1d,kg2d);
reporttest('RESCOEFF k1d and k2d satisfy (-,+) resonance condition',allall(bool(index)))

%swap k1 and k2
temp=kg2d;kg2d=kg1d;kg1d=temp;clear temp
bool=isres(1,-1,kg1d,kg2d);
reporttest('RESCOEFF swapped k1d and k2d satisfy (+,-) resonance condition',allall(bool(index)))
%output: k1,ths,thd, kg1s, kg2s, kg1d, kg2d
%\********************************************************


%/********************************************************
mi=nan*ones(length(thd(:,1)),3);
k2m=nan*ones(length(thd(:,1)),3);
for i=1:size(thd,2)
  temp=min(find(~isnan(ths(:,i))));
  if ~isempty(temp)
    mi(i,1)=temp;
    k2m(i,1)=k1(temp);
  end
  
  temp=max(find(~isnan(thd(:,i))));
  if ~isempty(temp)
    mi(i,2)=temp;
    k2m(i,2)=k1(temp);
  end
  
  temp=min(find(~isnan(thd(:,i))));
  if ~isempty(temp)
    mi(i,3)=temp;
    k2m(i,3)=k1(temp);
  end
end

k1m(:,1)=k1;
k1m(:,2)=k1;
k1m(1:end/2,2)=flipud(k1(end/2+1:end));
k2m(1:end/2,2)=flipud(k2m(end/2+1:end,1));
k2m(31,2)=k2m(30,2);  %This just repeats earlier value for plotting
%********************************************************
%output: k1,ths,thd, kg1s, kg2s, kg1d, kg2d

%/********************************************************
kg3s=kg1s+kg2s;

n=1.99;
KS=kfun(1,1,kg1s,kg2s)./(abs(kg1s).*abs(kg2s)).^n;
KD1=kfun(1,-1,kg3s,-kg1s)./(abs(kg3s).*abs(kg1s)).^n;
KD2=kfun(1,-1,kg3s,-kg2s)./(abs(kg3s).*abs(kg2s)).^n;
vswap(KS,KD1,KD2,nan+sqrt(-1)*nan,0);% imaginary nans result when real nans are rooted
vswap(KS,KD1,KD2,nan,0);
K1=KS-KD1-KD2;
reporttest('RESCOEFF net triad asymmertry is not everywhere positive for k^-1.99',minmin(K1)<=0)

n=2;
KS=kfun(1,1,kg1s,kg2s)./(abs(kg1s).*abs(kg2s)).^n;
KD1=kfun(1,-1,kg3s,-kg1s)./(abs(kg3s).*abs(kg1s)).^n;
KD2=kfun(1,-1,kg3s,-kg2s)./(abs(kg3s).*abs(kg2s)).^n;
vswap(KS,KD1,KD2,nan+sqrt(-1)*nan,0);% imaginary nans result when real nans are rooted
vswap(KS,KD1,KD2,nan,0);
K1=KS-KD1-KD2;
reporttest('RESCOEFF net triad asymmertry is everywhere positive for k^-2',maxmax(K1)>=0)

n=1.249;
KS=kfun(1,1,kg1s,kg2s)./(abs(kg1s).*abs(kg2s)).^n;
KD1=kfun(1,-1,kg3s,-kg1s)./(abs(kg3s).*abs(kg1s)).^n;
KD2=kfun(1,-1,kg3s,-kg2s)./(abs(kg3s).*abs(kg2s)).^n;
vswap(KS,KD1,KD2,nan+sqrt(-1)*nan,0);% imaginary nans result when real nans are rooted
vswap(KS,KD1,KD2,nan,0);
K1=KS.^2-KD1.^2-KD2.^2;% imaginary nans result when real nans are rooted
reporttest('RESCOEFF net triad energy change is not everywhere positive for k^-1.249',minmin(K1)<0)

n=1.25;
KS=kfun(1,1,kg1s,kg2s)./(abs(kg1s).*abs(kg2s)).^n;
KD1=kfun(1,-1,kg3s,kg1s)./(abs(kg3s).*abs(kg1s)).^n;
KD2=kfun(1,-1,kg3s,kg2s)./(abs(kg3s).*abs(kg2s)).^n;
vswap(KS,KD1,KD2,nan+sqrt(-1)*nan,0);% imaginary nans result when real nans are rooted
vswap(KS,KD1,KD2,nan,0);
K1=KS.^2-KD1.^2-KD2.^2;% imaginary nans result when real nans are rooted
reporttest('RESCOEFF net triad energy change is everywhere positive for k^-1.25',minmin(K1)>=0)


