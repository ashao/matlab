function[varargout]=morsewave(varargin)
% MORSEWAVE  Generalized Morse wavelets of Olhede and Walden (2002). 
%
%   PSI=MORSEWAVE(N,K,GAMMA,BETA,F) returns an N x K column vector PSI 
%   which contains time-domain versions of the first K generalized Morse
%   wavelets specified by GAMMA and BETA, concentrated at frequency F.
%
%   The frequency F is specifically the *radian* frequency at which the 
%   Fourier transform of the lowest-order (K=1) wavelet has its maximum
%   amplitude, assuming a unit sample rate.  If F has length L, PSI is of
%   size N x L x K, with the columns in order of decreasing frequency.
%
%   Note that the wavelets are centered at the midpoint in time, that is, 
%   row number ROUND(SIZE(PSI,1)/2).
%   _________________________________________________________________
%
%   Normalization
%
%   MORSEWAVE supports two kinds of normalization for the wavelets.
%
%   MORSEWAVE(..., 'bandpass') uses "bandpass normalization", the default
%   See Lilly and Gascard (2006).  This implies that the FFT of the lowest-
%   order (K=1) generalized Morse wavelet has a maximum value of 2. 
%
%   MORSEWAVE(..., 'energy') uses the unit energy normlization.  The time-
%   domain wavelet energy SUM(ABS(PSI).^2) is unity for all K. 
%
%   Note MORSEWAVE now uses bandpass-normalization by default.
%   _________________________________________________________________
%
%   Time and frequency versions
%
%   [PSI,PSIF]=MORSEWAVE(...) with two output arguments optionally returns
%   a frequency-domain version PSIF of the wavelets.  
%
%   PSIF is the same size as PSI.
%
%   PSIF=MORSEWAVE(...,'frequency') supresses the computation of the time-
%   domain wavelets.  This is more efficient if only PSIF is desired since
%   PSIF is actually computed first.
%
%   MORSEWAVE can be called with multiple string arguments.  These can be 
%   in any order after the five required arguments. 
%   _________________________________________________________________
%
%   Background
%
%   For further details on generalized Morse wavelets, see
%
%       Olhede and Walden (2002),  Generalized Morse Wavelets", IEEE 
%           Trans. Sig. Proc., 50 (11), 2661--2670.
%
%       Lilly and Olhede (2009).  Higher-order properties of analytic 
%           wavelets.  IEEE Trans. Sig. Proc., 57 (1), 146--160.
%   _________________________________________________________________
%
%   Usage: psi=morsewave(N,K,ga,be,f);
%          [psi,psif]=morsewave(N,K,ga,be,f);
%          [psi,psif]=morsewave(N,K,ga,be,f,'energy');
%
%   'morsewave --t' runs a test
%   'morsewave --f' generates a sample figure
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2011 J.M. Lilly and F. Rekibi
%                         --- type 'help jlab_license' for details  
 


%   _________________________________________________________________
%
%   The zero beta case
%
%   It is unlikely that you will need to read this section, which is 
%   why the comment is hidden.  It describes a feature mostly used for 
%   testing purposes.
%
%   For BETA equal to zero, the generalized Morse wavelets describe
%   a non-zero-mean function which is not in fact a wavelet. 
%
%   Only 'bandpass' normalization is supported for this case.
%
%   In this case the frequency speficies the half-power point of the
%   analytic lowpass filter.  
%
%   The frequency-domain definition of MOReSEWAVE is not necessarily 
%   a good way to compute the zero-beta functions, however.  You will
%   probably need to take a very small DT.
%   _________________________________________________________________


%   This behavior is depricated  
%
%   Sample rate
%
%   MORSEWAVE(DT,N,K,GAMMA,BETA,F) specifies the sample rate to be DT.
%   DT is optional with a default value of unity.
%   _________________________________________________________________


if strcmp(varargin{1},'--t')
  morsewave_test;return
end
if strcmp(varargin{1},'--f')
  morsewave_fig;return
end

str='bandpass';
fam='first';
dom='tim';

for i=1:3
if ischar(varargin{end})
    temp=varargin{end};
    if strcmp(temp(1:3),'fir')||strcmp(temp(1:3),'sec')
       fam=temp;
    elseif strcmp(temp(1:3),'ban')||strcmp(temp(1:3),'ene')
       str=temp;
    elseif strcmp(temp(1:3),'tim')||strcmp(temp(1:3),'fre')
       dom=temp;
    end
    varargin=varargin(1:end-1);
end
end

dt=1;
N=varargin{1};
K=varargin{2};
ga=varargin{3};
be=varargin{4};
fs=varargin{5};


str=str(1:3);

if be==0&&strcmp(str,'ene')
    str='ban';
    disp('For BETA=0, energy normalization is not defined.  Using bandpass normalization.')
end

for n=1:length(fs)
    [X(:,:,n),x(:,:,n)]=morsewave1(N,K,ga,be,abs(fs(n)),str,fam,dom);
    if fs(n)<0
        if ~isempty(x)
            x(:,:,n)=conj(x(:,:,n));
        end
        X(2:end,:,n)=flipdim(X(2:end,:,n),1);
    end
end

if size(X,3)>1
  if ~isempty(x)
    x=permute(x,[1 3 2]);
  end
  if nargout==2||~strcmp(dom(1:3),'tim')
    X=permute(X,[1 3 2]);
  end
end

if ~strcmp(dom(1:3),'tim')
    varargout{1}=X;
else
    varargout{1}=x;
    varargout{2}=X;
end

function[X,x]=morsewave1(N,K,ga,be,fs,str,fam,dom)
  
fo=morsefreq(ga,be);
fact=fs./fo;
om=2*pi*linspace(0,1-1./N,N)'./fact;
psizero=double((om.^be).*exp(-om.^ga));

psizero(1)=1/2*psizero(1); %Due to unit step function
%Ensure nice lowpass filters for beta=0;
%Otherwise, doesn't matter since wavelets vanishes at zero frequency


vswap(psizero,nan,0);

if strcmp(fam(1:3),'fir')
   X=morsewave_first_family(fact,N,K,ga,be,om,psizero,str); 
elseif strcmp(fam(1:3),'sec')
   X= morsewave_second_family(fact,N,K,ga,be,om,psizero,str); 
end

X=vswap(X,inf,0);

Xr=X.*rot(vrep(om,size(X,2),2).*(N+1)/2*fact); %ensures wavelets are centered 


if strcmp(dom(1:3),'tim')
    x=ifft(Xr);
else
    x=[];
end

function[X]=morsewave_first_family(fact,N,K,ga,be,om,psizero,str)

r=(2*be+1)./ga;
c=r-1;
L=0*om;
index=(1:round(N/2));

for k=0:K-1
  %Log of gamma function much better ... trick from Maltab's ``beta''
  if strcmp(str,'ene')
        %A=(pi*ga*(2.^r)*gamma(k+1)/gamma(k+r)).^(1/2) 
        A=double((pi*ga*(2.^r)*exp(gammaln(k+1)-gammaln(k+r))).^(1/2));
        coeff = sqrt(2./fact)*A;
  elseif strcmp(str,'ban')
        %A0=(pi*ga*(2.^r)/gamma(r)).^(1/2)
        %coeff= morsea(ga,be).*frac(A,A0)
        coeff=double(morsea(ga,be).*sqrt(exp(gammaln(r)+gammaln(k+1)-gammaln(k+r))));
        if be==0
            coeff=2;
        end
  end
    L(index)=laguerre(2*om(index).^ga,k,c);
    X(:,k+1)=coeff.*psizero.*L;%maxmax(coeff),maxmax(psizero),maxmax(L)
end


%  See Olhede and Walden, "Noise reduction in directional signals
%  using multiple Morse wavelets", IEEE Trans. Bio. Eng., v50, 51--57.
%  The equation at the top right of page 56 is equivalent to the
%  preceding expressions. Morse wavelets are defined in the frequency  
%  domain, and so not interpolated in the time domain in the same way
%  as other continuous wavelets.


function[X]=morsewave_second_family(fact,N,K,ga,be,om,psizero,str)

dt=1;
dom=om(2)-om(1);
a0=morsea(ga,be,'energy');
index=(1:round(N/2));
psi0=dt.*sqrt(1./fact).*a0.*psizero;

if K>3
    error('Sorry, can only compute the first 3 members of this family right now.');
end

phi=zeros(length(om),K);
    
for k=0:K-1
    ak=morsewave_ak(k,ga,be);
    for n=0:k
        cnk=morsewave_cnk(n,k,ga,be);
        phi(index,k+1)=phi(index,k+1)+ak.*cnk.*om(index).^(n.*ga-k).*psi0(index);
    end
end

%Ensure zero mean
if iseven(N)
    phi(1,:)=0;
end

psi=phi;
for k=0:K-1
    for n=0:k-1
        bnk=morsewave_bkl(k,n,ga,be);
        psi(:,k+1)=psi(:,k+1)-bnk.*phi(:,n+1);
    end
%   morsewave_atildek(k,ga,be)
    %psi(:,k+1)=morsewave_atildek(k,ga,be).*psi(:,k+1);
    psi(:,k+1)=psi(:,k+1)./sqrt(vsum(psi(:,k+1).^2,1))*sqrt(N);
end
X=psi;

%X=phi;


function[cnk]=morsewave_cnk(n,k,ga,be)
if k==0&&n==0
    cnk=1;
elseif k==1&&n==0||k==0&&n==1
     cnk=be;
elseif k==1&&n==1
     cnk=-ga;        
elseif k==2&&n==0||k==0&&n==2
     cnk=be.*(be-1);
elseif k==2&&n==1||k==1&&n==2
     cnk=-(ga.*(ga-1)+2.*be.*ga); 
elseif k==2&&n==2
     cnk=ga.^2;                     
end


function[ak]=morsewave_ak(k,ga,be)
akinv=0;

for n=0:k
    for p=0:k
        cnk=morsewave_cnk(n,k,ga,be);
        cpk=morsewave_cnk(p,k,ga,be);
        ratn=frac(morsea(ga,be,'energy'),morsea(ga,be+n*ga-k,'energy'));
        ratp=frac(morsea(ga,be,'energy'),morsea(ga,be+p*ga-k,'energy'));
        akinv=akinv+cnk.*cpk.*morseproj(ga,be+n*ga-k,be+p*ga-k).*ratn.*ratp;
    end
end
ak=sqrt(1./akinv);


function[bkl]=morsewave_bkl(k,l,ga,be)
bkl=0;
for n=0:k
    for p=0:l
        cnk=morsewave_cnk(n,k,ga,be);
        cpk=morsewave_cnk(p,l,ga,be);
        ratn=frac(cnk,morsea(ga,be+n*ga-k,'energy'));
        ratp=frac(cpk,morsea(ga,be+p*ga-l,'energy'));
        bkl=bkl+morseproj(ga,be+n*ga-k,be+p*ga-l).*ratn.*ratp;
    end
end
bkl=bkl.*morsea(ga,be,'energy').^2.*morsewave_ak(k,ga,be).*morsewave_ak(l,ga,be);
    

function[dnk]=morsewave_dnk(n,k,ga,be)

a1=morsewave_ak(1,ga,be);
a2=morsewave_ak(2,ga,be);
if k==0&&n==0
    dnk=1;
elseif k==1&&n==0
     dnk=-a1.*morsewave_cnk(0,1,ga,be);
elseif k==1&&n==1
     dnk=1;        
elseif k==2&&n==0
     dnk=(a2.*morsewave_cnk(0,2,ga,be).*a1-1).*morsewave_cnk(0,1,ga,be);
elseif k==2&&n==1
     dnk=-a2.*morsewave_cnk(1,2,ga,be).*a1;
elseif k==2&&n==2
     dnk=1;                     
end



%This is not yet working, and also, is really slow.
function[atildek]=morsewave_atildek(k,ga,be)

akinv=0;

for n=0:k
    for p=0:k
        dnk=morsewave_dnk(n,k,ga,be);
        dpk=morsewave_dnk(p,k,ga,be);
        %ratn=frac(morsea(ga,be,'energy'),morsea(ga,be+n*ga-k,'energy'));
        %ratp=frac(morsea(ga,be,'energy'),morsea(ga,be+p*ga-k,'energy'));
        akinv=akinv+dnk.*dpk.*morsewave_bkl(n,p,ga,be);%.*ratn.*ratp;
    end
end
atildek=sqrt(1./akinv);



function[]=morsewave_test
morsewave_test_unitenergy
morsewave_test_centering
morsewave_test_scorer
morsewave_test_cauchy
morsewave_test_gaussian
morsewave_test_dawson
morsewave_test_admiss
morsewave_test_fft
morsewave_test_notime

function[]=morsewave_test_unitenergy
fs=2*pi./logspace(log10(5),log10(40))'; 
N=1023;
w=morsewave(N,2,2,4,fs,'energy');
energy=vsum(abs(w(:,:,1)).^2,1);
reporttest('MORSEWAVE unit energy for unit sample rate, K=1',maxmax(energy-1)<1e-4)
energy=vsum(abs(w(:,:,2)).^2,1);
reporttest('MORSEWAVE unit energy for unit sample rate, K=2',maxmax(energy-1)<1e-4)


% dt=0.1;
% w=morsewave(dt,N,2,2,4,fs./dt,'energy');
% energy=vsum(abs(w(:,:,1)).^2,1);
% reporttest('MORSEWAVE unit energy for non-unit sample rate, K=1',maxmax(energy-1)<1e-4)
% energy=vsum(abs(w(:,:,2)).^2,1);
% reporttest('MORSEWAVE unit energy for non-unit sample rate, K=2',maxmax(energy-1)<1e-4)


function[]=morsewave_test_centering
fs=2*pi./logspace(log10(5),log10(40))'; 
N=1023;
w=morsewave(N,1,2,4,fs);
bool=0*fs;
for i=1:size(w,2)
   bool(i)=max(abs(w(:,i)))==abs(w(N/2+1/2,i));
end
reporttest('MORSEWAVE centered for odd N',all(bool))

N=1024;
w=morsewave(N,1,2,4,fs);
bool=0*fs;
for i=1:size(w,2)
   bool(i)=max(abs(w(:,i)))==abs(w(N/2,i)) || max(abs(w(:,i)))==abs(w(N/2+1,i));
end
reporttest('MORSEWAVE centered for even N',all(bool))


function[]=morsewave_test_scorer

dt=1;
t=(-50:dt:50)';
psi1=morsewave(length(t),1,3,0,morsefreq(3,1)./dt,'bandpass');
c=3.^(1/3);
psi2=(1./c).*scorer(sqrt(-1)*t./c,1000);
%figure,uvplot(psi1),
%figure,uvplot(psi2)

err=vsum(abs(psi1-psi2).^2,1)./vsum(abs(psi1).^2,1);
reporttest('MORSEWAVE for GAMMA=3 wavelet versus Scorer function expression',err<1e-1)

function[hi]=scorer(z,n)
%This is not a very good way to compute the Scorer functions.
%Need to set n really high to have to integal behave nicely.
%Only used for testing purposes.

z=z(:);
u=linspace(0,10,n)';
du=u(2)-u(1);
umat=osum(0*z,u);

%aiz=airy(0,z);
%biz=airy(2,z);
%gi=frac(1,pi).*vsum(sin(frac(tmat.^3,3)+oprod(z,t)),2).*dt;
hi=frac(1,pi).*vsum(exp(-frac(umat.^3,3)+oprod(z,u)),2).*du;


function[]=morsewave_test_cauchy


dt=.01;
t=(-25:dt:25)';
psi1=morsewave(length(t),1,1,0,morsefreq(1,1).*dt,'bandpass')./dt;
psi2=frac(1,pi).*frac(1,1-sqrt(-1)*t);

err=vsum(abs(psi1-psi2).^2,1)./vsum(abs(psi1).^2,1);
reporttest('MORSEWAVE for GAMMA=1 wavelet versus Cauchy function expression',err<1e-1)


function[]=morsewave_test_gaussian

dt=.1;
t=(-50:dt:50)';
psi1=morsewave(length(t),1,2,0,morsefreq(2,1).*dt,'bandpass')./dt;
psi2=frac(1,2*sqrt(pi)).*(exp(-(t/2).^2)+sqrt(-1)*dawson(t/2)*frac(2,sqrt(pi)));
err=vsum(abs(psi1-psi2).^2,1)./vsum(abs(psi1).^2,1);
reporttest('MORSEWAVE for GAMMA=2 wavelet versus Gaussian function expression',err<1e-1)

function[]=morsewave_test_dawson
dt=0.01;
t=(-15:dt:15)';

n=5;
herm=hermpoly(t(:)/2,n+1);
herm=herm(:,2:end);
g=exp(-frac(t.^2,4));

[psi1,psi2]=vzeros(length(t),5);
for k=1:5
    dk=dawsonderiv(t/2,k);
    coeffk=frac(1,4*sqrt(pi)).*morsea(2,k).*frac(sqrt(-1),2).^k;
    tic
    psi1(:,k)=morsewave(length(t),1,2,k,morsefreq(2,k).*dt,'bandpass')./dt;
    toc
    tic
    psi2(:,k)=coeffk*(g.*herm(:,k)+sqrt(-1)*(-1).^k.*dk*frac(2,sqrt(pi)));
    toc
    err=vsum(abs(psi1(:,k)-psi2(:,k)).^2,1)./vsum(abs(psi1(:,k)).^2,1);
    reporttest(['MORSEWAVE for GAMMA=2 derivatives matches Dawson expression for n=' int2str(k)],err<1e-3)
end


function[]=morsewave_test_admiss
ga1=(1:1:11);
be1=(1:1:10);
[ga,be]=meshgrid(ga1,be1);
vcolon(ga,be);

N=1000;
psi2om=zeros(N,length(ga));
dt=1/10;
om=(0:N-1)'./N;om(1)=1e-4;
for i=1:length(ga)
    [psii,psifi]=morsewave(N,1,ga(i),be(i),morsefreq(ga(i),be(i)).*dt,'bandpass');
    psi2om(:,i)=(psifi).^2./om;
end
cpsi1=vsum(psi2om,1).*(1/N);
cpsi1=reshape(cpsi1,length(be1),length(ga1));
cpsi2=morsea(ga,be).^2.*frac(1,ga.*2.^(2*be./ga)).*gamma(frac(2*be,ga));
cpsi2=reshape(cpsi2,length(be1),length(ga1));
        
%ga1=(1:1:11);
be1=(1:1:10);
[ga,be]=meshgrid(ga1,be1);
vcolon(ga,be);

N=1000;
psi2om=zeros(N,length(ga));
dt=1/10;
om=(0:N-1)'./N;om(1)=1e-4;
for i=1:length(ga)
    [psii,psifi]=morsewave(N,1,ga(i),be(i),morsefreq(ga(i),be(i)).*dt,'bandpass');
    psi2om(:,i)=(psifi).^2./om;
end
cpsi1=vsum(psi2om,1).*(1/N);
cpsi1=reshape(cpsi1,length(be1),length(ga1));
cpsi2=morsea(ga,be).^2.*frac(1,ga.*2.^(2*be./ga)).*gamma(frac(2*be,ga));
cpsi2=reshape(cpsi2,length(be1),length(ga1));
ga1=(1:1:11);
be1=(1:1:10);
[ga,be]=meshgrid(ga1,be1);
vcolon(ga,be);

N=1000;
psi2om=zeros(N,length(ga));
dt=1/10;
om=(0:N-1)'./N;om(1)=1e-4;
for i=1:length(ga)
    [psii,psifi]=morsewave(N,1,ga(i),be(i),morsefreq(ga(i),be(i)).*dt,'bandpass');
    psi2om(:,i)=(psifi).^2./om;
end
cpsi1=vsum(psi2om,1).*(1/N);
cpsi1=reshape(cpsi1,length(be1),length(ga1));
cpsi2=morsea(ga,be).^2.*frac(1,ga.*2.^(2*be./ga)).*gamma(frac(2*be,ga));
cpsi2=reshape(cpsi2,length(be1),length(ga1));
[p,dt,dom]=morsebox(ga,be);

reporttest('MORSEWAVE admissibility matches analytic expression',aresame(cpsi1,cpsi2,1e-2));




function[]=morsewave_test_fft
fs=1./logspace(log10(5),log10(40))'; 
N=1023;
[psi,Psi]=morsewave(N,2,2,4,fs,'bandpass');
reporttest('MORSEWAVE Fourier transform, N odd, FS positive, K=2',aresame(abs(fft(psi,[],1)),abs(Psi),1e-8))
[psi,Psi]=morsewave(N+1,2,2,4,fs,'bandpass');
reporttest('MORSEWAVE Fourier transform, N even, FS positive, K=2',aresame(abs(fft(psi,[],1)),abs(Psi),1e-8))
[psi,Psi]=morsewave(N,2,1,4,-fs,'bandpass');
reporttest('MORSEWAVE Fourier transform, N odd, FS negative, K=2',aresame(abs(fft(psi,[],1)),abs(Psi),1e-8))
[psi,Psi]=morsewave(N+1,2,2,4,-fs,'bandpass');
reporttest('MORSEWAVE Fourier transform, N even, FS negative, K=2',aresame(abs(fft(psi,[],1)),abs(Psi),1e-8))

function[]=morsewave_test_notime
fs=1./logspace(log10(5),log10(40))'; 
N=1023;
tic
[psi1,Psi1]=morsewave(N,3,2,4,fs,'bandpass');
etime(1)=toc;
tic
Psi2=morsewave(N,3,2,4,fs,'bandpass','freq');
etime(2)=toc;

disp(['MORSEWAVE was ' num2str(etime(1)./etime(2)) ' times faster with frequency-only computation.'])
reporttest('MORSEWAVE frequency-only versus default',aresame(Psi1,Psi2,1e-8))

function[]=morsewave_fig
morsewave_figure1;
morsewave_figure3;


function[]=morsewave_figure1

N=256*4;

be=5;
ga=2;
K=3;
fs=2*pi/8/4;

[x,X]=morsewave(N,K,ga,be,fs,'energy');
%[fmin,fmax,fc,fw] = morsefreq(1/10000,ga,be);
%fc=fc./ao;

f=(0:1:N-1)'./N;

figure
t=(1:length(x))'-length(x)/2;
ax=[-60 60 -maxmax(abs(x))*1.05 maxmax(abs(x))*1.05];
subplot 321
  uvplot(t,x(:,1));axis(ax)
  title('Morse wavelets, time domain')
subplot 323
  uvplot(t,x(:,2));axis(ax)
subplot 325
  uvplot(t,x(:,3));axis(ax)

ax=[0 120./N -maxmax(abs(X))*1.05 maxmax(abs(X))*1.05];
subplot 322
  plot(f,abs(X(:,1))),axis(ax),vlines(fs);
title('Morse wavelets, frequency domain')
subplot 324
  plot(f,abs(X(:,2))),axis(ax),vlines(fs);
subplot 326
  plot(f,abs(X(:,3))),axis(ax),vlines(fs);
  
  

function[]=morsewave_figure2
fs=2*pi*.05;N=1000;


[x2,X2]=morsewave(N,3,2,4,fs,'energy','first');
[y2,Y2]=morsewave(N,3,2,4,fs,'energy','second');

[x3,X3]=morsewave(N,3,3,4,fs,'energy','first');
[y3,Y3]=morsewave(N,3,3,4,fs,'energy','second');

[x4,X4]=morsewave(N,3,4,4,fs,'energy','first');
[y4,Y4]=morsewave(N,3,4,4,fs,'energy','second');

figure
subplot(231),
plot((1:length(X2))/N,X2./maxmax(X2(:,1))), xlim([0 .14]),vlines(fs,'k:'),hlines(0,'k:')
title('First Family | \gamma=2, \beta=4'),ylim([-1.1 1.1])
subplot(232),
plot((1:length(X2))/N,X3./maxmax(X3(:,1))), xlim([0 .14]),vlines(fs,'k:'),hlines(0,'k:')
title('First Family | \gamma=3, \beta=4'),ylim([-1.1 1.1])
subplot(233),
plot((1:length(X2))/N,X4./maxmax(X4(:,1))), xlim([0 .14]),vlines(fs,'k:'),hlines(0,'k:')
title('First Family | \gamma=4, \beta=4'),ylim([-1.1 1.1])
subplot(234),
plot((1:length(X2))/N,Y2./maxmax(Y2(:,1))), xlim([0 .14]),vlines(fs,'k:'),hlines(0,'k:')
title('Second Family | \gamma=2, \beta=4'),ylim([-1.1 1.1])
subplot(235),
plot((1:length(X2))/N,Y3./maxmax(Y3(:,1))), xlim([0 .14]),vlines(fs,'k:'),hlines(0,'k:')
title('Second Family | \gamma=3, \beta=4'),ylim([-1.1 1.1])
subplot(236),
plot((1:length(X2))/N,Y4./maxmax(Y4(:,1))), xlim([0 .14]),vlines(fs,'k:'),hlines(0,'k:')
title('Second Family | \gamma=4, \beta=4'),ylim([-1.1 1.1])

packcols(2,3);

if 0
    cd_figures
    orient landscape
    fontsize 14 14 14 14
    print -depsc morsefamilies.eps
end


function[]=morsewave_figure3

clear psi psif
for i=1:25
    [psi(:,i),psif(:,i)]=morsewave(1000,1,i,0,1/10,'bandpass');
end


