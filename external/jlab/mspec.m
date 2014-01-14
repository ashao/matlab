function[varargout]=mspec(varargin)
% MSPEC  Multitaper power and cross spectra.
% 
%   MSPEC implements spectral and cross-spectral analysis using the 
%   multitaper method for real or complex-valued data. 
%
%   MSPEC is to be run after calling SLEPTAP to compute the multitapers.
%   ______________________________________________________________________
%
%   One-sided power spectrum
%
%   [F,S]=MSPEC(X,PSI) returns the one-sided power spectrum of X at 
%   positive frequencies using data tapers PSI. 
%  
%   Input:   
%       X  --  M x N matrix containing N length M time series
%     PSI  --  M x K matrix of K data tapers
% 
%   Output:
%       F  --  [M/2] nonnegative *cyclic* frequencies
%       S  --  [M/2] x N one-sided power spectrum matrix
%                
%   In the above, [M/2] means: M/2 if M is even, and (M+1)/2 if M is odd.  
%
%   The spectra matrices are averages over the K "eigenspectra" computed 
%   with each of the K tapers, as discussed in Park et al. JGR 1987.
%
%   The one-sided spectrum S is normalized such that its integral over all 
%   frequencies F, SUM(S,1)*DF where DF is the frequency increment,
%   approximates the signal variance. 
%   ______________________________________________________________________
%  
%   Cross-spectra of real-valued data
%   
%   MSPEC can be used to compute the cross-spectrum of two real-valued
%   time series or sets of time series.
%
%   [F,SXX,SYY,SXY]=MSPEC(X,Y,PSI);   --- For cross-spectra
%        
%   Input:
%       X  --  M x N matrix containing N length M time series
%       Y  --  M x N matrix containing N length M time series
%     PSI  --  M x K matrix of K data tapers
% 
%   Output:
%       F  --  M/2 nonnegative frequencies
%     SXX  --  [M/2] x N one-sided spectra of X
%     SYY  --  [M/2] x N one-sided spectra of Y
%     SXY  --  [M/2] x N one-sided cross spectra of X and Y 
%   ______________________________________________________________________
%  
%   Rotary spectra of complex-valued data
%
%   MSPEC can also the so called "rotary spectra" of complex-valued 
%   time series or sets of time series.
%
%   [F,SPP,SNN,SPN]=MSPEC(Z,PSI);   --- For rotary spectra of Z=X+iY
%
%   Input:   
%       Z  --  M x N matrix containing N length M time series
%     PSI  --  M x K matrix of K data tapers
%   
%   Output:   
%        F  --  M/2 nonnegative frequencies
%      SPP  --  [M/2] x N positively rotating power spectrum matrix
%      SNN  --  [M/2] x N negatively rotating power spectrum matrix  
%      SPN  --  [M/2] x N rotary cross spectral matrix      
%
%   Note that the rotary spectra are defined such that SXX+SYY=SPP+SNN.
%  
%   The rotary spectra SPP and SNN are normalized such that the sum of 
%   their integrals over all frequencies F approximates the variance of Z. 
%   ______________________________________________________________________
%  
%   Sample rate
%
%   [F,S]=MSPEC(DT,...) specifies the sample interval to be used in the
%   calculation of the frequency array F. DT defaults to unity.
%
%   Spectral values depend linearly upon the same rate in order that the 
%   integral of the spectra over frequency approximate the variance.
%   ______________________________________________________________________
%   
%   Cross-spectra of complex-valued data
%
%   To compute the cross-spectra of two complex-valued time series or sets 
%   of time series Z1 and Z2, run MSPEC repeatedly.
%
%   [F,SP1P1,SP2P2,SP1P2]=MSPEC(Z1,Z2,PSI);  
%   [F,SN1N1,SN2N2,SN1N2]=MSPEC(CONJ(Z1),CONJ(Z2),PSI);  
%
%   The first call returns the spectra and cross-spectra of Z1 and Z2 at
%   positive frequencies, while the second returns their spectra and cross-
%   spectra at negative frequencies.  Finally
%
%   [F,SP1P1,SN2N2,SP1N2]=MSPEC(Z1,CONJ(Z2),PSI);  
%
%   returns the cross-spectra between the positive rotary components of Z1 
%   and the negative components of Z2.
%   ______________________________________________________________________
%
%   Adaptive spectra
%
%   MSPEC(...,LAMBDA,'adaptive'), where LAMBDA contains the eigenvalues of
%   the tapers as computed by SLEPTAP, alternately uses the "adaptive"
%   multitaper method of Thomson (1982).
% 
%   This implementation follows that of Park et al. (1987a), JGR.
%
%   For cross-spectra or for rotary spectra, the weights appearing in the
%   adaptive spectra are derived for the total spectrum of each signal 
%   compoment, that is for SXX+SYY or SPP+SNN as appropriate.  Then the
%   separate spectra and co-spectra are computed using identical weights.
%   ______________________________________________________________________
%
%   'mspec --t' runs some tests.
%   'mspec --f' generates some sample figures from Bravo mooring data.
%
%   See also: SLEPTAP, HERMFUN, MTRANS, MSVD.
%
%   Usage   [f,s]=mspec(x,psi);    
%           [f,s]=mspec(dt,x,psi);     
%           [f,spp,snn,spn]=mspec(z,psi);     
%           [f,sxx,syy,sxy]=mspec(x,y,psi);
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2011 J.M. Lilly --- type 'help jlab_license' for details        
  


%   The Cartesian and rotary spectra are then related by a unitary
%   transformation.  For details see
%
%       Lilly and Olhede (2010).  Bivariate instantaneous frequency and
%           bandwidth.  IEEE Trans. Sig. Proc.

if strcmp(varargin{1},'--t')
    mspec_test;return
end
if strcmp(varargin{1},'--f')
    mspec_fig;return
end

if isscalar(varargin{1})
  deltat=varargin{1};
  varargin=varargin(2:end);
else
  deltat=1;
end
  

%Sort out input arguments
lambda=1;  %This means use the average multitaper spectrum
if isstr(varargin{end})
    if strfind(varargin{end},'ada')
        lambda=varargin{end-1};
        varargin=varargin(1:end-2);
    end
end

x=varargin{1};
psi=varargin{end};
na=length(varargin);

if na==2
    if isreal(x)
        y=[];
    else
        y=conj(x);
    end
elseif na==3
    y=varargin{2};
end


%   This behavior is removed for now 
%   ______________________________________________________________  
%  
%   Cell array input / output
%
%   MSPEC generates cell array output given cell array input.
%
%   For example, let's say one has P different time series arrays, X1, X2,
%   ..., XP, having lengths M1, M2, ..., MP respectively.  
%
%   Put these into a cell array X{1}=X1, X{2}=X2, ..., X{P}=XP and use
%   PSI=SLEPTAP([M1, M2,..., MP]) to make a cell array of tapers.
%
%   [F,S]=MSPEC(X,PSI) then returns cell arrays F and S corresponding 
%   to the Fourier frequencies and spectra of the P arrays.  
%
%   The other argument forms given above also work.
%
% if iscell(x)&&(iscell(y)||isempty(y))&&iscell(psi)
%     if (length(x)~=length(psi)&&isempty(y)) ||  (length(x)~=length(psi)&&length(y)~=length(psi)) 
%         error('With cell array input, all input arguments must be cell arrays of the same size.')
%     end
%     for i=1:length(x)
%         if ~isempty(y)
%             cellout=mspec_one(deltat,x{i},y{i},psi{i},lambda);
%         else
%             if isreal(x{1})
%                 cellout=mspec_one(deltat,x{i},[],psi{i},lambda);
%             else
%                 cellout=mspec_one(deltat,real(x{i}),imag(x{i}),psi{i},lambda);
%             end
%         end      
%         for j=1:length(cellout)
%             varargout{j}{i}=cellout{j};
%         end
%     end
%elseif ~iscell(x)&&~iscell(y)&&~iscell(psi)

varargout=mspec_one(deltat,x,y,psi,lambda);


    
function[cellout]=mspec_one(dt,x,y,psi,lambda)

%In real cases, multiply by sqrt(2) to get one-sided spectrum
if isempty(y)        %One real-valued
     [f,mmatx]=mtrans(sqrt(2)*x,psi);
else
     if ~isreal(x)   %Two complex-valued
        [f,mmatx,mmaty]=mtrans(x,y,psi);
     else            %Two real-valued
        [f,mmatx,mmaty]=mtrans(sqrt(2)*x,sqrt(2)*y,psi); 
     end
end


if isempty(y) %One time series
     if lambda==1
         cellout{2}=avgspec(mmatx,mmatx).*dt;
     else
         var=squared(vstd(x,1)); 
         cellout{2}=adaptspec(abs(mmatx).^2,lambda,var).*dt;
     end
else         %Two time series
     if lambda==1
        cellout{2}=avgspec(mmatx,mmatx).*dt;
        cellout{3}=avgspec(mmaty,mmaty).*dt;
        cellout{4}=avgspec(mmatx,mmaty).*dt;
     else
        %For two time series one should do the adaptive spectra on both
        %with the same coefficients
        var=squared(vstd(x,1))+squared(vstd(y,1)); 
        
        [s,dk]=adaptspec(abs(mmatx).^2+abs(mmaty).^2,lambda,var);
        cellout{2}=squeeze(frac(sum(dk.^2.*abs(mmatx).^2.*dt,2),sum(abs(dk).^2,2)));
        cellout{3}=squeeze(frac(sum(dk.^2.*abs(mmaty).^2.*dt,2),sum(abs(dk).^2,2)));
        cellout{4}=squeeze(frac(sum(dk.^2.*mmatx.*conj(mmaty).*dt,2),sum(abs(dk).^2,2))); 
     end
end

cellout{1}=f./dt;

function[S]=avgspec(mmat1,mmat2)
eigspec=mmat1.*conj(mmat2);
S=squeeze(mean(eigspec,2));



function[s,dk]=adaptspec(eigspec,lambda,var)

s=squeeze(zeros(size(eigspec(:,1,:))));
dk=zeros(size(eigspec));
for i=1:size(eigspec,3)
    [s(:,i),dk(:,:,i)]=adaptspec_one(eigspec(:,:,i),lambda,var(i));
end


function[s,dk]=adaptspec_one(eigspec,lambda,var)

s=frac(1,2)*squeeze(eigspec(:,1)+eigspec(:,2));
lambda=lambda(:)';
tol=1e-3;

var=var+0*s;
i=0;
sold=s*0+1000;
while any(abs(s-sold)./sold>tol)&&i<20
	i=i+1;
	sold=s;
    bk=var*(ones(size(lambda))-lambda);  %Outer product
	dk=(s*real(sqrt(lambda)))./(s*lambda+bk);  %Outer products
	s=frac(sum(dk.^2.*eigspec,2),sum(abs(dk).^2,2));
end
if i~=20
    disp(['Adaptive spectral estimate took ' int2str(i) ' iterations.'])
else
    disp(['Adaptive spectral loop terminated at ' int2str(i) ' iterations.'])
end

function[]=mspec_test

tol=1e-10;
%[x,t,xo]=testseries_bravoninetyfour;

load bravo94
use bravo.rcm

num=yf2num(bravo.rcm.yearf);

[psi,lambda]=sleptap(length(cv),8);

[f,sxx,syy,sxy]=mspec(real(cv),imag(cv),psi);
[f,spp,snn,spn]=mspec(cv,psi);

S(1,1,:,:)=sxx;
S(2,2,:,:)=syy;
S(1,2,:,:)=sxy;
S(2,1,:,:)=conj(sxy);

SZ(1,1,:,:)=spp;
SZ(2,2,:,:)=snn;
SZ(1,2,:,:)=spn;
SZ(2,1,:,:)=conj(spn);

T=vrep(vrep(tmat,size(S,3),3),size(S,4),4);

SZ2=matmult(matmult(T,S,1),conj(permute(T,[2 1 3 4])),1);

reporttest('MSPEC for (x,y) and (z,z^*) are unitary transforms with matrix T',aresame(SZ,SZ2,1e-10))


tol=1e-10;
[x,t,xo]=testseries_bravoninetyfour;
[psi,lambda]=sleptap(length(x),8);

p0=vsum(abs(real(x)-vmean(real(x),1)).^2,1)./length(x);
[f,sp]=mspec(real(x),psi);
p1=(vsum(sp,1)).*(f(2)-f(1));
reporttest('MSPEC satisfies Parseval''s theorem to within 4% for real Bravo, unit sample rate',abs(p1-p0)./p0<4/100);

[f,sp]=mspec(t(2)-t(1),real(x),psi);
p1=(vsum(sp,1)).*(f(2)-f(1));
reporttest('MSPEC satisfies Parseval''s theorem to within 4% for real Bravo, non-unit sample rate',abs(p1-p0)./p0<4/100);

p0=vsum(abs(x-vmean(x,1)).^2,1)./length(x);
[f,sp,sn,spn]=mspec(x,psi);
p1=(vsum(sp,1)+vsum(sn,1)).*(f(2)-f(1));
reporttest('MSPEC satisfies Parseval''s theorem to within 4% for complex Bravo, unit sample rate',abs(p1-p0)./p0<4/100);

[f,sp,sn,spn]=mspec(t(2)-t(1),x,psi);
p1=(vsum(sp,1)+vsum(sn,1)).*(f(2)-f(1));
reporttest('MSPEC satisfies Parseval''s theorem to within 4% for complex Bravo, non-unit sample rate',abs(p1-p0)./p0<4/100);

[f,sx,sy]=mspec(real(x),imag(x),psi);
p1=(vsum(sx+sy,1)).*(f(2)-f(1));
reporttest('MSPEC satisfies Parseval''s theorem to within 4% for bivariate Bravo, unit sample rate',abs(p1-p0)./p0<4/100);

[f,sx,sy]=mspec(t(2)-t(1),real(x),imag(x),psi);
p1=(vsum(sx+sy,1)).*(f(2)-f(1));
reporttest('MSPEC satisfies Parseval''s theorem to within 4% for bivariate Bravo, non-unit sample rate',abs(p1-p0)./p0<4/100);

reporttest('MSPEC SXX+SYY=SPP+SNN for Bravo, non-unit sample rate',aresame(sx+sy,sp+sn,tol));


mspec_test_frequency;

function[]=mspec_test_frequency


T=10000;
fo=100./T;
x=cos([1:T]'*2*pi*fo);
[psi,lambda]=sleptap(T,1,1);
[f,sp,sn]=mspec(x+sqrt(-1)./1e10,psi);

%figure,plot(f,sp,'r.'),hold on,plot(f,sn,'bo'), xlog, ylog, vlines(fo)
[mp,jp]=max(sp);
[mn,jn]=max(sn);
bool=aresame(f(jp),fo)&&aresame(f(jn),fo);
reporttest('MSPEC frequency matches expected exactly, even number of points',bool);
%That's 5% of the Rayleigh frequency
%length(f)

T=10000-1;
fo=100/T;
x=cos([1:T]'*2*pi*fo);
[psi,lambda]=sleptap(T,1,1);
[f,sp,sn]=mspec(x+sqrt(-1)./1e10,psi);
%length(f)

%figure,plot(f,sp,'r.'),hold on,plot(f,sn,'bo'), xlog, ylog, vlines(fo)
[mp,jp]=max(sp);
[mn,jn]=max(sn);
bool=aresame(f(jp),fo)&&aresame(f(jn),fo);
reporttest('MSPEC frequency matches expected exactly, odd number of points',bool);
%That's 5% of the Rayleigh frequency



function[]=mspec_test_former

tol=1e-10;
[x,t,xo]=testseries_bravoninetyfour;
L=[200 512 1024];
[psi,lambda]=sleptap(L,8);

for i=1:3
    xcell{i}=x(1:L(i),1);
    [fcell{i},spcell{i},sncell{i},spncell{i}]=mspec(xcell{i},psi{i});
end

[fcell2,spcell2,sncell2,spncell2]=mspec(xcell,psi);

for i=1:3
    bool(i,1)=aresame(fcell{i},fcell2{i},tol);
    bool(i,2)=aresame(spcell{i},spcell2{i},tol);
    bool(i,3)=aresame(sncell{i},sncell2{i},tol);
    bool(i,4)=aresame(spncell{i},spncell2{i},tol);
end

 
reporttest('MSPEC cell array input and output',allall(bool))  


function[]=mspec_fig
load bravo94
x=bravo.rcm.cv;
vswap(x,nan,0);
[psi,lambda]=sleptap(length(x),16);
[f,sp,sn,spn]=mspec(x,psi);
[f,su,sv,suv]=mspec(real(x),imag(x),psi);

figure,plot(f,[sp sn]),xlog,ylog
title('Counterclockwise (blue) and clockwise (red) spectra'),
linestyle b b b b b b r r r r r r  

load bravo94
x=bravo.rcm.cv;
vswap(x,nan,0);
[psi,lambda]=sleptap(length(x),8);

for i=1:size(x,2)
  [f,Suu,Suu3,Cuu]=mspec(real(x(:,i)),real(x(:,3)),psi);
  gammauu(:,i)=Cuu./sqrt(Suu.*Suu3);
end
figure,
plot(f,abs(gammauu)),xlog,yoffset 1
title('Coherence of u(t) at each depth vs. u(t) at #3')


function[x,t,xo]=testseries_bravoninetyfour
load bravo94
use bravo
x=bravo.rcm.cv(:,3);
xo=vfilt(x,24,'nonans');
num=yf2num(bravo.rcm.yearf);
t=num-yf2num(floor(bravo.rcm.yearf(1)));


% [f,Cuv]=mspec(real(cv),imag(x),psi);
% [f,Suu]=mspec(real(cv),psi);
% [f,Svv]=mspec(imag(cv),psi);
% 
% gammauv=Cuv;
% for i=1:size(Suu,2)
%   gammauv(:,i)=Cuv(:,i)./sqrt(Suu(:,i).*Svv(:,3));
% end
% figure,
% 
% 
% plot(f,abs(gammauv)),xlog,yoffset 1
% title('Cross-spectrum of u(t) at each depth vs. u(t) at #3')

