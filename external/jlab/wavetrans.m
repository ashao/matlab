function[varargout]=wavetrans(varargin)
%WAVETRANS  Wavelet transform.
%
%   W=WAVETRANS(X,PSI) computes the wavelet transform W of a dataset X
%   using wavelet maxtrix PSI. X is a column vector, and the columns PSI 
%   are time-domain wavelets at different scales.
%
%   X and PSI may both contain multiple components, in which case
%   PSIF(:,:,k) specifies the kth wavelet and X(:,n) specifies the nth data
%   component.  If there are K wavelets at J frequencies and N data points 
%   in M data vectors, W is of size N x J x M x K.  Note that W is always
%   squeezed to remove singleton dimensions.
%
%   If the wavelet is too long for the signal, complex INFs are returned.
%   ___________________________________________________________________
%
%   Boundary conditions
%
%   W=WAVETRANS(..., STR), where STR is a string, optionally specifies the
%   boundary condition to be imposed at the edges of the time series.  
%   Valid options for STR are 
%
%         STR = 'periodic' for periodic boundary conditions 
%         STR = 'zeros' for zero-padding beyond the endpoints 
%         STR = 'mirror' for reflecting the time series at both ends
%
%   The default value of STR is 'periodic', which means endpoints of the 
%   time series are implicitly joined to make a periodic signal. All 
%   boundary conditions take into account potential blocks of missing data,
%   marked by NaNs, at beginning and end of each column.  
%   ___________________________________________________________________
%
%   Generalized Morse wavelets
%
%   WAVETRANS(X,{K,GAMMA,BETA,F,STR},...), with a cell array input, uses 
%   the first K Generalized Morse Wavelets specified by the parameters 
%   GAMMA and BETA.  F is the frequency array and STR is an optional 
%   argument determine the normalization. See MORSEWAVE for details.  
%
%   Calling WAVETRANS in this form is numerically more efficient than
%   constructing the wavelets first and then calling WAVETRANS. 
%
%   If F is an array it should be sorted in decreasing order, from highest
%   frequency to lowest frequency.
%
%   For details on the generalized Morse wavelets, see MORSEWAVE and
%
%       Lilly and Olhede (2009).  Higher-order properties of analytic 
%           wavelets.  IEEE Trans. Sig. Proc., 57 (1), 146--160.
%
%   For general purpose use, set GAMMA=3 and choose BETA to be no smaller
%   than one.  Increase BETA to make your wavelet have more 'wiggles'.
%
%   With the generalized Morse wavelets, use MORSESPACE to easily choose
%   the frequency bins.
%   ___________________________________________________________________
%
%   Missing data
%
%   The data X may contain blocks of NANs at the beginning and/or end of 
%   each column, marking the absence of data.  In this case only the 
%   data series is taken to correspond to the block of finite data values,
%   and the boundary conditions are applied accordingly. The corresponding
%   portions of the transform matrix W are then also set to NANs. No NANs
%   may occur in the interior of the data series.
%   ___________________________________________________________________
%
%   Multiple input series
%
%   [W1,W2,...,WN]=WAVETRANS(X1,X2,...,XN,...) also works, where the XN are
%   all datasets of the same size.
%   ___________________________________________________________________
%
%   Detrending
%
%   Note that the data X is detrended before transforming.  This feature
%   is suppressed by WAVETRANS(...,[STR], 'nodetrend').
%   ___________________________________________________________________
%
%   Complex-valued data
%
%   The wavelet transform is normalized differently for complex-valued data
%   than for real-valued data.
%
%   If WX and WY are the wavelet transforms of two real-valued signals, X 
%   and Y, then 
% 
%        WP=WAVETRANS(X+iY,PSI)   = (1/SQRT(2))*(WX + i WY)
%        WN=WAVETRANS(X-iY,PSI)   = (1/SQRT(2))*(WX - i WY)
%
%   defines the positive and negative rotary transforms WP and WN.  
%
%   The factors of SQRT(2) are included such that the transform power is
%   unchanged, that is, ABS(WX).^2+ABS(WY).^2 = ABS(WP).^2+ABS(WN).^2.
%
%   [WP,WN]=VECTMULT(TMAT,WX,WY) alternatively gives the rotary tranforms 
%   WP and WN from a unitary transformation of WX and WY with matrix TMAT.
%   ___________________________________________________________________
%
%   See also MORSESPACE, RIDGEWALK, WAVESPECPLOT.
%
%   Usage:  w=wavetrans(x,psi);
%           w=wavetrans(x,{k,gamma,beta,f,str});
%           w=wavetrans(x,{k,gamma,beta,f,str},str);
%           [wx,wy]=wavetrans(x,y,{k,gamma,beta,f,str},str);
%
%   'wavetrans --t' runs some tests.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2011 J.M. Lilly --- type 'help jlab_license' for details

if strcmp(varargin{1},'--t')
   wavetrans_test;return
end


bool=zeros(nargin,1);
for i=1:nargin
    bool(i)=isstr(varargin{i});
end
if ~allall(bool==0);
    N=find(bool,1,'first')-2;
else
    N=nargin-1;
end
argcell=varargin(N+1:end);
for j=1:N
    x=varargin{j}; 
    if iscell(x)
        for i=1:length(x)
            disp(['WAVETRANS transforming time series ' int2str(i) ' of ' int2str(length(x)) '.'])
            T{i}=wavetrans_one(x{i},argcell);
        end
    else
        T=wavetrans_one(x,argcell);
    end
    varargout{j}=T;
end
    
function[T]=wavetrans_one(x,argcell)

%Unitary transform normalization
if ~isreal(x)
    x=x./sqrt(2);
end

w=argcell{1};
bdetrend=1;
str='periodic';

if ischar(argcell{end})
  if strcmp(argcell{end}(1:3),'nod')
     bdetrend=0;
     if ischar(argcell{end-1})
        str=argcell{end-1};
     end  
  else
     str=argcell{end};
  end
end


x0=x;
M0=size(x0,1);


if isreal(x0)
    x=timeseries_boundary(x0,str,bdetrend);
else
    x=           timeseries_boundary(real(x0),str,bdetrend)+...
        sqrt(-1)*timeseries_boundary(imag(x0),str,bdetrend);
end
M=size(x,1);
N=size(x,2);
W=[];



if ~iscell(w)
    if size(w,1)>M0 && strcmp(str,'periodic')  
        disp('Data length must exceed the filter length---returning INFs.')
    elseif size(w,1)>3*M0
        disp('Data length must exceed one-third the filter length---returning INFs.')
    end
end


%/********************************************************
if iscell(w)
    if length(w)==4
        w{5}=[];
    end
    wp=w;
    %W=morsewave(M,wp{1},wp{2},wp{3},wp{4},wp{5},'frequency','second');
    W=morsewave(M,wp{1},wp{2},wp{3},wp{4},wp{5},'frequency','first');
    K=size(W,3);
    L=size(W,2);
else
    K=size(w,3);
    L=size(w,2);

    %Generate a frequency-domain wavelet matrix of same size as data
    if size(w,1)<M   %Some subtlety here for even/odd or odd/even
       wnew=zeros(M,L,K);
       index=(1:size(w,1))'+floor((M-size(w,1))./2);
       wnew(index,:,:)=w;
       w=wnew;
    elseif size(w,1)>M
       w=nan*w(1:M,:,:,:);  %If wavelet is too long, just truncate
    end
    
    W=fft(w);
    f=linspace(0,1-1./M,M)';
    f=vrep(f,L,2);
    f=vrep(f,K,3);
    W=W.*rot(-2*pi*f*(M+1)/2).*sign(1/2-f); %ensures wavelets are centered
    %Note, the sign function you need when the wavelets are real-valued
    W=reshape(W,M,L,K);
end
%\********************************************************

W=conj(W);
X=fft(x);
T=nan*ones(M0,L,N,K);


if M0~=M
    index=M0+1:M0*2;
else 
    index=(1:M0);
end

if 0
%Cleverer but not faster
Xmat=vrep(vrep(permute(X,[1 3 2]),L,2),K,4);
Wmat=vrep(permute(W,[1 2 4 3]),N,3);
T=ifft(Xmat.*Wmat);
T=T(index,:,:,:,:);
end

if 1
for k=1:K
  for n=1:N;
     %Xmat=vrep(X(:,n),L,2);
     Xmat=osum(X(:,n),zeros(L,1));
     Ttemp=ifft(Xmat.*W(:,:,k));
     T(:,:,n,k)=Ttemp(index,:);
  end
end
end

T=squeeze(T);

if isreal(x0)&&~iscell(w)
    if isreal(w)
        T=real(T);  %Strip imaginary part if real wavelet and real signal 
    end
end

if ~anyany(isfinite(T))
    T=inf*(1+sqrt(-1))*ones(size(T));
end


%/********************************************************
%Set missing data back to NANs
for i=1:size(x,2)
  index=find(isnan(x0(:,i)));
  if~isempty(index)
       Ttemp=T(:,:,i);
       Ttemp(index,:)=nan*(1+sqrt(-1));
       T(:,:,i)=Ttemp;
    end
end
%\********************************************************

function[]=wavetrans_test
wavetrans_test_centered;
wavetrans_test_sizes;
wavetrans_test_complex;
wavetrans_test_boundary;
wavetrans_test_tooshort;

function[]=wavetrans_test_tooshort
load npg2006
use npg2006

%Decide on frequencies
fs=2*pi./(logspace(log10(10),log10(100),50)');

%Compute wavelet transforms using generalized Morse wavelets
psi=morsewave(length(cx)+10,1,2,4,fs,'bandpass');
wx=wavetrans(real(cx),psi,'periodic');
reporttest('WAVETRANS returns INFs when signal is too short, periodic',allall(~isfinite(wx)))
psi=morsewave(3*length(cx)+10,1,2,4,fs,'bandpass');
wx=wavetrans(real(cx),psi,'mirror');
reporttest('WAVETRANS returns INFs when signal is too short, mirror',allall(~isfinite(wx)))


function[]=wavetrans_test_boundary
load npg2006
use npg2006

%Decide on frequencies
fs=2*pi./(logspace(log10(10),log10(100),50)');

%Compute wavelet transforms using generalized Morse wavelets
wx=wavetrans(real(cx),{1,2,4,fs,'bandpass'},'mirror');
wx2=wavetrans(real(cx),{1,2,4,fs,'bandpass'},'periodic');
wx3=wavetrans(real(cx),{1,2,4,fs,'bandpass'},'zeros');

res1=abs(wx-wx2);
res2=abs(wx-wx3);

reporttest('WAVETRANS mirror and periodic boundary conditions match in interior',allall(res1(200:920)<1e-2))
reporttest('WAVETRANS mirror and zero boundary conditions match in interior',allall(res2(200:920)<1e-2))


function[]=wavetrans_test_complex
load npg2006
use npg2006

%Decide on frequencies
fs=2*pi./(logspace(log10(10),log10(100),50)');

%Compute wavelet transforms using generalized Morse wavelets
[wx,wy]=wavetrans(real(cx),imag(cx),{1,2,4,fs,'bandpass'},'mirror');
[wp,wn]=wavetrans(cx,conj(cx),{1,2,4,fs,'bandpass'},'mirror');
[wp2,wn2]=vectmult(tmat,wx,wy);


reporttest('WAVETRANS complex-valued input',aresame(wp,wp2,1e-6)&&aresame(wn,wn2,1e-6))
function[]=wavetrans_test_sizes

x=testseries_lillypark_array;
N=size(x,1);
M=size(x,2);
%Calculate wavelet matrix
J=5;
K=3;
fs=1./(logspace(log10(20),log10(600),J)');
psi=morsewave(N,K,2,4,fs,'bandpass');
%Compute wavelet transforms
wx=wavetrans(x,psi,'mirror'); 
[N2,J2,M2,K2]=size(wx);

bool=aresame([N,J,K,M],[N2,J2,K2,M2]);

reporttest('WAVETRANS output matrix has size N x J x M x K',bool)

function[]=wavetrans_test_centered

J=4;
ao=logspace(log10(5),log10(40),J);
x=zeros(2^10-1,1);t=(1:length(x))';
[w,f]=morsewave(length(x),1,2,4,ao);
x(2^9,1)=1;
y=wavetrans(x,w);
clear maxi
for i=1:size(y,2);
  maxi(i)=find(abs(y(:,i))==max(abs(y(:,i))));
end
b(1)=max(abs(maxi-2^9)<=1);
reporttest('WAVETRANS Morsewave transform has peak at delta-function',b(1))

function[x,t]=testseries_lillypark_array
[x1,t]=testseries_lillypark;
x1=anatrans(x1);
x1=x1(1:10:end);
t=t(1:10:end);
dom=pi/6;
p1=[1 2 3 3 2.5 1.5]'.*rot(dom*(0:5)');
p1=p1./sqrt(p1'*p1);
x=real(oprod(x1,p1));
x=10*x./maxmax(x);

function[x,t,xo]=testseries_lillypark
t=(1:4096)';
om=zeros(size(t));
x2=0*om;
index=1000:3000;
om(index)=(t(index))*2/100000 - .019999;
x=sin(om.*t);
x2(index)=-sin((t(index)-1600)*2*pi*1/2800);
x=x.*x2;
x(3500:3515)=-1;
x(3516:3530)=1;
t=linspace(-100,100,length(x))';
xo=x2;





