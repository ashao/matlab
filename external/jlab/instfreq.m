function[varargout]=instfreq(varargin)
%INSTFREQ  Instantaneous frequency, bandwidth, and their generalizations.
%
%   OMEGA=INSTFREQ(X), where X is an analytic signal, computes the
%   instantaneous *radian* frequency OMEGA using the first central 
%   difference and assuming a unit sample rate. 
%
%   [OMEGA,UPSILON,XI]=INSTFREQ(X) also returns the instantaneous bandwidth
%   UPSILON and the instantaneous curvature XI.  
%
%   X is an array with the first dimension being "time".  Thus, X can be a 
%   matrix of analytic signals oriented as column vectors, or a 2- or 3-D 
%   wavelet transform such as output by WAVETRANS.
%
%   The output arrays are the same size as X. 
%
%   The instantaneous frequency, bandwidth, and curvature are defined as
%
%          OMEGA  = d/dt Im ln X = d/dt arg X
%         UPSILON = d/dt Re ln X = d/dt ln abs X 
%            XI   = d^2/dt^2 abs X / abs X + i d^2/dt^2 arg X
%                 = UPSILON^2 + d/dt UPSILON + i d/dt OMEGA      
%  
%   where i=SQRT(-1) as usual.
%   _________________________________________________________________
%
%   Multivariate signals
% 
%   INSTFREQ can also handle multivariate signals.
%
%   OMEGA=INSTFREQ(X,DIM), where X has more than one entry along dimension
%   DIM, will compute the joint instantaneous frequency OMEGA.  
%
%   The joint instantaneous frequency is computed by a power-weighted 
%   average of the component instantaneous frequencies along dimension DIM.
%   OMEGA therefore has the same size as X, except along dimension DIM 
%   where OMEGA has only one entry.
%
%   [OMEGA,UPSILON,XI]=INSTFREQ(X,DIM) also works, where UPSILON and XI
%   are now the joint instantaneous bandwidth and curvature, respectively.   
%   _____________________________________________________________________
%   
%   Sample interval
%
%   INSTFREQ(DT,...) uses sample interval DT, where DT is a scalar, for 
%   computing time derivatives.  DT=1 is the default.
%   _____________________________________________________________________
%
%   First and last points
%
%   The first and last points must be treated differently, as the central 
%   difference is not defined there.  Three different methods can be used.
%
%   INSTFREQ(...,STR) specifies the method: STR= 'endpoint' (the default),
%   'periodic', or 'nans'.  See VDIFF for details.  
%   _____________________________________________________________________
%
%   'instfreq --f' generates some sample figures.
%
%   See also JOINTFREQ.
%
%   Usage: om=instfreq(x);
%          [om,up,c]=instfreq(dt,x);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2011 J.M. Lilly --- type 'help jlab_license' for details


%   _____________________________________________________________________
%
%   Higher-order modulation functions
%
%   [OMEGA,RHO1,...RHON]=INSTFREQ(X) also outputs the higher-order 
%   instantaneous modulation functions.  RHO1 is identical to the
%   bandwidth, and RHO2 is called the curvature.
%
%   For details see 
%
%      Lilly and Olhede (2010).  On the analytic wavelet transform.
%
%   Note that the modulation functions are defined in their non-normalized 
%   form, that is, not divided by powers of the instananeous frequency, 
%   unlike in Lilly and Olhede (2010).
%
%   _____________________________________________________________________


if strcmp(varargin{1}, '--f')
    instfreq_figure,return
elseif strcmp(varargin{1}, '--f1')
    instfreq_figure1,return
elseif strcmp(varargin{1}, '--f2')
    instfreq_figure2,return
end

if length(varargin{1})==1
    dt=varargin{1};
    x=varargin{2};
else
    dt=1;
    x=varargin{1};
end

str='endpoint';
dim=[];
for i=1:2
    if isstr(varargin{end})
        str=varargin{end};
        varargin=varargin(1:end-1);
    else
        if length(varargin{end})==1
            dim=varargin{end};
            varargin=varargin(1:end-1);
        end
    end
end    



if isreal(x)
    error('X should be an analytic signal, which implies it should be complex-valued.')
end


%Note: there are different ways to reasonably define the instantaneous
%frequency as a vdiff of the signal.  These differ greatly actually. 
%Diffing the unwrapped angle seems to give the best results.  If you take 
%the Im of the diffed log, you will break the wavelet phase algorithm!

%You definitely don't want to just do diff(x).  This differences the real
%and imaginary parts separately, which are themselves rapidly varying.

if strcmp(str(1:3),'per')
    om=vdiff(unwrap(angle([x(end,:);x;x(1,:)])),1,str)./dt;
    om=om(2:end-1,:);
else
    om=vdiff(unwrap(angle(x)),1,str)./dt;
end

up=vdiff(log(abs(x)),1,str)./dt; 
eta=om-sqrt(-1)*up;

varargout{1}=om;
nmax=nargout-1;

%I prefer this version, though difference is minor
%if nmax>=1
%   varargout(3)=frac(1,abs(x)).*vdiff(vdiff(abs(x),1,str),1,str)./dt./dt+sqrt(-1)*vdiff(om,1,str)./dt; 
%end

if nmax>=1
    etadiff=cell(nmax,1);   
    polyargs=cell(nmax,1);
    bcell=cell(nmax,1);

    etadiff{1}=vdiff(eta,1,str)./dt;
    polyargs{1}=up;

    for n=2:nmax
        etadiff{n}=vdiff(etadiff{n-1},1,str)./dt;
        polyargs{n}=sqrt(-1)*etadiff{n-1};
    end
    varargout(2:nmax+1)=bellpoly(polyargs);
end

%When DIM is input, this tells me to output joint quantities
if ~isempty(dim)
    if nargout==1
        varargout{1}=jointfreq(x,varargout{1},dim);
    elseif nargout==2
        [varargout{1},varargout{2}]=jointfreq(x,varargout{1},varargout{2},dim);
    elseif nargout==3
        [varargout{1},varargout{2},varargout{3}]=jointfreq(x,varargout{1},varargout{2},varargout{3},dim);
    end
end
    
    

function[]=instfreq_figure
instfreq_figure1
instfreq_figure2
 
function[]=instfreq_figure1
dt=0.5;
t=(-200:dt:200)';

x=morsewave(length(t),1,24,2,2*pi/20.*dt);
%x=10*morsexpand(10000,t,24,2,2*pi/20);
x=x./maxmax(abs(x));
%x=10*morsexpand(200,t,7,2,2*pi/20);
%x(abs(t)>100)=nan;
[om,rho1,rho2,rho3,rho4]=instfreq(dt,x);
eta=om-sqrt(-1)*rho1;


figure,
subplot(231),plot(t,abs(x)),hold on,uvplot(t,x); ylim([-1.1 1.1]),linestyle 2k k k-- k:
title('Analytic Signal \psi_{24,2}(t)'),ytick(-.8:.4:.8),fixlabels([0 -1]), xlim([-100 100])
subplot(232),uvplot(t,eta);ylim([-.45 .45]),linestyle k k-- k:
title('Complex Instantaneous Frequency \eta(t)'),xlim([-100 100])
subplot(233),plot(t,real(rho1)./om),ylim([-.25 .25]),linestyle k k:
title('Modulation Function #1 \rho_1(t)'),xlim([-100 100])
subplot(234),plot(t,abs(rho2)./om.^2),hold on,uvplot(t,rho2./om.^2),ylim([-.1 .1]),linestyle 2k k k-- k:
title('Modulation Function #2 \rho_2(t)'),ytick(-.15:.05:.15),xlim([-100 100])
%subplot(235),uvplot(t,rho1.*rho2./om.^3),ylim([-.1 .1]),linestyle  k k-- k:
%title('Cross-Term 3'),ytick(-.15:.05:.15),
subplot(235),plot(t,abs(rho3)./om.^3),hold on,uvplot(t,rho3./om.^3),ylim([-.1 .1]),linestyle 2k k k-- k:
title('Modulation Function #3 \rho_3(t)'),ytick(-.15:.05:.15),xlim([-100 100])

subplot(236),plot(t,abs(rho4)./om.^4),hold on,uvplot(t,rho4./om.^4),ylim([-.1 .1]),hlines(0),linestyle 2k k k-- k:
title('Modulation Function #4 \rho_4(t)'),ytick(-.15:.05:.15),xlim([-100 100])


for i=1:6
    subplot(2,3,i),vlines([-22 22],'k:')
    hlines(0,'k:'),
    xlim([-75 75]),xtick(-60:20:60),fixlabels([0 -2]), 
   % if i==3,hlines([-.2 .2],'k--'),end
end

letterlabels(4)

% omm=powermean(real(eta),x,1);
% mu3=[imag(rho2).*rho1 (real(eta)-omm).*(2*rho1.^2-real(rho2)) (real(eta)-omm).^3];
% mu3=[mu3 vsum(mu3,2)];
%  figure,plot(t,mu3.*[abs(x).^2 abs(x).^2 abs(x).^2 abs(x).^2])
% 
%  
% mu2=[rho1.^2 (real(eta)-omm).^2];
% mu2=[mu2 vsum(mu2,2)];
%  figure,plot(t,mu2.*[abs(x).^2 abs(x).^2 abs(x).^2 ])


function[]=instfreq_figure2

dt=0.5;
t=(-100:dt:100)';
x=10*morsexpand(2000,t,24,2,2*pi/20);
x=x./maxmax(abs(x));

[om,rhon{1},rhon{2},rhon{3},rhon{4},rhon{5},rhon{6}]=instfreq(dt,x);
eta=om-sqrt(-1)*rhon{1};


tmid=round(length(t)/2);

xo=x(tmid);
omo=om(tmid);


xhat=zeros(length(x),7);
xhat(:,1)=xo.*rot(t.*omo);
for n=1:6
    xhat(:,n+1)=xhat(:,1).*frac(1,factorial(n)).*(t.^n).*rhon{n}(tmid);
end

xhat=cumsum(xhat,2);

figure
subplot(121),plot(t,real(x)),hold on,plot(t,real(xhat(:,3:2:end))),xlim([-35 35]),hlines(0)
ylim([-1.1 1.1]),ytick(-.8:.4:.8),xtick(-60:20:60),fixlabels([0 -1]),linestyle 2k k k-- k-. k: 
title('Demodulate Expansion of \Re\{\psi_{24,2}(t)\} about t=0'),vlines([-22 22],'k:')
vlines(0,'k:')
%plot(t(tmid),real(x(tmid)),'*')

subplot(122),plot(t,imag(x)),hold on,plot(t,imag(xhat(:,3:2:end))),xlim([-35 35]),hlines(0)
ylim([-1.1 1.1]),ytick(-.8:.4:.8),xtick(-60:20:60),fixlabels([0 -1]),linestyle 2k k k-- k-. k:
title('Demodulate Expansion of \Im\{\psi_{24,2}(t)\} about t=0'),vlines([-22 22],'k:')
vlines(0,'k:')
%plot(t(tmid),imag(x(tmid)),'*')

letterlabels(4),packcols(1,2)



