function[zp,zn]=anatrans(x,str)
%ANATRANS  Analytic part of signal.
%
%   Z=ANATRANS(X) returns the analytic part of the real-valued signal X,
%   which is a column vector or a matrix with 'time' in columns. 
%
%   Z is defined for real X such that X = REAL(Z) = 1/2*(Z+CONJ(Z)).
%
%   [ZP,ZN]=ANATRANS(Z) returns the analytic part ZP and conjugate-analytic 
%   part ZN of the complex-valued signal Z.  
%
%   ZP and ZN are then defined such that Z=ZP+CONJ(ZN). Both ZP and ZN are
%   analytic signals.
% 
%   ANATRANS does not detrend, so if you time series has a mean value or 
%   a linear trend you should call DETREND first. 
%
%   The analytic part of a real-valued signal and that of a complex-valued 
%   signal are thus defined differently by a factor of two, following the 
%   convention of Lilly and Gascard (2006) and Lilly and Olhede (2010).
%   ___________________________________________________________________
%
%   Boundary conditions
%
%   ANATRANS(X STR), where STR is a string, optionally specifies the
%   boundary condition to be imposed at the edges of the time series.  
%   Valid options for STR are 
%
%         STR = 'periodic' for periodic boundary conditions 
%         STR = 'zeros' for zero-padding beyond the endpoints 
%         STR = 'mirror' for reflecting the time series at both ends
%
%   The default value of STR is 'mirror', as this removes `edge effects' 
%   which otherwise tend to occur near the ends of the time series.
%   ___________________________________________________________________
%
%   'anatrans --t' runs some tests.
%
%   Usage: z=anatrans(x);
%          z=anatrans(x,'mirror');
%          [zp,zn]=anatrans(z);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2011 J.M. Lilly --- type 'help jlab_license' for details        



if strcmp(x, '--t')
    anatrans_test,return
end

% antrans is basically the same as Matlab's Hilbert:
% if isreal(x)
%     zp=x+sqrt(-1)*hiltrans(x);
% else
%     zp=frac(1,2)*(x+sqrt(-1)*hiltrans(x));
%     zn=frac(1,2)*(x-sqrt(-1)*hiltrans(x));
% end

if nargin==1
    str='mirror';
end

M0=size(x,1);
x=timeseries_boundary(x,str,false);
M=size(x,1);

mx=vmean(x,1);
mx=vrep(mx,M,1);
x=x-mx;
Z=2*fft(x);
Z(round(end/2):end,:)=0;
zp=ifft(Z);
zp=zp+mx;

if M0~=M
    index=M0+1:M0*2;
    x=x(index,:,:);          
    zp=zp(index,:,:);          
end
    
if ~isreal(x)
    zp=frac(1,2)*zp;
    zn=conj(x-zp);
end


function[]=anatrans_test

load solomon 
use solomon

z=anatrans(x,'periodic');

res=frac(vsum(abs(real(z)-x).^2,1),vsum(abs(z).^2,1));
reporttest('ANATRANS departure of real part from real-valued original less than 1/1000, Solomon Islands',allall(res<1/1000))

[zp,zn]=anatrans(z,'periodic');
res=frac(vsum(abs(z-zp).^2,1),vsum(abs(z).^2,1));
reporttest('ANATRANS positive rotary part recovers input analytic signal, Solomon Islands',allall(res<1/1000))
res=vsum(abs(zn).^2,1);
reporttest('ANATRANS negative rotary part negligible for input analytic signal, Solomon Islands',allall(abs(zn)<1e-8))


