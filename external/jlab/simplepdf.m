function[f]=simplepdf(x,mu,sig,flag)
%SIMPLEPDF  Gaussian, uniform, Cauchy, and exponential pdfs.
%
%   F=SIMPLEPDF(X,MU,SIG,'gaussian') computes a Gaussian pdf with mean
%   MU and standard deviation SIG.
%  
%   F=SIMPLEPDF(X,MU,SIG,'boxcar') computes a uniform pdf with mean MU
%   and standard deviation SIG.
%  
%   F=SIMPLEPDF(X,XO,ALPHA,'cauchy') computes a Cauchy pdf with location
%   parameter XO and scale parameter ALPHA.
%
%   F=SIMPLEPDF(X,BETA,'exponential') computes an exponential pdf with
%   scale parameter, hence mean and standard deviation, equal to BETA.
%
%   'simplepdf --f' generates a sample figure
%
%   Usage: f=simplepdf(x,mu,sig,'gaussian');
%          f=simplepdf(x,mu,sig,'boxcar');
%          f=simplepdf(x,xo,alpha,'cauchy');
%          f=simplepdf(x,beta,'exponential');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001--2008 J.M. Lilly --- type 'help jlab_license' for details    
  
warning('off','MATLAB:divideByZero')
  
if strcmp(x,'--f')
  simplepdf_fig;
  return
end
dx=x(2)-x(1);

if nargin==3
    flag=sig;
end

if nargin<3&&strcmp(flag,'exponential')||nargin<4&&~strcmp(flag,'exponential')
    error('Not enough input arguments.')
end

if strcmp(flag,'gaussian')
  f=exp(-(x-mu).^2./2./sig.^2)./sig./sqrt(2*pi);
elseif strcmp(flag,'boxcar')
  f=0*x;
  ia=min(find(x-mu>-3.4641*sig/2))-1;
  ib=min(find(x-mu>3.4641*sig/2));
  f(ia:ib)=1;
  f=f./vsum(f*dx,1);
elseif strcmp(flag,'cauchy')
  alpha=sig;
  f=frac(alpha./pi,(x-mu).^2 + alpha.^2);
elseif strcmp(flag,'exponential')
  f=frac(1,mu).*exp(-abs(x)./mu);
end


warning('on','MATLAB:divideByZero')

function[]=simplepdf_fig

x=(-100:.1:100)';
mu=25;
sig=10;
f=simplepdf(x,mu,sig,'gaussian');
%[mu2,sig2]=pdfprops(x,f);
figure,plot(x,f),vlines(mu,'r')
%a=conflimit(x,f,95);
%vlines(mu+a,'g'),vlines(mu-a,'g')
title('Gaussian with mean 25 and standard deviation 10')
