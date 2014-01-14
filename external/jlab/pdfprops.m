function[mu,sigma,skew,kurt]=pdfprops(x,fx)
%PDFPROPS  Mean and variance associated with a probability distribution.
%
%   [MU,SIGMA]=PDFPROPS(X,FX), given a probability distribution
%   function FX over values X, returns the mean MU and the standard
%   deviation SIGMA. Each column of X must have uniform spacing.
%
%   The statistics are computed using a trapezoidal integration.
%   FX is multiplied by a constant so that it integrates to one.
%
%   [MU,SIGMA,SKEW,KURT]=PDFPROPS(X,FX) also retuns the skewness and 
%   the kurtosis, which are the third and fourth central moments, 
%   respectively normalized by the third and fourth powers of the 
%   standard deviation.  
%
%   'pdfprops --t' runs a test.
%
%   Usage:  [mu,sigma]=pdfprops(x,fx); 
%           [mu,sigma,skew,kurt]=pdfprops(x,fx);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001--2009 J.M. Lilly --- type 'help jlab_license' for details    
  
if strcmp(x,'--t')
    pdfprops_test;
    return
end

if jisrow(x)
  x=x(:);
end
if jisrow(fx)
  fx=fx(:);
end

if size(x,2)==1
   x=x*ones(size(fx(1,:)));
end


%dx=vrep(x(2,:)-x(1,:),size(x,1),1);
dx=x(2,:)-x(1,:);
fx=fx./vrep(trapint(fx,dx),size(x,1),1);
mu=real(trapint(fx.*x,dx));
%mu=trapint(fx.*x,dx);
murep=vrep(mu,size(x,1),1);
sigma=sqrt(trapint((x-murep).^2.*fx,dx));
if nargout>=3
   skew=trapint((x-murep).^3.*fx,dx);
   skew=skew./sigma.^3;
end
if nargout==4
   kurt=trapint((x-murep).^4.*fx,dx);
   kurt=kurt./sigma.^4;
end


% for i=1:size(fx,2)
%     if trapint(fx(:,i),dx(i))~=1
%         %disp('Normalizing FX to unit area.')
%         fx(:,i)=fx(:,i)./trapint(fx(:,i),dx(i));
%     end
% end
% 
% for i=1:size(fx,2)
%   mu(i,1)=real(trapint(fx(:,i).*x(:,i),dx(i)));
%   sigma(i,1)=sqrt(trapint((x(:,i)-mu(i,1)).^2.*fx(:,i),dx(i)));
% end
% 


function[y]=trapint(f,dx)
%Trapezoidal integration

fa=f;
fb=vshift(fa,1,1);
fa(1,:)=0;
fb(1,:)=0;
fa(end,:)=0;
fb(end,:)=0;
y=vsum(frac(fa+fb,2),1).*dx;
vswap(y,0,1e-10);

function[]=pdfprops_test
x=(-30:.001:30)';

mu0=2;
sigma0=5;

f=simplepdf(x,mu0,sigma0,'gaussian');  %f=simplepdf(x,mu,sig,flag)  
[mug,sigmag,skewg,kurtg]=pdfprops(x,f);

f=simplepdf(x,mu0,sigma0,'boxcar');  %f=simplepdf(x,mu,sig,flag)  
[mu,sigma]=pdfprops(x,f);
tol=1e-3;

bool(1)=aresame(mu,mu0,tol).*aresame(sigma,sigma0,tol);
bool(2)=aresame(mug,mu0,tol).*aresame(sigmag,sigma0,tol);
bool(3)=aresame(skewg,0,tol).*aresame(kurtg,3,tol);

reporttest('PDFPROPS with uniform pdf', bool(1));
reporttest('PDFPROPS with Gaussian pdf', bool(2));
reporttest('PDFPROPS Gaussian skewness=0, kurtosis=3', bool(3));

% %/********************************************************
% x=(-10:.001:10)';
% f=simplepdf(x,0,2,'gaussian');
% f(end/2:end)=2*f(end/2:end);
% f(1:end/2)=0;
% f=f./sum(f)./0.001;
% plot(x,cumsum(f*.001))
% %********************************************************
