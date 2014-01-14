function[yn,mu,sigma]=pdfconv(x,y,n)
%PDFCONV  Convolution of a probability distribution with itself.
%
%   YN=PDFCONV(X,Y,N) given a probability distribution function Y over
%   values X, returns a matrix containing the pdfs resulting from
%   convolving Y with itself [0,1,2,...N-1] times.
%
%   [YN,MU,SIGMA]=PDFCONV(...) also returns an array of means MU and
%   standard deviations SIGMA associated with each column of YN.
%
%   'pdfconv --t' runs a test
%   'pdfconv --f' generates a sample figure
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001, 2004 J.M. Lilly --- type 'help jlab_license' for details    
  
if strcmp(x,'--t')
   pdfconv_test;
   return
end
if strcmp(x,'--f')
   pdfconv_fig;
   return
end

vcolon(x,y);

dx=x(2)-x(1);

y=y./vsum(y*dx,1);
yn=y;
yni=y;
for i=1:n-1
    yni=conv(yni,y);
    %xtemp=(x(1):dx/(i+1):x(end))';
    xtemp=linspace(x(1),x(end),length(yni))';
    %size(xtemp),size(yni),size(x)
    yn(:,i+1)=interp1(xtemp,yni,x);
    yn(:,i+1)=yn(:,i+1)./vsum(yn(:,i+1)*dx,1);
end

if 0
    for i=1:n-1
        yni=conv(yni,y);
        xtemp=(x(1):dx/(i+1):x(end))';
        yn(:,i+1)=interp1(xtemp,yni,x);
        yn(:,i+1)=yn(:,i+1)./sumnan(yn(:,i+1)*dx);
    end
end

[mu,sigma]=pdfprops(x*ones(size(yn(1,:))),yn);

function[x1,fn]=pdfconv_test
x1=(-15:.01:10)';
x2=(10:.01:40)';

f1=simplepdf(x1,-2,3,'gaussian');  %f=simplepdf(x,mu,sig,flag)  
%f2=simplepdf(x2,25,3,'gaussian');  %f=simplepdf(x,mu,sig,flag)  

dx=x1(2)-x1(1);
%fz=conv(f1,f2)';
%z=(1:length(f)]'*dx;
%z=z-z(1)+x1(1)+x2(1);
[fn,mu,sigma]=pdfconv(x1,f1,10);

sigma0=3./sqrt((1:length(sigma)));
mu0=-2+0*sigma0;
tol=0.005;

bool(1)=aresame(mu,mu0,tol);
bool(2)=aresame(sigma,sigma0,tol);

reporttest('PDFCONV has correct mean', bool(1));
reporttest('PDFCONV has correct standard deviation',bool(2));

function[]=pdfconv_fig
[x1,fn]=pdfconv_test;
figure,plot(x1,fn)
title('Convolution of a Gaussian with itself; mean=-2, \sigma= 3')


