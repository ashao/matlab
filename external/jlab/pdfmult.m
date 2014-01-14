function[fz]=pdfmult(xi,yi,fx,fy,zi)
%PDFMULT  Probability distribution from multiplying two random variables.
%
%   FZ=PDFMULT(XI,YI,FX,FY,ZI), given two probability distribution
%   functions FX and FY defined over XI and YI, returns the pdf FZ
%   corresponding to Z=X*Y over values ZI.
%
%   PDFMULT uses PDFDIVIDE.
%
%   Usage: yn=pdfmult(xi,yi,fx,fy,zi);
%
%   'pdfmult --f' generates a sample figure.  
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001--2009 J.M. Lilly --- type 'help jlab_license' for details    
  
if strcmp(xi,'--f')
    pdfmult_fig;
    return
end

vcolon(xi,yi,fx,fy,zi);
fyi=jinterp(yi,fy,zi);
fxi=jinterp(xi,fx,zi);
vswap(fyi,nan,0);
vswap(fxi,nan,0);
fyinv=pdfinv(zi,fyi);
fz=pdfdivide(zi,zi,fxi,fyinv,zi);
[mu,sigma]=pdfprops(zi,fz);

function[]=pdfmult_fig

dx=0.1;
dy=0.05;
dz=0.025;
s1=1;
s2=2;
xi=(-10:dx:10)';
yi=(-10:dy:10)';
zi=(-10:dz:10)';

fx=simplepdf(xi,0,s1,'gaussian');
fy=simplepdf(yi,0,s2,'gaussian');

fz0=s1.*s2./pi./(s2.^2.*xi.^2+s1.^2);
fz0=fz0./vsum(fz0*dx,1);
fz=pdfmult(xi,yi,fx,fy,zi);

figure,plot(zi,fz),hold on,plot(xi,fx)
plot(yi,fy)
linestyle default
title('RV with green pdf multiplied by RV with red pdf equals RV with blue pdf')
x1=randn(100000,1)*s1;
y1=randn(100000,1)*s2;
[fz1,n]=hist(x1.*y1,(-10:.1:10));
plot(n,fz1/10000,'.')

text(4,0.4,'Green and red are Gaussian')
text(4,0.35,'Blue is disribution of product')
text(4,0.30,'Dots are from a random trial')

axis([-10 10 0 0.45])
