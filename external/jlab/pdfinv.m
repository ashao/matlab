function[fz]=pdfinv(yi,fy)
%PDFINV  Probability distribution of the inverse of a random variable.
%
%   YN=PDFMULT(YI,FY) given a probability distribution functions FY
%   defined over YI, returns the pdf of the inverse random variable 1/Y.
%
%   'pdfinv --t' runs a test
%   'pdfinv --f' generates a sample figure
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information (C)
%   2001, 2004 J.M. Lilly --- type 'help jlab_license' for details
  
if strcmp(yi,'--t')
  pdfinv_test;
  return
end

if strcmp(yi,'--f')
  pdfinv_fig;
  return
end

%tol=1e-10;
%index=find(yi==0);
%if ~isempty(index)
%  yi(index)=1e-10;
%end

index=1:round(length(yi)/2)-1;
N=round(length(yi)/2)-1;

%index=(N:-1:1  length(yi):-1:N+1);
%fz=pdfchain(yi(index),fy(index),1./yi(index),yi);
warning('off','MATLAB:divideByZero')
fz=pdfchain(yi,fy,1./yi,yi);
warning('on','MATLAB:divideByZero')

vswap(fz,nan,0);

function[]=pdfinv_test
  
alpha=2;

dy=0.01;
yi=(-40:dy:40)';
fx=simplepdf(yi,0,alpha,'cauchy');
fy=simplepdf(yi,0,1./alpha,'cauchy');
fy2=pdfinv(yi,fx);

tol=1e-3;
bool=vmean(abs(fy-fy2).^2,1)<tol;
reporttest('PDFINV for Cauchy, Papoulis special case p. 94',bool)

function[]=pdfinv_fig

s2=2;
dy=0.01;
yi=(-40:dy:40)';
fy=simplepdf(yi,0,s2,'gaussian');

fz=pdfinv(yi,fy);

y1=randn(100000,1)*s2;
[fz1,n]=hist(1./y1,(-11:.1:11));

figure,
plot(yi,fz)
hold on
plot(n,fz1/10000,'.'),xlim([-10 10])
title('PDF of the inverse of a Gaussian RV')
text(4,0.30,'Dots are from a random trial')
