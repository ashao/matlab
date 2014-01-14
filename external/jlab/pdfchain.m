function[fy]=pdfchain(x,fx,g,yi)
%PDFCHAIN  The "chain rule" for probabilty density functions.
%
%   FY=PDFCHAIN(X,FX,G,Y), where FX is a probability density function
%   defined over values X, and G is some function of X, returns the
%   probility density function of random variable G(X) over values Y.  
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001, 2004 J.M. Lilly --- type 'help jlab_license' for details    
  
warning('off','MATLAB:divideByZero')

vcolon(g,x,fx);

index=find(~isnan(x)&~isnan(fx)&~isnan(g));

x=x(index);
fx=fx(index);
g=g(index);

[g,index]=sort(g);
vindex(x,fx,index,1);

gprime=vdiff(g,1);
%gprime=gprime*ones(size(fx(1,:)));
%x=x*ones(size(fx(1,:)));
%g=g*ones(size(fx(1,:)));

fyx=fx./abs(gprime);
fyx(1)=0;
fyx(end)=0;
index=find(~isfinite(fyx)|fyx==0);
index=index(2:end-1);
if ~isempty(index)
   fyx(index)=fyx(index-1)./2+fyx(index+1)./2;
end
%figure,plot(fyx)
fy=jinterp(g,fyx,yi);

dy=yi(2)-yi(1);
fy=fy./(ones(size(fy(:,1)))*vsum(fy*dy,1));

warning('on','MATLAB:divideByZero')

