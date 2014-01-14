function[fz,mu,sigma]=pdfadd(xi,yi,fx,fy,zi)
%PDFADD  Probability distribution from adding two random variables.
%  
%   YN=PDFADD(XI,YI,FX,FY,ZI), given two probability distribution
%   functions FX and FY defined over values XI and YI respectively,
%   returns the pdf FZ corresponding to Z=X+Y over values ZI.
%
%   'pdfadd --f' generates a sample figure
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001, 2004 J.M. Lilly --- type 'help jlab_license' for details    
 
if strcmp(xi,'--f')
    pdfadd_fig;
    return
end

dx=xi(2)-xi(1);
dy=yi(2)-yi(1);
if abs(dx-dy)>1e-6*dx
  error('Sorry, dx and dy must be the same')
end

vcolon(xi,yi,fx,fy,zi);

fz=conv(fx,fy)';
z=(1:length(fz))'*dx;
z=z-z(1)+xi(1)+yi(1);

fz=jinterp(z,fz,zi);
dz=zi(2)-zi(1);
fz=fz./vsum(fz(:)*dz,1);

[mu,sigma]=pdfprops(zi,fz);

  
function[]=pdfadd_fig

dx=0.1;
dy=0.1;
s1=1;
s2=2;
xi=(-10:dx:10)';
xi=(30:dx:50)';
yi=(-10:dy:10)';
zi=(0:dy:60)';

fx=exp(-(xi-43).^2./s1.^2);
%fx(round(end/2):end)=0;
fx=fx./vsum(fx*dx,1);

fy=exp(-yi.^2./s2.^2);fy=fy./vsum(fy*dy,1);

fz0=conv(fx,fy(1:2:end));
fz0=fz0./vsum(fz0*dx,1);

%size(fz0)
fz=pdfadd(xi,yi,fx,fy,zi);
%size(fz)

figure,plot(xi,fx),hold on,plot(yi,fy),plot(zi,fz)
linestyle default
title('RV with green pdf plus RV with red pdf equals RV with blue pdf')
%plot(xi,fz0,'r')

mu=pdfprops(xi,fx);
vlines(mu)



