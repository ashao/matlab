function[a]=conflimit(x,y,alpha,mu)
%CONFLIMIT Computes confidence limits for a probability density function.
%
%   A=CONFLIMIT(X,FX,ALPHA) finds the distance A about the mean such
%   that for the pdf FX distributed over X, MEAN(FX) +/- A encloses
%   ALPHA percent of the area of FX.  X must be uniformly spaced.
%  
%   A=CONFLIMIT(X,FX,ALPHA,X0) finds the distance A about the value
%   X=X0 rather than the mean.  This is useful, for instance, with
%   one-sided probability density functions where one is interested in
%   the distance above zero.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001, 2004 J.M. Lilly --- type 'help jlab_license' for details    
  
  
%   Not currently supported:  
%
%   [A,B]=CONFLIMIT() for an asymmetric probability density FX
%   function gives the distance A below the mean and distance B above
%   the mean, within which interval ALPHA percent of the area of FX is
%   concentrated.   

if strcmp(x,'--t'),return,end
if strcmp(x,'--f')
  conflimit_fig;return
end

if size(x,2)==1
  x=osum(x,0*y(1,:)');
end

vswap(y,nan,0);

if nargin==3
    for i=1:size(y,2)
      mu=pdfprops(x(:,i),y(:,i));
    end
end

for i=1:size(y,2)
    dx=x(2,i)-x(1,i);

    yi=zeros(length(y(:,i)),1);
    [temp,mi]=min(abs(x(:,i)-mu(:,i)));
    i1=mi:length(y(:,i));
    i2=1:mi-1;
    yi(1:length(i1),1)=y(i1,i);
    yi(2:length(i2)+1,1)=yi(2:length(i2)+1,1)+flipud(y(i2,i));
    yi(:,1)=cumsum(yi(:,1))*dx;
    
%     Trapcumsum version is not worth it
%     i1=mi-1:length(y);
%     i2=1:mi+1;
%     yi(1:length(i1),1)=y(i1,i);
%     yi(1:length(i2),1)=yi(1:length(i2),1)+flipud(y(i2,i));
%     yi(:,1)=trapcumsum(yi(:,1),dx);

    ii=min(find(yi>=alpha./100));
    if ~isempty(ii)
        a(i,1)=dx*(ii-1);
    else
        a(i,1)=0;
    end
end

function[]=conflimit_fig

x=(-10:.1:10)';
f=simplepdf(x,0,1,'gaussian');
a=conflimit(x,f,95);
figure,plot(x,f)
vlines([a -a])
title('Unit variance Gaussian, 95% confidence')
