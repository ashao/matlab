function[S,B,B1]=dmspec(L,K,X)
% DMSPEC  Computes spectrum and bispectrum for discrete wave modes.
%  
%   S=DMSPEC(L,K,X), where L is a scalar, K is a array of complex-
%   valued wavenumbers, X is an array of the same size as K containing
%   the corresponding Fourier coefficients, returns the discrete-mode
%   spectrum S, an array of the same size as X and K.
%
%   [S,B,B1]=DMSPEC(L,K,X) optionally returns the discrete-mode 
%   bispectrum B and the 'reduced' bispectrum B1.  
%
%   See also ISCOMPAT, DMSTD, DMSKEW, DMASYM.
%
%   Usage: S=dmspec(L,k,X);
%
%   'dmspec --t' runs a test. 
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002--2006 J.M. Lilly --- type 'help jlab_license' for details    

if strcmp(L, '--t')
  dmspec_test, return
end


tol=1e-10;
  
S=frac(1,L^2)*abs(X.^2);

if nargout>1
  
  if ~iscompat(X,K)
    error('The size of X must be compatible with the size of K.')
  end
      
  sizeK=size(K);  
  K=K(:);
  N1=length(K);
  
  %resize X to match K
  sizeX=size(X);
  X=X(:);
  N2=length(X)./N1;
  X=reshape(X,[N1 N2]);
  
  K3=osum(K,K);
  [ii,jj]=lookup(K,K3(:),tol);

  %extract those pairs for which K3 is included in K.
  %note I'm using a numerical tolerance 
  
  [i1,i2]=ind2sub(size(K3),jj);
 
  B=zeros(N1*N1,N2); 
  index=sub2ind(size(K3),i1,i2);
  B(index,:)=frac(1,L^2).*X(i1,:).*X(i2,:).*conj(X(ii,:));
  B=reshape(B,N1,N1,N2);
  %for i=1:length(i1)
  %B(i1(i),i2(i),:)=frac(1,L^2).*X(i1(i),:).*X(i2(i),:).*conj(X(ii(i),:));
  %end
end

if nargout==3
  B1=frac(1,L^2)*sum(B,1);
  B1=permute(B1,[2 3 1]);
end


function[]=dmspec_test
%/********************************************************
k=[3.667*exp(sqrt(-1)*47*pi/360);3.667*exp(-sqrt(-1)*47*pi/360)];
k(3)=k(1)+k(2);
ahat=[.687 .487]/10/2;%convert to cm
mu=[0 0 -pi/2];  
dmu=mu(1)+mu(2)-mu(3); %dmu is +pi/2

[B12,B13,B23]=rescoeff(k(1),k(2),k(3));
B12=B12/2;B13=B13/2;B23=B23/2;

m=B23.*ahat(2).^2./(B13.*ahat(1).^2);

[g,T]=gc_params;
s3=sqrt(g.*abs(k(3))+T.*abs(k(3)).^3);
p = 2*pi/s3;

t=p*(0:.1:32)';
a=triadevolve(ahat,k,t);

clear sig asym
sig(:,1)=sqrt(2)*(a(:,1).^2+a(:,2).^2+a(:,3).^2).^(1/2);
asym(:,1)=12.*sin(dmu).*prod(a,2)./sig(:,1).^3;

for i=1:3
  dphi(:,i)=mu(i)+pi/2*(1+sign(-a(:,i)));
end
%\********************************************************

%/***************************************************
%Compute bispec, asymmetry for mcg
clear MB1
L=1;
k1=[k;-flipud(k)];
for i=1:length(a)
  X=conj(abs(a(i,:)).*exp(sqrt(-1).*dphi(i,:)))';
  X=[X;conj(flipud(X))]; 
  [Stemp,Btemp,B1temp]=dmspec(L,k1,X);
  sig(i,2)=dmstd(L,k1,Stemp);
  asym(i,2)=dmasym(L,k1,B1temp,sig(i,2),1);
  MB1(i,:)=conj(B1temp');
  clear Stemp Btemp B1temp
end
tol=1e-10;
b1=maxmax(abs(sig(:,1)-sig(:,2)))<tol;
b2=maxmax(abs(asym(:,1)-asym(:,2)))<tol;
reporttest('DMSPEC variance vs. direct expression for M65 triad',b1)
reporttest('DMSPEC asymmetry vs. direct expression for M65 triad',b2)
%\***************************************************
  

