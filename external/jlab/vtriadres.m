function[ck2,ck3]=vtriadres(k1,k2i)
%VTRIADRES  Returns resonant wave triads given a sum wavenumber.
%
%   [K1,K2]=VTRIADRES(K3,K) returns the complex-valued wave vectors
%   K1 and K2 satisfying the resonance condition for gravity-
%   capillary wave triads.  Denoting the last dimension of K1 and
%   K2 by 'j', the three "columns" correspond to three resonance
%   conditions:
%  
%         j=1    OM(K1)+OM(K2)=OM(K1+K2)
%         j=2    OM(K1)-OM(K2)=OM(K1-K2)
%         j=3    OM(K2)-OM(K1)=OM(K2-K1)
%
%   All wavenumbers have units of rad cm^-1. K3 must be purely real.
%   K1 and K2 are in complex form, that is, K = Kx + i*Ky.  Values 
%   for gravity and surface tension are specified in GC_PARAMS, 
%   and the dispersion relation is given in OM.  
%
%   See also TRIADRES, TRIADEVOLVE.
%
%   Usage: [k1,k2]=vtriadres(k3,k);
%  
%   'vtriadres --f' makes a sample figure  
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002--2006 J.M. Lilly --- type 'help jlab_license' for details    

if strcmp(k1,'--f')  
   vtriadres_fig, return
end

  
%in the following k1 and k3 are switched
  
thp=triadres(k1,k2i,-1);   %recall this angle will be between K1 and K2
thm=triadres(k1,k2i,1); 

k2i=k2i*ones(1,length(k1));
k1=ones(length(k2i),1)*k1;

ck2=zeros(size(k2i,1),size(k2i,2),3); 
ck3=zeros(size(k2i,1),size(k2i,2),3); 
oms=zeros(size(k2i,1),size(k2i,2),3); 
k1s=zeros(size(k2i,1),size(k2i,2),3); 

%note it's the signs of 1 and 2 that determine thp or thm

%type 1: k1=k3+k2   
ck2(:,:,1)=k2i.*exp(sqrt(-1).*thp);
ck3(:,:,1)=k1-ck2(:,:,1);

%type 2: k1=-k3+k2
ck2(:,:,2)=k2i.*exp(sqrt(-1).*thp);
ck3(:,:,2)=-k1+ck2(:,:,2);

%type 3: k1=k3-k2
ck2(:,:,3)=k2i.*exp(sqrt(-1).*thm);
ck3(:,:,3)=k1+ck2(:,:,3);

oms(:,:,1)= om(ck3(:,:,1))+om(ck2(:,:,1));
oms(:,:,2)=-om(ck3(:,:,2))+om(ck2(:,:,2));
oms(:,:,3)= om(ck3(:,:,3))-om(ck2(:,:,3));

k1s(:,:,1)=real(ck3(:,:,1)+ck2(:,:,1));
k1s(:,:,2)=real(-ck3(:,:,2)+ck2(:,:,2));
k1s(:,:,3)=real(ck3(:,:,3)-ck2(:,:,3));

for i=1:3
  index=find(abs(oms(:,:,i)-om(k1))>.1);
  temp=squeeze(ck2(:,:,i));
  temp(index)=nan;
  ck2(:,:,i)=temp;
  temp=squeeze(ck3(:,:,i));
  temp(index)=nan;
  ck3(:,:,i)=temp;
end


ck2(end+2:2*end+1,:,:)=flipdim(conj(ck2),1);
ck3(end+2:2*end+1,:,:)=flipdim(conj(ck3),1);
ck2((end-1)/2+1,:,:)=nan*ck2((end-1)/2+1,:,:);
ck3((end-1)/2+1,:,:)=nan*ck3((end-1)/2+1,:,:);

ck2=squeeze(ck2);
ck3=squeeze(ck3);


function[]=vtriadres_fig
  

%/********************************************************
k3=[16  8   4     2  1.5  1     1/2     1/4 ].*kmin; 
k2i=(.00000001:.01:10)'*k3;

clear k1 k2
for i=1:length(k3)
    [k1(:,:,i),k2(:,:,i)]=vtriadres(k3(i),k2i(:,i));
end

%k3mat=ones(size(k1(:,1)))*k3+sqrt(-1)*0.000000001;

k1p=squeeze(k1(:,1,:)); %sum interaction
k2p=squeeze(k2(:,1,:)); %sum interaction
k1d=squeeze(k1(:,2,:)); %difference interaction +/-
k2d=squeeze(k2(:,2,:)); %difference interaction +/-

%k1d(find(abs(imag(k1d))<1e-4))=nan;plot(real(k1d(:,8)),'*')
%/********************************************************


%/********************************************************
n=3;
figure
plot(k1p(:,n)./kmin,'k','linewidth',2); hold on
plot(k1d(:,n)./kmin,'k--','linewidth',2); 
plot(k2d(:,n)./kmin,'k--','linewidth',2); 
plot(-k1p(:,n)./kmin,'k--')
plot(-k1d(:,n)./kmin,'k--')
plot(-k2d(:,n)./kmin,'k--')
axis equal
axis([-10 10 -10 10])

%jarrow([0,0],[4,0],0.25,pi/5,2,[1 1 1],'-');
%jarrow([0,0],[4,0],0.2,pi/10,2,[0 0 0],'-');

jj=40;
a(1)=k1p(jj,n)./kmin;
a(2)=k3(n)./kmin-a(1);

jj=200;
a(3)=k1d(jj,n)./kmin;
a(4)=a(3)-k3(n)./kmin;


%plot([0 a(1)],'k-')
%plot([0 a(3)],'k-')
%plot([a(3) a(4) 4+0*sqrt(-1) a(3)],'k-')

plot(a(1),'wo','markersize',8,'markerfacecolor',[1 1 1]*0,'linewidth',1.5)
plot(a(2),'wo','markersize',8,'markerfacecolor',[1 1 1]*0,'linewidth',1.5)
plot(a(3),'w^','markersize',8,'markerfacecolor',[1 1 1]*0,'linewidth',1.5)
plot(a(4),'w^','markersize',8,'markerfacecolor',[1 1 1]*0,'linewidth',1.5)

plot(a(1),'ko','markersize',6,'markerfacecolor',[1 1 1],'linewidth',1)
plot(a(2),'ko','markersize',6,'markerfacecolor',[1 1 1],'linewidth',1)
plot(a(3),'k^','markersize',6,'markerfacecolor',[1 1 1],'linewidth',1)
plot(a(4),'k^','markersize',6,'markerfacecolor',[1 1 1],'linewidth',1)
vlines(0,'k:')
hlines(0,'k:')

%plot([0 a(1) 4+0*sqrt(-1) a(2) 0],'-')

%plot([0 a(3) a(4) 0 a(3) 4+0*sqrt(-1) a(4) 0],'-')
title('Resonance curves for (4,0) \times k_{min}')
xlabel('Eastward wavenumber k/k_{min}')
ylabel('Northward wavenumber l/k_{min}')
set(gca,'xtick',(-10:2:10))
set(gca,'ytick',(-10:2:10))
grid off
axis([-10 10 -10 10]*8/10)

fontsize jpofigure
set(gcf,'paperposition',[2 2 3.5 3.5])
%print -dpsc rescurves.ps 
%!gv rescurves.ps &

