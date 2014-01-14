function[varargout]=makefigs_asilomar(str)
%MAKEFIGS_ASILOMAR Makes figures for Lilly and Olhede (2009b).
%
%   MAKEFIGS_ASILOMAR  Makes all figures for 
%
%                       Lilly & Olhede (2009b)
%                   "Wavelet ridge estimation of 
%           jointly modulated multivariate oscillations"
%     In press at Asilomar Conference on Signals, Systems, and Computers.      
%
%   Type 'makefigs_asilomar' at the matlab prompt to make all figures 
%   for this paper and print them as .eps files into the current directory.
%  
%   Type 'makefigs_asilomar noprint' to supress printing to .eps files.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009--2011 J.M. Lilly --- type 'help jlab_license' for details


if nargin==0
  str='print';
end

if strcmp(str,'--t')
     return
elseif strcmp(str,'--f')
     makefigs_asilomar('noprint');return
end

%/******************************************************************
load ebasn
use ebasn

len=cellength(lat);
index=find(len>200);
lato=30;
lono=-25;

id=id(index);num=num(index);
lat=lat(index);lon=lon(index);
p=p(index);t=t(index);s=s(index);

clear ga be fs
for i=1:length(lat)
    if i==25|i==26
        %Special treatment for 36 (148000, now 25) and 37 (149000, now 26), 
        %       which are not getting ridges otherwise
        ga(i)=3;
        be(i)=8;
        fs{i}=morsespace(ga(i),be(i),{0.05,pi},2*pi/10,8);
    else
        ga(i)=3;
        be(i)=3;
        fs{i}=morsespace(ga(i),be(i),{0.05,2*pi/3},2*pi/100,8);
    end
end

%Compute wavelet transforms using generalized Morse wavelets
clear wx wy cx ir jr xr fr
for i=1:length(lat)
    i
    cx{i}=latlon2xy(lat{i},lon{i},lato,lono);
    cx{i}=fillbad(cx{i});
    
    wx{i}=wavetrans(real(cx{i}),{1,ga(i),be(i),fs{i},'bandpass'},'mirror');
    wy{i}=wavetrans(imag(cx{i}),{1,ga(i),be(i),fs{i},'bandpass'},'mirror');
    
    dt=num{i}(2)-num{i}(1);
   % [ir{i},jr{i},xr{i},fr{i}]=ridgewalk(dt,wx{i},wy{i},fs{i},{2*morseprops(ga(i),be(i)),0,'pha'});  
    [ir{i},jr{i},xr{i},fr{i}]=ridgewalk(dt,wx{i},wy{i},fs{i},{2*morseprops(ga(i),be(i)),0,'amp'});  
end

%for i=1:length(ir),length(find(isnan(ir{i}))),end

clear xhat fhat xres lathat lonhat latres lonres kap lam the phi fbar rm vm vg
for i=1:length(lat)
    [xhat{i},fhat{i}]=ridgemap([length(cx{i}) 2],xr{i},fr{i},ir{i});
    
    %Small overlap of two weak ridges for i=9 and i=23 
    if size(xhat{i},3)~=1
        fhat{i}=squeeze(powermean(fhat{i},xhat{i},3));
        xhat{i}=squeeze(powermean(xhat{i},xhat{i},3));
    end
    [kap{i},lam{i},the{i},phi{i}]=ellparams(xhat{i}(:,1),xhat{i}(:,2));  
    fbar{i}=powermean(fhat{i},xhat{i},2);
    xhat{i}=real(xhat{i}(:,1))+sqrt(-1)*real(xhat{i}(:,2));
    xres{i}=cx{i}-vswap(xhat{i},nan,0);
    [latres{i},lonres{i}]=xy2latlon(xres{i},lato,lono);
    rm{i}=ellrad(kap{i},lam{i});
    vm{i}=ellvel(kap{i},lam{i},the{i},phi{i});
    vg{i}=sign(lam{i}).*fbar{i}.*kap{i}.^2./rm{i};
end

clear cv cvres cvhat
for i=1:length(lat)
    cv{i}=latlon2uv(num{i},lat{i},lon{i});
    cv{i}=fillbad(cv{i});
    cvres{i}=latlon2uv(num{i},latres{i},lonres{i});
    cvres{i}=fillbad(cvres{i});
    cvhat{i}=cv{i}-cvres{i};
end

jj=find(cellength(ir)>0);
%\******************************************************************


%/******************************************************************
figure
ax=[-30.5 -19.5 18 36.5];

subplot(131),h=cellplot(lon,lat);linestyle -h h D
hold on,h=cellplot(lon(jj),lat(jj));linestyle -h h k

ylabel('Latitude')
xlabel('West Longitude'),latratio(30),axis(ax)
title('Float Trajectories')

subplot(132)
xlabel('West Longitude'),title('Elliptical Signals')

for i=1:length(lat)
     index=periodindex(num{i}(2)-num{i}(1),fbar{i},1);
     ar=[frac(360,2*pi*radearth*cosd(30)) frac(360,2*pi*radearth) ];
     h=ellipseplot(2*kap{i},lam{i},the{i},lonres{i}+sqrt(-1)*latres{i},ar,'index',index);hold on
end

latratio(30),axis(ax)
linestyle k D

subplot(133),h=cellplot(lonres,latres);linestyle -h h D
hold on,h=cellplot(lonres(jj),latres(jj));linestyle -h h k
xlabel('West Longitude'),latratio(30),axis(ax),title('Residuals')

for i=1:3
    subplot(1,3,i)
    xtick([-30:2.5:-20])
    xtl=get(gca,'xticklabel');
    set(gca,'xticklabel',xtl(:,2:end))
end
letterlabels(1)
packcols(1,3)


orient landscape
fontsize 12 10 10 10
set(gcf,'paperposition',[1 1 9 7])

if strcmp(str,'print')
   print -depsc ebasn-trajectories.eps
end
%\******************************************************************



%/*************************************************    
figure
ci=(0:5:65);
ii=24;
numo=datenum(1986,1,1)-1;
[h,hl]=wavespecplot(num{ii}-numo,cv{ii},2*pi./fs{ii},sqrt(abs(wx{ii}).^2+abs(wy{ii}).^2),1,ci);
linestyle -h hl k k--
axes(h(1)),ylim([-18 18]),ylabel('Current Speed (cm/s)'),title('Bivariate Ridge Method Example')
text(-90,15,'(a)')

axes(h(2)),caxis([0 40]),colormap gray,flipmap,ylim([3.6 60]),hold on
plot(num{ii}-numo,2*pi./fbar{ii},'w','linewidth',4)
plot(num{ii}-numo,2*pi./fbar{ii},'k','linewidth',2)

xlabel('Day of Year 1986'),ylabel('Period in Days')
set(gca,'ytick',2.^(2:.5:5.5))
set(gca,'yticklabel',[' 4';'  ';' 8';'  ';'16';'  ';'32';'  '])
inticks
text(-90,4.5,'(b)')


orient landscape
fontsize 12 10 10 10
set(gcf,'paperposition',[1 1 9 5])
if strcmp(str,'print')
    print -deps ebasn-example.eps
end
%\*************************************************




