function[varargout]=makefigs_multivariate(str)
%MAKEFIGS_MULTIVARIATE  Makes figures for Lilly and Olhede (2012), in press at ITSP.
%
%   MAKEFIGS_MULTIVARIATE  Makes all figures for 
%
%                       Lilly & Olhede (2012)
%           "Analysis of Modulated Multivariate Oscillations"
%          In press at IEEE Transactions on Signal Processing
%
%   Type 'makefigs_multivariate' at the matlab prompt to make all figures 
%   for this paper and print them as .eps files into the current directory.
%  
%   Type 'makefigs_multivariate noprint' to supress printing to .eps files.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2011 J.M. Lilly --- type 'help jlab_license' for details


if nargin==0
  str='print';
end

if strcmp(str,'--t')
     return
elseif strcmp(str,'--f')
     makefigs_multivariate('noprint');return
end

%/******************************************************************
load ebasn
use ebasn

len=cellength(lat);
index=find(len>200);
lato=30;
lono=-25;



vindex(id,num,lat,lon,p,t,s,index,1);


clear ga be fs p
for i=1:length(lat)
    if i==25|i==26
        %Special treatment for 36 (148000, now 25) and 37 (149000, now 26), 
        %       which are not getting ridges otherwise
        ga(i)=3;
        be(i)=8;
        fs{i}=morsespace(ga(i),be(i),{0.05,pi},2*pi/100,8);
    else
        ga(i)=3;
        be(i)=3;
        fs{i}=morsespace(ga(i),be(i),{0.05,2*pi/3},2*pi/100,8);
    end
    ppsi(i)=morseprops(ga(i),be(i));
end

%Compute wavelet transforms using generalized Morse wavelets
dt=num{1}(2)-num{1}(1);
clear wx wy cx ir jr xr fr
for i=1:length(lat)
    i
    cx{i}=fillbad(latlon2xy(lat{i},lon{i},lato,lono));   
    %cv=latlon2uv(num{i},lat{i},lon{i});
    %cx{i}=cumsum(cv).*dt*frac(24*3600,100*1000);  %wow that's clever
    [wx{i},wy{i}]=wavetrans(real(cx{i}),imag(cx{i}),{1,ga(i),be(i),fs{i},'bandpass'},'mirror');
    [ir{i},jr{i},xr{i},fr{i},br{i},cr{i}]=ridgewalk(dt,wx{i},wy{i},fs{i},{2*ppsi(i),0,'amp'});  
end
F
%for i=1:length(ir),length(find(isnan(ir{i}))),end

clear xhat fhat xres lathat lonhat latres lonres kap lam the phi fbar bhat chat bbar cbar
for i=1:length(lat)
    %Small overlap of two weak ridges for i=9 and i=23 
    [xhat{i},fhat{i},bhat{i},chat{i}]=ridgemap([length(cx{i}) 2],xr{i},fr{i},br{i},cr{i},ir{i},'collapse');
    [d1{i},d2{i}]=mom2dev(xhat{i},fhat{i},bhat{i},chat{i},2);
    fbar{i}=jointfreq(xhat{i},fhat{i},2);
    d2{i}=frac(1,2)*squared(frac(ppsi(i),vrep(fbar{i},2,2))).*d2{i};
    %[fbar{i},bbar{i},cbar{i}]=jointfreq(xhat{i},fhat{i},bhat{i},chat{i},2);
    [kap{i},lam{i},the{i},phi{i}]=ellparams(xhat{i}(:,1),xhat{i}(:,2));  
    xhat{i}=real(xhat{i}(:,1))+sqrt(-1)*real(xhat{i}(:,2));
    xres{i}=cx{i}-vswap(xhat{i},nan,0);
    [kapres{i},lamres{i},theres{i},phires{i}]=ellparams(anatrans(real(xres{i})),anatrans(imag(xres{i})));
    [kaperr{i},lamerr{i},theerr{i},phierr{i}]=ellparams(d2{i}(:,1),d2{i}(:,2));  
    [latres{i},lonres{i}]=xy2latlon(xres{i},lato,lono);
    %  bbar{i}=bbar{i}./fbar{i};
    %  err{i}=cbar{i}./fbar{i}.^2.*morseprops(ga(i),be(i)).^2/2;  
end



clear cv cvres cvhat
for i=1:length(lat)
    cv{i}=latlon2uv(num{i},lat{i},lon{i});
    cv{i}=fillbad(cv{i});
    cvres{i}=latlon2uv(num{i},latres{i},lonres{i});
    cvres{i}=fillbad(cvres{i});
    cvhat{i}=cv{i}-cvres{i};
    temp=vdiff(real(d2{i}(:,1))+sqrt(-1)*real(d2{i}(:,2)),1)*1000*100/3600/24;
    temp(~isfinite(temp))=0+sqrt(-1)*1e-10;
    cverr{i}=temp;
end
jj=find(cellength(ir)>0);
%\******************************************************************

%/******************************************************************
figure
ax=[-30.5 -19.5 18 36.5];

subplot(141),h=cellplot(lon,lat);linestyle -h h D
hold on,h=cellplot(lon(jj),lat(jj));linestyle -h h k
h1=plot(lon{24},lat{24});linestyle -h h1 5C
h1=plot(lon{24},lat{24});linestyle -h h1 1k


ylabel('Latitude')
xlabel('West Longitude'),latratio(30),axis(ax)
title('Float Trajectories')



subplot(142)
xlabel('West Longitude'),title('Signal Ellipses')

for i=1:length(lat)
    if anyany(~isnan(kap{i}))
        dt=round(frac(2*pi,fbar{i}));
        clear index
        index(1)=dt(find(~isnan(dt),1,'first'));

        while index(end)<length(kap{i})-dt(find(~isnan(dt),1,'last'))
            index(end+1)=index(end)+dt(index(end));
        end
        index=nonnan(index);
        index=index(~isnan(kap{i}(index)));
        if ~isempty(index)
            ar=[frac(360,2*pi*radearth) frac(360,2*pi*radearth*cosd(30))];
            h=ellipseplot(2*kap{i},lam{i},the{i},lonres{i}+sqrt(-1)*latres{i},ar,'index',index);hold on
        end
    end
end
latratio(30),axis(ax)
linestyle k D


subplot(143),h=cellplot(lonres,latres);linestyle -h h D
hold on,h=cellplot(lonres(jj),latres(jj));linestyle -h h k
xlabel('West Longitude'),latratio(30),axis(ax),title('Residuals')
h1=plot(lonres{24},latres{24});linestyle -h h1 5C
h1=plot(lonres{24},latres{24});linestyle -h h1 1k



subplot(144)
xlabel('West Longitude'),title('Bias Ellipses')

for i=1:length(lat)
    if anyany(~isnan(kaperr{i}))
        dt=round(frac(2*pi,fbar{i}));
        clear index
        index(1)=dt(find(~isnan(dt),1,'first'));

        while index(end)<length(kaperr{i})-dt(find(~isnan(dt),1,'last'))
            index(end+1)=index(end)+dt(index(end));
        end
        index=nonnan(index);
        index=index(~isnan(kaperr{i}(index)));
        if ~isempty(index)
            ar=[frac(360,2*pi*radearth*cosd(30)) frac(360,2*pi*radearth)];
            h=ellipseplot(2*kaperr{i},lamerr{i},theerr{i},lonres{i}+sqrt(-1)*latres{i},ar,'index',index);hold on
        end
    end
end
latratio(30),axis(ax)
linestyle k D


for i=1:4
    subplot(1,4,i)
    xtick([-30:2.5:-20])
    xtl=get(gca,'xticklabel');
    set(gca,'xticklabel',xtl(:,2:end))
    %fixlabels([-1,0])
end
letterlabels(4)
packcols(1,4)


orient landscape
fontsize 12 10 10 10
set(gcf,'paperposition',[1 1 11 7])

if strcmp(str,'print')
   print -depsc ebasn-trajectories.eps
end
%\******************************************************************


%/******************************************************************
figure
ci=(0:5:65);
ii=24;
numo=datenum(1986,1,1)-1;
[h,hl]=wavespecplot(num{ii}-numo,cv{ii},2*pi./fs{ii},sqrt(abs(wx{ii}).^2+abs(wy{ii}).^2),1,ci);
linestyle -h hl k k--
axes(h(1)),ylim([-18 18]),ylabel('Current Velocity (cm/s)'),title('Multivariate Ridge Method Example')
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
fontsize 16 14 14 14
set(gcf,'paperposition',[1 1 9 5])
if strcmp(str,'print')
    print -deps ebasn-example.eps
end
%\******************************************************************

