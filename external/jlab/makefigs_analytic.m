function[]=makefigs_analytic(str)
%MAKEFIGS_ANALYTIC  Makes figures for Lilly and Olhede (2010b).
%
%   MAKEFIGS_ANALYTIC  Makes all figures for 
%
%                       Lilly & Olhede (2010b)
%               "On the Analytic Wavelet Transform"
%   IEEE Transactions on Information Theory, 56 (8), 4135--4156
%
%   Type 'makefigs_analytic' at the matlab prompt to make all figures for
%   this paper and print them as .eps files into the current directory.
%  
%   Type 'makefigs_analytic noprint' to supress printing to .eps files.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009--2010 J.M. Lilly --- type 'help jlab_license' for details
 
if nargin==0
  str='print';
end

if strcmp(str,'--f')
     makefigs_analytic('noprint');return
end


%/*************************************************
figure
ga=[3 3 3 3 nan 1/3.5 1 6 18 nan];be=[1.5 3 7 81 9./ga(5:end)];
[p,skew,kurt]=morseprops(ga,be);

t=1:5000;t=t-mean(t);
[psi,psif]=vzeros(5000,5,'nan');
for i=1:length(be)
    [psi(:,i),psif(:,i)]=morsewave(5000,1,ga(i),be(i),2*pi/40);
end

textstr{1}='(1.5,3)';
textstr{2}='(3,3)';
textstr{3}='(7,3)';
textstr{4}='(81,3)';
textstr{5}='';
textstr{6}='(31.5,0.29)';
textstr{7}='(9,1)';
textstr{8}='(1.5,6)';
textstr{9}='(0.5,18)';
textstr{10}='';

clear ha
for i=1:10
    ha(i)=subplot(2,5,i);
    if i~=5&&i~=10
        psi(:,i)=psi(:,i)./maxmax(abs(psi(:,i)));
        uvplot(t./p(i).*(2*pi/40),psi(:,i)),hold on,plot(t./p(i).*(2*pi/40),abs(psi(:,i))),linestyle k E-- 2k
        ylim([-1 1.05]),xlim([-4.5 4.5]),
        if i==3
            %text(2,1.5,'Examples of Generalized Morse Wavelets','fontsize',14)
            title('Examples of Generalized Morse Wavelets')
        %    pos=get(get(gca,'title'),'position');
        %    set(get(gca,'title'),'position',pos-[0 .01 0])
        end
        text(-4.2,-.9,textstr{i})
        fixlabels([0 -1]),hlines(0,'k:'),
        if i<4,hlines([-1 1.05],'k'),end
        if i==1,vlines(-4.5,'k'),end
        if i==3,vlines(4.5,'k'),end
        set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
    else
        ylim([0 2.02]),xlim([0 2.5]),hold on
    end
   %title(textstr{i})
end

f=(0:5000-1)./126;
h=subplot(2,5,5);plot(f,psif(:,1:4)),ylim([0 2.02]),xlim([0 2.5]), linestyle E-- 2k k 2E--
ytick(0:.5:2),fixlabels([-1 -1])
boxoff,set(gca,'xtick',[],'ytick',[]),vlines(1,'k:'),
%title('Frequency Domain')
legend('a','b','c','d')

h=subplot(2,5,10);plot(f,psif(:,6:10)),ylim([0 2.02]),xlim([0 2.5]), linestyle E-- 2k k 2E--
ytick(0:.5:2),xlabel('Frequency / Peak Frequency'),fixlabels([-1 -1])
boxoff,set(gca,'xtick',[],'ytick',[]),vlines(1,'k:')
%title('Frequency Domain')
legend('f','g','h','i')

letterlabels(ha,1);

packboth(2,5)

orient landscape
fontsize jpofigure
set(gcf,'paperposition',[1 1 9 4.5])

if strcmp(str,'print')
   %cd_figures
   print -depsc analytic-morsies.eps
end
%\******************************


%/************************************************************
N=20;
ga1=[(1:1:21)];
be1=[(1:1:21)];
[ga,be]=meshgrid(ga1,be1);
[p,skew,kurt]=morseprops(ga,be);
dcell=morsederiv(N,ga,be);

x=zeros([size(ga) N]);
for n=1:length(dcell);
   if iseven(n)
       fact=n/2;
   else
       fact=(n-1)/2;
   end
   x(:,:,n)=dcell{n}./((p.^2).^fact)./factorial(n);
end

figure
for i=1:size(x,1)
    for j=1:size(x,2)
           if ga1(j)<6
                plot(2:size(x,3),squeeze(abs(x(i,j,2:end))),'ko','markersize',5,'markerfacecolor','k'),hold on,ylog
           else 
                plot(2:size(x,3),squeeze(abs(x(i,j,2:end))),'k.','color',[.6 .6 .6]),hold on,ylog
           end
           ylim([10^-6 5]),xlim([1.8 20.5]),xtick([2 3 4 5 10 15 20]),xlog
    end
end

for i=1:size(x,1)
    for j=1:size(x,2)  
           if ga1(j)==3
%                 h=plot(4:size(x,3),squeeze(abs(x(i,j,4:end)))); linestyle -h h 2w
                 h=plot(4:size(x,3),squeeze(abs(x(i,j,4:end)))); linestyle -h h F

           end
    end
end

hlines(1,'k:'),title('Morse Wavelet Decay with \beta>1, 1<\gamma< 6')
xlabel('Derivative Number at Peak Frequency'),ylabel('Normalized Magnitude')

fontsize jpofigure
set(gcf,'paperposition',[1 1 3.5 3.5])

if strcmp(str,'print')
   print -depsc analytic-decay.eps
end
%\********************************************************  



%/*************************************************************************
load ebasn
use ebasn
num=num{33};lat=lat{33};lon=lon{33};p=p{33};t=t{33};
vindex(num,lat,lon,1:549,1);

num=num-datenum(1986,1,1)+1;
fc=corfreq(vmean(lat(:),1))*24/2/pi;

%corfreq(vmean(lat(:),1))*24%=0.76 cycles per day

%vindex(num,lat,lon,1:400,1);
cx=latlon2xy(lat,lon,vmean(lat(:),1),vmean(lon(:),1));
cv=latlon2uv(num,lat,lon);

dt=(num(2)-num(1));
mlat=vmean(lat(:),1);

fmax=2*pi*corfreq(mlat)*frac(dt*24,2); %One cycle per 2 inertial periods = 1 cycle per 2.6 day
fmin=2*pi*corfreq(mlat)*frac(dt*24,40);%One cycle per 40 inertial periods = 1 cycle per 53 days
fs=morsespace(3,3,fmax,fmin,8);

ga=[3 3 3];
be=[1.5 3 7];
p=morseprops(ga, be);
c=squared(2./squared(p));

clear ir jr xr fr ar br cr xra fra bra cra wx
    
for i=1:3
      wx(:,:,i)=wavetrans(real(cx),{1,ga(i),be(i),fs,'bandpass'},'mirror');
    [ir{i},jr{i},xr{i},fr{i},br{i},cr{i}]=ridgewalk(dt,wx(:,:,i),fs,{2*morseprops(ga(i),be(i)),0,'amp'});  
    [xra(:,i),fra(:,i),bra(:,i),cra(:,i)]=ridgemap(length(cx),xr{i},fr{i},br{i},cr{i},ir{i},'collapse');
end


figure
M=4;N=3;
ci=(0:5:65)/2;
  
for i=1:size(xra,2)
    subplot(M,N,i)
    plot(num,vdiff([  real(cx-xra(:,i)) real(xra(:,i)) real(cx)],1))
    yoffset 25,
    linestyle  1.5k  G k
    axis([-95 num(end) -10 70 ])
    if i==1,ylabel('Velocity (cm/s)'),end
    title(['Ridge Analysis with \gamma=' int2str(ga(i)) ', \beta=' num2str(be(i))])
    
    subplot(M,N,i+N)
    contourf(num,2*pi./fs,abs(wx(:,:,i)'),ci),caxis([1 25]),colormap gray,flipmap,flipy,ylog,hold on,nocontours
    plot(num(nonnan(ir{i})),2*pi./fs(nonnan(jr{i})),'w','linewidth',4)
    plot(num(nonnan(ir{i})),2*pi./fs(nonnan(jr{i})),'k','linewidth',2)
    ytick([2 4 8 16 32])
    if i==1,ylabel('Period in Days'),end
    
    subplot(M,N,i+2*N)
    plot(num,fra(:,i),'k'),hold on,plot(num,bra(:,i),'k--')
    axis([-95  num(end) -.5 1.9]),ytick([0:.5:1.5])
    if i==1,ylabel('\omega(t) & \upsilon(t)'),fixlabels([0 -1]),end
    hlines(0,'k:')
    
    subplot(M,N,i+3*N)
    uvplot(num,frac(cra(:,i),fra(:,i).^2)),hold on,ytick([-.2:.1:.2]),fixlabels([0 -1])
    plot(num,frac(bra(:,i),fra(:,i)).^2),%plot(num,[abs(frac(cra(:,i),fra(:,i).^2)) -abs(frac(cra(:,i),fra(:,i).^2))])
    linestyle k G 2k     
    axis([-95  num(end) -.225 .225]),hlines([-1 1]*c(i),'k:')
    if i==1,ylabel('\rho_2(t) & \rho_1^2(t)'),end  
    
    xlabel('Day of Year 1986')
end

for i=1:M*N
    h(i)=subplot(M,N,i);%text(,['(' setstr(i+96) ')' ])
    axes(h(i))
    xtick([-100:100:500])
    set(gca,'xticklabelmode','auto')
end
letterlabels(h,1);


packboth(4,3)
set(gcf,'paperposition',[1 1 9 4.5])
fontsize jpofigure
if strcmp(str,'print')
   print -depsc analytic-transforms.eps
end
%\*************************************************************************

% 
% clear ir2 jr2 xr2 fr2 ar2 br cr2 xra2 fra2 bra2 cra2 
% 
% 
% xrao=xra;
% for i=1:3
%     wx2(:,:,i)=wavetrans(real(xra(:,i)),{1,ga(i),be(i),fs,'bandpass'},'mirror');
%     [ir2{i},jr2{i},xr2{i},fr2{i},br2{i},cr2{i}]=ridgewalk(dt,wx2(:,:,i),fs,{2*morseprops(ga(i),be(i)),0,'amp'});  
%     [xra2(:,i),fra2(:,i),bra2(:,i),cra2(:,i)]=ridgemap(length(cx),xr2{i},fr2{i},br2{i},cr2{i},ir2{i},'collapse');
% end
% %xra3=xra2./(1+1/2*vrep(p.^2,length(xra),1).*cra./fra.^2);
% xra=xra2;
% fra=fra2;
% bra=bra2;
% cx=xrao;
% figure
% M=2;N=3;
%   
% vsum(abs(xra-xra2).^2,1)./vsum(abs(xra).^2,1)
% %0.0389    0.0475    0.0849
% 
% 
% % vsum(abs(real(xra-xra2)).^2,1)./vsum(abs(real(xra)).^2,1)
% 
% 
% for i=1:size(xra,2)
%    
%     subplot(M,N,i)
% 
%   
%     
%     subplot(M,N,i+N)
% 
%     uvplot(num,frac(cra(:,i),fra(:,i).^2)./c(i)),hold on,
%     plot(num,frac(bra(:,i),fra(:,i)).^2./c(i)),%plot(num,[abs(frac(cra(:,i),fra(:,i).^2)) -abs(frac(cra(:,i),fra(:,i).^2))])
%     linestyle k G 2k     
%     axis([-95  num(end) -4.9 4.9]),hlines([-1 1],'k:')
%     if i==1,ylabel('\rho_2(t) & \rho_1^2(t) x  P_\psi^4/4'),end  
%     xtick([0:100:400])
%     
%     %subplot(M,N,i+2*N)
%     %plot(num,vdiff([real(xra(:,i))-real(xra2(:,i)) real(xra(:,i))-real(xra3(:,i))],1))
%     %linestyle k G 1.5k
%     %axis([-95 num(end) -2 2 ])
%     %if i==1,ylabel('Velocity (cm/s)'),end
%     %xtick([0:100:400])
%     %title(['Ridge Analysis with \gamma=' int2str(ga(i)) ', \beta='
%     %num2str(be(i))]
%     xlabel('Day of Year 1986')
% 
% end
% 
% 
% for i=1:M*N
%     h(i)=subplot(M,N,i);%text(,['(' setstr(i+96) ')' ])
% end
% letterlabels(h,1);
% 
% 
% packboth(4,3)
% set(gcf,'paperposition',[1 1 9 3])
% fontsize jpofigure
% if strcmp(str,'print')
%    print -depsc analytic-transforms2.eps
% end
% 
% 
% 
% 
% for i=1:M*N
%     h(i)=subplot(M,N,i);%text(,['(' setstr(i+96) ')' ])
% end
% letterlabels(h,1);
% 
% 
% packboth(4,3)
% set(gcf,'paperposition',[1 1 9 6])
% fontsize jpofigure
% if strcmp(str,'print')
%    print -depsc analytic-transforms.eps
% end

 
    
%     
%  
% %/*************************************************
% load ebasn
% 
% use ebasn
% num=num{33};lat=lat{33};lon=lon{33};p=p{33};t=t{33};
% vindex(num,lat,lon,1:549,1);
% 
% num=num-datenum(1986,1,1)+1;
% fc=corfreq(vmean(lat(:),1))*24/2/pi;
% 
% %corfreq(vmean(lat(:),1))*24%=0.76 cycles per day
% 
% %vindex(num,lat,lon,1:400,1);
% cx=latlon2xy(lat,lon,vmean(lat(:),1),vmean(lon(:),1));
% cv=latlon2uv(num,lat,lon);
% 
% dt=(num(2)-num(1));
% mlat=vmean(lat(:),1);
% 
% fmax=2*pi*corfreq(mlat)*frac(dt*24,2); %One cycle per 2 inertial periods = 1 cycle per 2.6 day
% fmin=2*pi*corfreq(mlat)*frac(dt*24,40);%One cycle per 40 inertial periods = 1 cycle per 53 days
% fs=morsespace(3,3,fmax,fmin,8);
% %fs=morsespace(ga,be,{.2,fmax},fmin,8);
% 
% ga=[3 3 3];
% be=[1.5 3 7];
% p=morseprops(ga, be);
% clear ir jr xr fr ar br cr xra fra bra cra wx
%     
% for i=1:3
%     wx(:,:,i)=wavetrans(real(cx),{1,ga(i),be(i),fs,'bandpass'},'mirror');
%     [ir{i},jr{i},xr{i},fr{i},br{i},cr{i}]=ridgewalk(dt,wx(:,:,i),fs,{2*morseprops(ga(i),be(i)),0,'amp'});  
%     [xra(:,i),fra(:,i),bra(:,i),cra(:,i)]=ridgemap(length(cx),xr{i},fr{i},br{i},cr{i},ir{i},'collapse');
% end
% 
% clear ir2 jr2 xr2 fr2 ar2 br2 cr2 xra2 fra2 bra2 cra2 wx
% 
% for i=1:3
%     wx(:,:,i)=wavetrans(real(xra(:,i)),{1,ga(i),be(i),fs,'bandpass'},'mirror');
%     [ir2{i},jr2{i},xr2{i},fr2{i},br2{i},cr2{i}]=ridgewalk(dt,wx(:,:,i),fs,{2*morseprops(ga(i),be(i)),0,'amp'});  
%     [xra2(:,i),fra2(:,i),bra2(:,i),cra2(:,i)]=ridgemap(length(cx),xr2{i},fr2{i},br2{i},cr2{i},ir2{i},'collapse');
% end
% 
% 
% 
% vmedian(abs((xra-xra2)./vswap(xra,0,nan)).^2,1)
% %ans =   0.0235    0.0139    0.0215
% vmean(abs((xra-xra2)./vswap(xra,0,nan)).^2,1)
% %ans = 0.0404    0.0355    0.0565
% 
% figure
% M=4;N=3;
% ci=(0:5:65)/2;
% 
% c=squared(2./squared(p));
% 
%   
% for i=1:size(xra,2)
%     subplot(M,N,i)
%     plot(num,vdiff([real(cx) real(xra(:,i)) real(cx-xra(:,i))],1))
%     yoffset 25,
%     linestyle k G 1.5k
%     axis([-95 num(end) -10 70 ])
%     if i==1,ylabel('Velocity (cm/s)'),end
%     xtick([0:100:400])
%     title(['Ridge Analysis with \gamma=' int2str(ga(i)) ', \beta=' num2str(be(i))])
%     
%     subplot(M,N,i+N)
%     contourf(num,2*pi./fs,abs(wx(:,:,i)'),ci),caxis([1 25]),colormap gray,flipmap,flipy,ylog,hold on,nocontours
%     plot(num(nonnan(ir{i})),2*pi./fs(nonnan(jr{i})),'w','linewidth',4)
%     plot(num(nonnan(ir{i})),2*pi./fs(nonnan(jr{i})),'k','linewidth',2)
%     ytick([2 4 8 16 32])
%     if i==1,ylabel('Period in Days'),end
%     xtick([0:100:400])
%     
%     subplot(M,N,i+2*N)
%     plot(num,fra(:,i),'k'),hold on,plot(num,bra(:,i),'k--')
%     axis([-95  num(end) -.5 1.9]),ytick([0:.5:1.5])
%     if i==1,ylabel('Frequency & Bandwidth'),fixlabels([0 -1]),end
%     hlines(0,'k:')
%     xtick([0:100:400])
%     
%     subplot(M,N,i+3*N)
% 
%     uvplot(num,frac(cra(:,i),fra(:,i).^2)./c(i)),hold on,
%     plot(num,frac(bra(:,i),fra(:,i)).^2./c(i)),%plot(num,[abs(frac(cra(:,i),fra(:,i).^2)) -abs(frac(cra(:,i),fra(:,i).^2))])
%     linestyle k G 2k     
%     axis([-95  num(end) -4.9 4.9]),hlines([-1 1],'k:')
%     if i==1,ylabel('\rho_2(t) & \rho_1^2(t) x  P_\psi^4/4'),end  
%     xlabel('Day of Year 1986')
%     xtick([0:100:400])
% end
% 
% 
% for i=1:M*N
%     h(i)=subplot(M,N,i);%text(,['(' setstr(i+96) ')' ])
% end
% letterlabels(h,1);
% 
% 
% packboth(4,3)
% set(gcf,'paperposition',[1 1 9 6])
% fontsize jpofigure
% if strcmp(str,'print')
%    print -depsc analytic-transforms.eps
% end
% 
% 
% for i=1:3
%         %vmean(abs(frac(cra(:,i),fra(:,i).^2)),1),vmedian(abs(frac(cra(:,i),fra(:,i).^2)),1)
%         vmean(abs(frac(cra(:,i),fra(:,i).^2))./c(i),1),vmedian(abs(frac(cra(:,i),fra(:,i).^2))./c(i),1)
% end




% 
% %/*************************************************
% load ebasn
% 
% use ebasn
% num=num{33};lat=lat{33};lon=lon{33};p=p{33};t=t{33};
% vindex(num,lat,lon,1:549,1);
% 
% num=num-datenum(1986,1,1)+1;
% fc=corfreq(vmean(lat(:),1))*24/2/pi;
% 
% %corfreq(vmean(lat(:),1))*24%=0.76 cycles per day
% 
% %vindex(num,lat,lon,1:400,1);
% cx=latlon2xy(lat,lon,vmean(lat(:),1),vmean(lon(:),1));
% cv=latlon2uv(num,lat,lon);
% 
% dt=(num(2)-num(1));
% mlat=vmean(lat(:),1);
% 
% fmax=2*pi*corfreq(mlat)*frac(dt*24,2); %One cycle per 2 inertial periods = 1 cycle per 2.6 day
% fmin=2*pi*corfreq(mlat)*frac(dt*24,40);%One cycle per 40 inertial periods = 1 cycle per 53 days
% fs=morsespace(3,3,fmax,fmin,8);
% 
% be=[1:0.25:10];
% ga=3+0*be;
% 
% p=morseprops(ga, be);
% c=squared(2./squared(p));
% 
% clear ir jr xr fr ar br cr xra fra bra cra wx
%     
% for i=1:length(be)
%     i
%     wx(:,:,i)=wavetrans(real(cx),{1,ga(i),be(i),fs,'bandpass'},'mirror');
%     [ir{i},jr{i},xr{i},fr{i},br{i},cr{i}]=ridgewalk(dt,wx(:,:,i),fs,{2*morseprops(ga(i),be(i)),0,'amp'});  
%     [xra(:,i),fra(:,i),bra(:,i),cra(:,i)]=ridgemap(length(cx),xr{i},fr{i},br{i},cr{i},ir{i},'collapse');
% end
% 
% rq=ridgequantity(wx(:,:,3),fs,'amp');
% [omx,bx,cx]=instfreq(wx);
% 
% for i=1:length(be)
%     i
%     [xr{i},fr{i},br{i},cr{i}]=ridgeinterp(wx(:,:,3),fs,rq,ir{3},jr{3},wx(:,:,i),omx(:,:,i),bx(:,:,i),cx(:,:,i));
%     [xra(:,i),fra(:,i),bra(:,i),cra(:,i)]=ridgemap(length(cx),xr{i},fr{i},br{i},cr{i},ir{3},'collapse');
% 
% end
% 
% %    Usage:  xi=ridgeinterp(w,fs,rq,ir,jr,x);
% %pcolor(abs(cra')./abs(fra').^2.*vrep((p').^4,length(xra),2)),flipy,shading flat

