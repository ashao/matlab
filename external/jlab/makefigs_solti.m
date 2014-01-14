function[]=makefigs_solti(str)
%MAKEFIGS_SOLTI  Makes all figures for Lilly and Lettvin (2004)a.
%
%   MAKEFIGS_SOLTI  Makes all figures for 
%
%                  Lilly & Lettvin (2004)
%    "The 'switch-on' problem for truncated linear operators"
%              Signal Processing #84, 763--784
%
%   Type 'makefigs_solti' at the matlab prompt to make all figures for
%   this paper and print them as .eps files into the current directory.
%  
%   Type 'makefigs_solti noprint' to supress printing to .eps files.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2004 J.M. Lilly --- type 'help jlab_license' for details        
  
%   Tested 10.02.2004 JML with just JLAB and SOLTI directories

if nargin==0
  str='print';
end
  
if strcmp(str,'--t')
     return
elseif strcmp(str,'--f')
     makefigs_solti('noprint');return
end


om0=1;
dt=.025;

dom=0.05;

om1=(-6:dom:6);
t1=2.^(0:.1:4)';

om=ones(size(t1))*om1;
t=t1*ones(size(om1));

d=dk(t,om);
f=fk(t,om);

figure
subplot(121)

ci=(.1:.1:.9);
[h,hc]=contour(om1,log2(t1),d,ci,'k');hold on
ci=-ci;
[h,hc]=contour(om1,log2(t1),d,ci,'k--');


ci=(1:.5:3);
[h,hc]=contour(om1,log2(t1),d,ci,'k');hold on
set(hc,'linewidth',1.5)
ci=-ci;
[h,hc]=contour(om1,log2(t1),d,ci,'k--');hold on
set(hc,'linewidth',1.5)

axis([-3 3 1 4])
grid off

subplot(122)
ci=(.1:.1:.9);
[h,hc]=contour(om1,log2(t1),f,ci,'k');hold on

ci=(.02:.02:.08);
[h,hc]=contour(om1,log2(t1),f,ci,'k:');hold on

ci=(1:.5:3);
[h,hc]=contour(om1,log2(t1),f,ci,'k');hold on
set(hc,'linewidth',1.5)
ci=-ci;
[h,hc]=contour(om1,log2(t1),f,ci,'k--');hold on
set(hc,'linewidth',1.5)
axis([-3 3 1 4])
grid off
set(gca,'yticklabel',[])
subplot(121),title('Continuous Dirichlet kernel D(\omega,t)')
set(gca,'ytick',(0:.5:4))
set(gca,'xtick',(-3:1:2.5))
set(gca,'xticklabelmode','auto')
xlabel('Frequency (\omega)')
ylabel('Log_2(Time)')
fixlabels([0 -1])
h=vlines((0:.5:3),'k');
set(h,'color',0.7*[1 1 1])
subplot(122),title('Continuous Fejer kernel F(\omega,t)')
set(gca,'xtick',(-3:1:3))
set(gca,'ytick',(0:.5:4))
set(gca,'yticklabel',[])
xlabel('Frequency (\omega)')
h=vlines((0:.5:3),'k');
set(h,'color',0.7*[1 1 1])

subplot(121)
pos=get(gca,'position');
pos(2)=.15;
pos(4)=.775;
set(gca,'position',pos);

subplot(122)
pos=get(gca,'position');
pos(2)=.15;
pos(4)=.775;
set(gca,'position',pos);



fontsize jpofigure
set(gcf,'paperposition',[2 2 8 3.5])
if strcmp(str,'print')
  print -depsc solti_figure1.eps
end
%print -depsc kernelsplan.eps 
%!gv kernelsplan.eps &
%\********************************************************



%/********************************************************
dt=.025;
dom=0.05;
om0=1;
om1=(-6:dom:6);
t1=2.^(0:.05:4)';

om=ones(size(t1))*om1;
t=t1*ones(size(om1));

c2=abs(ck(t1,om1,om0)).^2;
fm=fk(t,om-om0);
fp=fk(t,om+om0);
g=frac(pi,2*om0^2).*t.*(fm+fp);
hx=frac(pi,2*om0^2).*hk(t,om,om0);



figure
subplot(131)

ci=(.1:.1:1)*10;
[h,hc]=contour(om1,log2(t1),c2,ci,'k');hold on

ci=(20:20:100);
[h,hc]=contour(om1,log2(t1),c2,ci,'k');hold on
set(hc,'linewidth',2)

axis([-2.5 2.5 1.5 4]),grid off

subplot(132)

ci=(.1:.1:1)*10;
[h,hc]=contour(om1,log2(t1),g,ci,'k');hold on

ci=(20:20:100);
[h,hc]=contour(om1,log2(t1),g,ci,'k');hold on
set(hc,'linewidth',2)

axis([-2.5 2.5 1.5 4]),grid off
set(gca,'yticklabel',[])


subplot(133)

ci=(.1:.1:1)*10;
[h,hc]=contour(om1,log2(t1),hx,ci,'k');hold on
ci=-ci;
[h,hc]=contour(om1,log2(t1),hx,ci,'k--');hold on

ci=(20:10:100);
[h,hc]=contour(om1,log2(t1),hx,ci,'k');hold on
set(hc,'linewidth',2)

axis([-2.5 2.5 1.5 4]),grid off
set(gca,'yticklabel',[])

subplot(131),title('|C(\omega, t; \omega_o)|^2 with \omega_o=1')
set(gca,'ytick',(0:.5:4))
set(gca,'xtick',(-3:1:2.5))
set(gca,'xticklabelmode','auto')
xlabel('Frequency (\omega)')
ylabel('Log_2(Time)')
fixlabels([0 -1])
h=vlines((0:.5:3),'k');
set(h,'color',0.7*[1 1 1])
subplot(132),title('Fejer kernel portion of |C(\omega, t;  \omega_o)|^2')
set(gca,'xtick',(-3:1:3))
set(gca,'ytick',(0:.5:4))
set(gca,'yticklabel',[])
xlabel('Frequency (\omega)')
h=vlines((0:.5:3),'k');
set(h,'color',0.7*[1 1 1])
subplot(133),title('H-kernel portion of |C(\omega, t; \omega_o)|^2')
set(gca,'xtick',(-3:1:3))
set(gca,'ytick',(0:.5:4))
set(gca,'yticklabel',[])
xlabel('Frequency (\omega)')
h=vlines((0:.5:3),'k');
set(h,'color',0.7*[1 1 1])

packcols(1,3)

fontsize jpofigure
set(gcf,'paperposition',[2 2 12 3])
if strcmp(str,'print')
  print -depsc solti_figure3.eps
end
%print -depsc doublekernelsplan.eps 
%!gv doublekernelsplan.eps &

%\********************************************************

t1=2.^[1.5 2.5]';
om=ones(size(t1))*om1;
t=t1*ones(size(om1));
d=dk(t,om);
f=fk(t,om);

%/********************************************************
om1=(0:.5:3);
t1=(0:.01:16)';
om=ones(size(t1))*om1;
t=t1*ones(size(om1));
d=dk(t,om);
f=fk(t,om);


figure,
subplot(121),plot(t1,d,'k'),grid off,axis([0 16 -.2 4])
yoffset 0.25
xlabel('Time') 
ylabel('D(\omega,t) or F(\omega,t)')
title('Dirichlet kernel D(\omega,t) at fixed \omega')
hlines((0:.25:1.5),'k:')
fixlabels([0 -1])
set(gca,'xtick',(0:2:14))
set(gca,'xticklabelmode','auto')

subplot(122),plot(t1,f,'k'),grid off,axis([0 16 -.2 4])
yoffset 0.25
set(gca,'yticklabel',[])
xlabel('Time') 
title('Fejer kernel F(\omega,t) at fixed \omega')
hlines((0:.25:1.5),'k:')
set(gca,'xtick',(0:2:16))

packcols(1,2)
fontsize jpofigure
set(gcf,'paperposition',[2 2 6 3.5])
if strcmp(str,'print')
  print -depsc solti_figure2.eps
end
%print -depsc kernelslices.eps 
%!gv kernelslices.eps &

%\********************************************************


%/********************************************************
dt=.025;
dom=0.05;
om0=1;
om1=(0:.5:2.5);
t1=(0:.01:16)';

om=ones(size(t1))*om1;
t=t1*ones(size(om1));

c2=abs(ck(t1,om1,om0)).^2;
fm=fk(t,om-om0);
fp=fk(t,om+om0);
g=frac(pi,2*om0^2).*t.*(fm+fp);
hx=frac(pi,2*om0^2).*hk(t,om,om0);

figure,
subplot(131),h=plot(t1,c2,'k');
set(h(3),'linewidth',2)
set(h(1),'linestyle','--')
grid off,axis([0 16 -2 75])
yoffset 5
xlabel('Time') 
ylabel('Kernel value (offset)')
title('|C(\omega, t; 1)|^2 with at fixed \omega')
hlines((0:5:30),'k:')
fixlabels([0 -1])
set(gca,'xtick',(0:2:14))
set(gca,'xticklabelmode','auto')

subplot(132),h=plot(t1,g,'k');
set(h(3),'linewidth',2)
set(h(1),'linestyle','--')
grid off,axis([0 16 -2 75])
yoffset 5
set(gca,'yticklabel',[])
xlabel('Time') 
title('Fejer kernel portion')
hlines((0:5:30),'k:')
set(gca,'xtick',(0:2:14))
set(gca,'xticklabelmode','auto')
subplot(133),h=plot(t1,hx,'k');grid off,axis([0 16 -2 75])
set(h(3),'linewidth',2)
set(h(1),'linestyle','--')
yoffset 5
set(gca,'yticklabel',[])
xlabel('Time') 
title('H-kernel portion')
hlines((0:5:30),'k:')
set(gca,'xtick',(0:2:16))

packcols(1,3)

fontsize jpofigure
set(gcf,'paperposition',[2 2 6 3.5])
if strcmp(str,'print')
  print -depsc solti_figure4.eps
end
%print -depsc doublekernelslices.eps 
%!gv doublekernelslices.eps &

%\********************************************************
