function[a,sigt,sigo]=morsebox(ga,be)
%MORSEBOX  Heisenberg time-frequency box for generalized Morse wavelets.
%
%   A=MORSEBOX(GAMMA,BETA) returns the time-bandwidth product of the
%   a generalized Morse wavelet, that is, the area of its Heisenberg box.
%
%   [A,SIGT,SIGO]=MORSEBOX(GAMMA,BETA) also returns the box width in time 
%   SIGT and in frequency SIGO, with A=SIGT*SIGO. 
%
%   Both SIGT and SIGO are non-dimensionalized with respect to the (radian)
%   peak frequency MORSEFREQ, as in Lilly and Olhede (2009) given below.
%
%   Not for all values of BETA and GAMMA are these standard derivations
%   defined as real-valued quantities.  For locations where these would 
%   be imaginary, their values are set to NAN.
%
%   For details see
%
%       Lilly and Olhede (2009).  Higher-order properties of analytic 
%           wavelets.  IEEE Trans. Sig. Proc., 57 (1), 146--160.
%
%   See also MORSEFREQ, MORSEWAVE.
%
%   'morsebox --t' runs a test.
%   'morsebox --f' generates a sample figure.
%
%   Usage: a=morsebox(ga,be);
%          [a,sigt,sigo]=morsebox(ga,be);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2011 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(ga, '--t')
    morsebox_test,return
end
if strcmp(ga, '--f')
    morsebox_fig,return
end

be(be==1)=1+1e-10;  %Correction for unstable ratios and beta=1

[mo,no]=morsemom(0,ga,be);
[m1,n1]=morsemom(1,ga,be);
[m2,n2]=morsemom(2,ga,be);
[m3,n3]=morsemom(3,ga,be);


[ma,na]=morsemom(0,ga,be-1);
[mb,nb]=morsemom(0,ga,be-1+ga);
[mc,nc]=morsemom(0,ga,be-1+ga/2);


rata=frac(morsea(ga,be),morsea(ga,be-1)).^2;
ratb=frac(morsea(ga,be),morsea(ga,be-1+ga)).^2;
ratc=frac(morsea(ga,be),morsea(ga,be-1+ga/2)).^2;


om = morsefreq(ga,be);

sigo=frac(1,om).*sqrt(frac(n2,no)-frac(n1,no).^2);
index=find(frac(n2,no)-frac(n1,no).^2<=0);
if ~isempty(index),
    sigo(index)=nan;
end
sig2a=  rata.*be.^2.*frac(na,no);
sig2b=  ratb.*ga.^2.*frac(nb,no);
sig2c=2*ratc.*be.*ga.*frac(nc,no);
sigt=real(om.*sqrt(sig2a+sig2b-sig2c));
%Note: take real part because of small numerical noise for 0<beta<1
%Imaginary part is much much smaller than real part

index=find(sig2a+sig2b-sig2c<=0);
if ~isempty(index)
    sigt(index)=nan;
end

sigt(be<=1/2)=inf;

%frac(be,sqrt(2.^(1./ga)));
a=sigt.*sigo;



function[]=morsebox_test
 
t=(-500:500)';
n=0;
ga1=(2:2:10);
be1=(2:2:10);
clear sigmat sigmao
for i=1:length(ga1)
    for j=1:length(be1)
        n=n+1;
        fs=2*pi/20;
        psi=morsewave(length(t),1,ga1(i),be1(j),fs,'bandpass');
        [mut,sigmat(j,i)]=pdfprops(t,abs(psi).^2); 
        sigmat(j,i)=sigmat(j,i).*fs;
    
        fs=2*pi/5;
        psi=morsewave(length(t),1,ga1(i),be1(j),fs,'bandpass');
        [muo,sigmao(j,i)]=pdfprops(2*pi*fftshift(t./length(t)),abs(fft(psi)).^2); 
        sigmao(j,i)=sigmao(j,i)./(fs);
    end
end


[ga,be]=meshgrid(ga1,be1);
[a,sigt,sigo]=morsebox(ga,be);

reporttest('MORSEBOX numerical trials for SIGT',aresame(sigt,sigmat,1e-5))
reporttest('MORSEBOX numerical trials for SIGO',aresame(sigo,sigmao,1e-3))


function[]=morsebox_fig

ga1=(1/3:.1:11);
be1=(1/3:.1:10);

[ga,be]=meshgrid(ga1,be1);
[fm,fe,fi,cf] = morsefreq(ga,be);
a=morsebox(ga,be);

figure
contourf(ga1,be1,a,(.5:.01:.6))
hold on, contour(ga1,be1,a,(.500:.002:.51),'k:')
colormap gray,flipmap,
axis([1 10 1 10])
xtick(1:10),ytick(1:10)
ax=gca;
hc=colorbar;hc=discretecolorbar(hc,[.5 .6],(.5:.01:.6));
axes(hc),hlines((.500:.002:.51),'k:')
axes(ax)

contour(ga1,be1,(fm-fe)./(2*pi),[ 0 0],'k','linewidth',2)
contour(ga1,be1,(fm-fi)./(2*pi),[ 0 0],'k','linewidth',2)
contour(ga1,be1,cf./(2*pi),[ 0 0],'k','linewidth',2)
caxis([.5 .6])
vlines(3,'k--')
plot(ga1,12./ga1,'k')
 
title('Morse Wavelet Area and Transitions')
xlabel('Gamma Parameter')
ylabel('Beta Parameter')
plot([1+sqrt(-1)*0 10+sqrt(-1)*9/2],'k','linewidth',3)
plot([1+sqrt(-1)*0 10+sqrt(-1)*9/2],'w--','linewidth',2)




