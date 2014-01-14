function [fm,fe,fi,cf] = morsefreq(ga,be)
%MORSEFREQ  Frequency measures for generalized Morse wavelets. [with F. Rekibi]
%
%   [FM,FE,FI]=MORSEFREQ(GAMMA,BETA) calculates three different measures of
%   the frequency of the lowest-order generalized Morse wavelet specified 
%   by parameters GAMMA and BETA.
%
%   FM is the modal or peak, FE is the "energy" frequency, and FI is the 
%   instantaneous frequency at the wavelet center.
%
%   [FM,FE,FI,CF]=MORSEFREQ(GAMMA,BETA) also computes the curvature CF of 
%   the instantaneous frequency at the wavelet center. 
%
%   Note that all frequency quantities here are *radian* as in cos(omega t)
%   and not cyclic as in cos(2 pi f t).
% 
%   The input parameters must either be matrices of the same size,
%   or some may be matrices and the others scalars.   
%
%   For BETA=0, the "wavelet" becomes an analytic lowpass filter, and FM 
%   is not defined in the usual way.  Instead the 1/2 power point of the
%   analytic bandpass filter is returned for FM. 
%
%   For details see
%
%       Lilly and Olhede (2009).  Higher-order properties of analytic 
%           wavelets.  IEEE Trans. Sig. Proc., 57 (1), 146--160.
%
%   See also MORSEBOX, MORSEWAVE.
%
%   Usage: fm = morsefreq(ga,be);
%          [fm,fe,fi] = morsefreq(ga,be);  
%          [fm,fe,fi,cf] = morsefreq(ga,be);  
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2011 J. M. Lilly and F. Rekibi
%                         --- type 'help jlab_license' for details    

%   'morsefreq --f' generates a sample figure, but I'm hiding this since
%   it's uncommented

if strcmp(ga,'--f')
  morsefreq_fig;return
end

fm=frac(be,ga).^frac(1,ga); 
fm(be==0)=(log(2)).^frac(1,ga(be==0)); %Half-power point

%Oher possibilities for zero beta case
%    fm=frac(be+1,ga).^frac(1,ga);%Peak of next wavelet
%    fm=frac(ga-1,ga).^frac(1,ga);%Most rapidly decreasing point


if nargout>1
    fe=frac(1,2.^frac(1,ga)).*frac(gamma(frac(2*be+2,ga)),gamma(frac(2*be+1,ga)));
end

if nargout>2
    fi=frac(gamma(frac(be+2,ga)),gamma(frac(be+1,ga)));
end

if nargout>2
    [m2,n2,k2]=morsemom(2,ga,be);
    [m3,n3,k3]=morsemom(3,ga,be);
 
    cf=-frac(k3,sqrt(k2.^3));
end



if 0
%    C = morsecfun(A,ga,be);

%   MORSEFREQ uses the formula of Olhede and Walden (2002),
%   "Generalized Morse Wavelets", the first formula on page 2665.
%

    r=(2*be+1)./ga;
    coeff=gamma(r+(1./ga)).*2.^(-1./ga);
    fmin=coeff./(2*pi.*gamma(r).*(C+sqrt(C.^2-1)).^(1./ga));
    fmax=coeff./(2*pi.*gamma(r).*(C-sqrt(C.^2-1)).^(1./ga));

    if nargout >4
        fo=frac(1,2).*(fmax+fmin);
    end
    if nargout >5
        fw=frac(1,2).*(fmax-fmin);
    end
end


%function[]=morsefreq_test
 


function[]=morsefreq_fig
ga=[2 3 4 8];
be=[(1:.1:10) (10:.5:100)];

[ga,be]=meshgrid(ga,be);
be(be<frac(ga-1,2))=nan;

[fm,fe,fi,cf]=morsefreq(ga,be);
p=frac(sqrt(be.*ga),pi);


%figure,plot(p,cf./(frac(be+1,ga).^(1./ga)));
%hold on,plot(p,cf./fm,'k');

figure
subplot(221)
plot(p,fe./fm),linestyle k 2k k-- k-. 
title('Energy Frequency / Peak Frequency')
ylim([0.90 1.15]),xlim([0.25 3/2]*2),
hlines(1,'k:'),xtick((0:10)/2);ytick(0.8:0.05:1.2),
xlabel('Duration / Period Ratio'),ylabel('Frequency Ratio')
fixlabels([-1 -2])

subplot(222)
plot(p,fi./fm),linestyle k 2k k-- k-. 
title('Instantaneous / Peak Frequency')
ylim([0.90 1.15]),xlim([0.25 3/2]*2),
hlines(1,'k:'),xtick((0:10)/2);ytick(0.8:0.05:1.2),
xlabel('Duration / Period Ratio'),ylabel('Frequency Ratio')
fixlabels([-1 -2])

subplot(223)
plot(p,cf/2/pi),linestyle k 2k k-- k-. 
title('Instantaneous Frequency Curvature')
ylim([-.15 .15]),xlim([0.25 3/2]*2),
hlines(0,'k:'),xtick((0:10)/2);
ytick(-.15:.05:.15)
xlabel('Duration / Period Ratio'),ylabel('Dimensionless Curvature')
fixlabels([-1 -2])

[a,dt,dom]=morsebox(ga,be);
subplot(224)
plot(p,a),linestyle k 2k k-- k-. 
title('Heisenberg Box Area')
ylim([.5 .55]),xlim([0.25 3/2]*2),
hlines(0,'k:'),xtick((0:10)/2);
ytick(.5:.01:.6)
xlabel('Duration / Period Ratio'),ylabel('Area')
fixlabels([-1 -2])



%/********************************************************************
%Mumerical computation of wavelet properties
%m=(0.3:.05:2);
be=1:.1:6;
f=2*pi*(0.3:.05:2)./2;

fs=2*pi/100;

t=(-10000:10000)';

clear fm fe fi cf p a 
disp('Sorry, this might take a while...')
for i=1:length(f)    
    psi=morlwave(length(t),f(i),fs);
    [fm(i),fe(i),fi(i),cf(i),p(i),a(i)]=morsefreq_numerical(t,psi,fs);
end
%h=gcf;figure,plot(m,fm./fs,'+'),figure(h)

clear fm2 fe2 fi2 cf2 p2 a2 
for i=1:length(be)    
    psi=morsewave(length(t),1,3,be(i),fs,'bandpass');
    [fm2(i),fe2(i),fi2(i),cf2(i),p2(i),a2(i)]=morsefreq_numerical(t,psi,fs);
end

%\********************************************************************

letterlabels(2)

subplot(221)
plot(p2,fe2./fm2,'k.'),plot(p,fe./fm,'k:'),plot(p,fe./fm,'k+')

subplot(222)
plot(p2,fi2./fm2,'k.'),plot(p,fi./fm,'k:'),plot(p,fi./fm,'k+')

subplot(223)
plot(p2,cf2,'k.'),plot(p,real(cf),'k:'),plot(p,real(cf),'k+')

subplot(224)
plot(p2,a2,'k.'),plot(p,a,'k:'),plot(p,a,'k+')


function[fm,fe,fi,cf,p,a]=morsefreq_numerical(t,psi,fs)
warning('off','MATLAB:divideByZero')
om_axis=vshift(2*pi*fftshift(t./length(t)),-1,1);
vswap(om_axis,0,1e-10);

W=fft(fftshift(psi(:,:,1)));
psi=psi.*2./maxmax(abs(W));

[maxpsi,index]=max(abs(fft(psi)));
om=om_axis(index);    
%om=2*pi*fs;
fm=om./(2*pi);

[mut,sigma]=pdfprops(t,psi.*rot(-t.*fm*2*pi)); 
p=real(sigma).*(fm)*2;

dpsi=vdiff(imlog(psi),1);      

omi=dpsi((length(t)+1)/2,1);
fi=omi/(2*pi);

ddom=vdiff(vdiff(dpsi,1),1);

[mut,sigmat]=pdfprops(t,abs(psi).^2); 
[ome,sigmao]=pdfprops(om_axis,abs(fft(psi)).^2);
fe=ome/(2*pi);
cf=frac(ome.^3,ome.^2).*ddom((length(t)+1)/2,1)./sigmao.^3;
a=sigmao.*sigmat;

warning('on','MATLAB:divideByZero')
