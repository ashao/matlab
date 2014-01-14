function[y]=slidetrans(x,w,fs,str)
%SLIDETRANS  Sliding-window ('moving-window') Fourier transform.
%   
%   Y=SLIDETRANS(X,PSI,F) computes the sliding-window Fourier transform of
%   the signal X using window PSI at frequencies F. X and W are column
%   vectors.  Y is a matrix of size LENGTH(X) by LENGTH(FS).
%
%   The frequencies F are cyclic frequencies as in COS(2 PI F T).   
%
%   PSI may also be a matrix, in which case Y is a 3-D array of size 
%   LENGTH(X) by LENGTH(FS) by SIZE(W,2).
%
%   SLIDETRANS(...,STR) sets the endpoint boundary conditions, with STR
%   equal to 'mirror', 'periodic', or 'zeros'.  The default is 'periodic'. 
%   For more details, see the boundary condition discussion in WAVETRANS.
%
%   A good choice for window is the first "Slepian" taper; see SLEPTAP.
%
%   SLIDETRANS follows the same normalization as WAVETRANS when the input
%   signal X is complex-valued, so see that function for details.
%   _________________________________________________________________
%
%   Definition
%
%   The sliding window transform is defined as
%
%      y(t,f) = int psi^*(u-t) exp[-2 pi i f (u-t)] x(t) du.
%
%   When a sinusoid is transformed at its own frequency, the rate of change
%   of the phase of the sliding window transform recovers the frequency of
%   the sinusoid.  This is the same as with the wavelet transform.
%
%   This differs from the definition of Mallat (1999), p. 69, by a unit
%   amplitude factor of exp[2 pi i t].  In Mallat's defintion, the
%   transform of a sinusoid at the sinusoid's frequency has constant phase.
%   Both definitions have the same modulus.
%   _________________________________________________________________
%
%   'slidetrans --t' runs a test.
%   'slidetrans --f' generates a sample figure.
%
%   Usage: y=slidetrans(x,psi,f);
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2004--2011 J.M. Lilly --- type 'help jlab_license' for details  

if strcmp(x,'--f')
  slidetrans_fig;return
end
if strcmp(x,'--t')
  slidetrans_test;return
end

if nargin<4
    str='periodic';
end

M=size(x,1);
N=size(x,2);

fs=fs(:);
Mf=length(fs);
Mw=size(w,1);
K=size(w,2);

t=(0:M-1)';
tw=(0:Mw-1)';
tw=tw-mean(tw);


%Make a wavelet
w=vrep(permute(w,[1 3 2]),Mf,2);
phasor=vrep(rot(oprod(2*pi*tw,fs)),K,3);
psi=w.*phasor;

%psi=zeros(Mw,Mf,K);
%for j=1:Mf
%  for k=1:K
%    psi(:,j,k)=w(:,k).*rot(2*pi*tw*fs(j));
%  end
%end

y=wavetrans(x,psi,str);

%for j=1:Mf
%  for k=1:K
   % y(:,j,k)=y(:,j,k).*rot(-2*pi*t*fs(j)) ;
%  end
%end
  
function[]=slidetrans_test

N=501;
w=hermfun((-N:N)'./(N/4),0);
w=w./sum(w);
M=3001;
t=(0:M-1)';
x=rot(2*pi*t./70);
y=slidetrans(x,w,1./70);
om=vdiff(unwrap(angle(y)),1);
bool(1)=aresame(om(N:M-N),2*pi/70+0*(N:M-N)',1e-4);
bool(2)=aresame(abs(y(N:M-N)).^2,1/2+0*(N:M-N)',1e-5);
reporttest('SLIDETRANS complex sinusoid',all(bool))

function[]=slidetrans_fig
  
M=3000;
t=(0:M-1)';%t=t-M/2;
N=500;
w=hermfun((-N:N)'./(N/4),0);
w=w./sqrt(w'*w);
fs=(1:30)./1000;
clear x
x(1:M/2,1)=sin(2*pi*t(1:M/2)./70/3);
x(M/2:M,1)=sin(2*pi*t(M/2:M)./70);
y=slidetrans(x,w,fs,'zeros');
%y=slidetrans(x,w,2*pi/70/3);
h=wavespecplot(t,x,1./fs,abs(y),1/2);
hlines(70*3),hlines(70)
% 
% 
% x=testseries(6);
% N=2000;
% w=hermfun([-N:N]'./(N/4),4);
% for i=1:size(w,2)
%   w(:,i)=w(:,i)./sqrt(w(:,i)'*w(:,i));
% end
% fs=(0:0.5:30]./length(x);
% y=slidetrans(x,w,fs);
% t=(0:length(x)-1)';
% jimage(t,fs,abs(y)'),shading interp,flipy
% jimage(t,fs,mean(abs(y(:,:,1:4)),3)'),shading interp,flipy
% contourf(t,fs,abs(y)',20),nocontours,flipy
