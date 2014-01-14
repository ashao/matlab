function[a]=dmasym(L,K,B,sig,n)
% DMASYM   Computes asymmetry from discrete-mode formulation.
%  
%   A=DMASYM(L,K,B1,SIG,N), where L is a scalar, K is a array of complex-
%   valued wavenumbers, B1 is an array of the same size as K containing
%   the corresponding values of the reduced bispectrum, and SIG is the 
%   standard deviation, returns the asymmetry along direction ARG(N).  
%  
%   B1 and K may also be of different sizes, but with the size of S
%   "compatible" with that of K. In this case DMASYM will return a 
%   row vector whose length is NUBSLABS(B1,K), e.g. the number of 
%   time steps.  SIG must then have the same length (or length 1).
%
%   The asymmetry A has size NUBSLABS(B1,K) x LENGTH(N).
%
%   See also ISCOMPAT, DMSPEC, DMSKEW, DMSTD.
%
%   Usage: A=dmasym(L,k,B1,sig,N);
%  
%   'dmasym --t' runs seme tests. 
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002--2006 J.M. Lilly --- type 'help jlab_license' for details    

  
if strcmp(L, '--t')
  dmasym_test,return
end

if ~iscompat(B,K)
   error('The size of B1 must be compatible with that of K.  See ISCOMPAT.')
end 

N=numslabs(B,K);
K=K(:);  
M=length(K);
B=reshape(B,M,N);
K=vrep(K,N,2);
%B and K are now the same size

sig=col2row(sig);
if length(sig)==N
   sig=vrep(sig,M,1);
elseif ~isscalar(sig)
   error('SIG must be of the same length as NUMSLABS(B1,K).')  
end

for i=(1:length(n))
  fact=sign(cdot(K,n(i))).*frac(3,sig.^3.*L^2);  %3 not 6 b/c I use the sgn not gamma
  a(i,:)=squeeze(vsum(imag(B).*fact,1));
end

a=real(a)';  %strip imaginary part


function[]=dmasym_test

L=1;

K=[3.667*exp(sqrt(-1)*47*pi/360);3.667*exp(-sqrt(-1)*47*pi/360)];
K(3)=K(1)+K(2);
a=[0.59 0.34 0.34]';
phi=[0,0,-pi/2]';

X=a.*exp(sqrt(-1)*phi);
K=[-K;K];
X=[conj(X);X];
[S,B,B1]=dmspec(L,K,X);

bool=allall(abs(imag(B1))-2*(2*pi)^2*abs(a(1)*a(2)*a(3))./(2*pi)^2<1e-10);
bool=bool.*(allall(real(B1)<1e-10));
reporttest('DMSPEC vs mcgtriad for single time slab',bool)

sig=dmstd(L,K,S);
bool=allall(abs(sig-sqrt(2*sum(a.^2)))<1e-10);
reporttest('DMSTD vs mcgtriad for single time slab',bool)

asym=dmasym(L,K,B1,sig,1);
bool=abs(asym-12*sig^(-3)*abs(a(1)*a(2)*a(3))*sin(phi(1)+phi(2)-phi(3)))<1e-10;
reporttest('DMASYM vs mcgtriad for single time slab',bool)

X=[X X];
[S,B,B1]=dmspec(L,K,X);

bool=allall(abs(imag(B1))-2*(2*pi)^2*abs(a(1)*a(2)*a(3))./(2*pi)^2<1e-10);
bool=bool.*(allall(real(B1)<1e-10));
reporttest('DMSPEC vs mcgtriad for multiple time slabs',bool)

sig=dmstd(L,K,S);
bool=allall(abs(sig-sqrt(2*sum(a.^2)))<1e-10);
reporttest('DMSTD vs mcgtriad for multiple time slabs',bool)

asym=dmasym(L,K,B1,sig,1);
asym0=12*sig(1).^(-3)*abs(a(1)*a(2)*a(3))*sin(phi(1)+phi(2)-phi(3)); 
asym1=[asym0 asym0]';
reporttest('DMASYM vs mcgtriad for multiple time slabs',aresame(asym1,asym,1e-10))

asym=dmasym(L,K,B1,sig,rot([0 pi/2 pi]));
asym1=[asym0 asym0; 0 0 ; -asym0 -asym0]';
reporttest('DMASYM vs mcgtriad for multiple slabs and angles',aresame(asym1,asym,1e-10))
