function[a]=dmskew(L,K,B,sig)
% DMSKEW   Computes skewness from discrete-mode formulation.
%  
%   MU=DMSKEW(L,K,B1,SIG), where L is a scalar, K is a array of complex-
%   valued wavenumbers, B1 is an array of the same size as K containing
%   the corresponding values of the reduced bispectrum, and SIG is the 
%   standard deviation, returns the skewness.
%  
%   B1 and K may also be of different sizes, but with the size of S
%   "compatible" with that of K. In this case DMASKEW will return a 
%   row vector whose length is NUBSLABS(B1,K), e.g. the number of 
%   time steps.  SIG must then have the same length (or length 1).
%   
%   See also ISCOMPAT, DMSPEC, DMSTD, DMASYM.
%
%   Usage: mu=dmskew(L,k,B1,sig);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002--2006 J.M. Lilly --- type 'help jlab_license' for details    

  
if ~iscompat(B,K)
  error('The size of B must be compatible with that of K.  See ISCOMPAT.')
end

N=numslabs(B,K);
if length(sig)==1
  sig=sig+zeros(N,1);
end

if length(sig)~=N
  error('SIG must be of the same length as NUMSLABS(B1,K).')
end

B=reshape(B,numel(K),N);
sig=row2col(sig);
fact=frac(1,sig.^3*L^2);
a=zeros(N,1);

for i=1:N
   a(i)=sumsum(B(:,i)*fact(i));
end

a=real(a);  %Because of residual imaginary part due to numerical noise
