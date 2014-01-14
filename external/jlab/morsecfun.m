function[c]=morsecfun(a,ga,be)
%MORSECFUN  Morse wavelet "C"-function.
%
%   C=MORSECFUN(A,GAMMA,BETA) returns the value of the generalized
%   Morse wavelet "C"-parameter in terms of the area of concentration A
%   and the GAMMA and BETA parameters. 
%
%   The input parameters must either be matrices of the same size, or
%   some may be matrices and the others scalars.  
%
%   MORSECFUN uses the formula of Olhede and Walden (2002),
%   "Generalized Morse Wavelets", for the area of "D_{C,BE,GA}" given
%   at the bottom of page 2664.
% 
%   Usage: C = morsecfun(A,ga,be);
%  
%   'morsecfun --t' runs a test
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004 J.M. Lilly --- type 'help jlab_license' for details        
  
if strcmp(a,'--t')
    morsecfun_test;return
end

%/********************************************************
%Sort out sizes
M=1;
N=1;
if length(a(:))~=1
  M=size(a,1);
  N=size(a,2);
end
if length(ga(:))~=1
  M=size(ga,1);
  N=size(ga,2);
end
if length(be(:))~=1
  M=size(be,1);
  N=size(be,2);
end
if length(a(:))==1
  a=a+zeros(M,N);
end
if length(ga(:))==1
  ga=ga+zeros(M,N);
end
if length(be(:))==1
  be=be+zeros(M,N);
end
%\********************************************************
  
r=((2*be)+1)./ga;
fact=2.*pi.*gamma(r+1-(1./ga)).*gamma(r+(1./ga))./(ga.*gamma(r).^2);
c=a./fact+1;

function[b]=morsecfun_test
be=1;
ga=1;
A=[10 150]';

C= morsecfun(A,ga,be);
A2=morsearea(C,ga,be);

tol=1e-10;
b=aresame(A,A2,tol);
reporttest('MORSECFUN inverts MORSEAREA',b);
