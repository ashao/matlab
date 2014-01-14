function[sig]=dmstd(L,K,S)
% DMSTD Computes standard deviation from discrete-mode formulation.
%  
%   SIG=DMSTD(L,K,S), where L is a scalar, K is a array of complex-
%   valued wavenumbers, S is an array of the same size as K containing
%   the corresponding values of the spectrum, returns the standard
%   deviation.
%
%   S and K may also be of different sizes, but with the size of S
%   "compatible" with that of K. In this case DMSTD will return a 
%   row vector whose length is NUBSLABS(S,K), e.g. the number of 
%   time steps.  SIG must then have the same length (or length 1).
%
%   See also ISCOMPAT, DMSPEC, DMASYM, DMSKEW.
%
%   Usage: sig=dmstd(L,k,S);
%
%   'dmstd --t' runs a test. 
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002--2006 J.M. Lilly --- type 'help jlab_license' for details    

  
if strcmp(L, '--t')
  dmstd_test,return
end

if ~iscompat(S,K)
  error('The size of S must be compatible with that of K.  See ISCOMPAT.')
end

S=reshape(S,numel(K),numslabs(S,K));

fact=frac(1,L^2);

%Note two factors of (2 pi)^2 cancel each other
sig=sqrt(sum(fact.*S,1));  

function[]=dmstd_test
reporttest('DMSTD one wave of unit amplitude', aresame(dmstd(1,[3 3],[1 1]),sqrt(2),1e-6));



