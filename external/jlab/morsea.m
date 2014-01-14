function[a]=morsea(ga,be,str)
%MORSEA  Returns the generalized Morse wavelet amplitude "a".
%
%   A=MORSEA(GAMMA,BETA) returns the generalized Morse wavelet 
%   amplitude, called "A_{BETA,GAMMA}" by Lilly and Olhede (2009).
%
%   By default, A is chosen such that the maximum of the frequency-
%   domain wavelet is equal to 2, the ``bandpass normalization.''
%
%   A=MORSEA(GAMMA,BETA,'energy') instead returns the coefficient
%   giving the wavelet unit energy.  
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006--2011 J.M. Lilly --- type 'help jlab_license' for details

if strcmp(ga,'--t')
      morsea_test;return
end
    
if nargin==2
    str='ban';
end

if strcmp(str(1:3),'ban')
    om=morsefreq(ga,be);          
    a=frac(2,(om.^be).*exp(-om.^ga));
elseif strcmp(str(1:3),'ene')
    a=sqrt(frac(2 *pi*ga.*2.^frac(2*be+1,ga),gamma(frac(2*be+1,ga)))); 
end

function[]=morsea_test

ga1=(2:1:9);
be1=(1:1:10);
[ga,be]=meshgrid(ga1,be1);
om=morsefreq(ga,be);

dom=0.01;
omgrid=permute((0:dom:20)',[3 2 1]);
omgrid=vrep(omgrid,length(ga1),2);
omgrid=vrep(omgrid,length(be1),1);

omgrid=omgrid.*vrep(om,size(omgrid,3),3);
a=morsea(ga,be,'energy');

begrid=vrep(be,size(omgrid,3),3);
gagrid=vrep(ga,size(omgrid,3),3);
agrid=vrep(a,size(omgrid,3),3);

psi=agrid.*omgrid.^begrid.*exp(-omgrid.^gagrid);
psiint=vsum(psi.^2,3).*dom.*om./(2*pi);

reporttest('MORSEA unit energy', allall(abs(psiint-1)<1e-2))

