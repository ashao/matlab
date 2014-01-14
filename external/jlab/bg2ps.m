function[p,s,k]=bg2ps(be,ga)
%BG2PS  Convert Morse wavelet beta and gamma 
%
%   BG2PS
%
%   'bg2ps --t' runs a test.
%
%   Usage: [p,skew,kurt]=bg2ps(be,ga);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(be, '--t')
    bg2ps_test,return
end
p=sqrt(be.*ga/2);
s=frac(ga-3,sqrt(2)*p);
k=3-s-frac(1,p.^2);

function[]=bg2ps_test
 
%reporttest('BG2PS',aresame())
