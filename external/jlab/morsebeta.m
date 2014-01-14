function[be]=morsebeta(ga,p)
%MORSEBETA  Returns value of beta given gamma and time-bandcenter product.
%
%   BETA=MORSEBETA(GAMMA,P) where P is a time-bandcenter product value
%   returns the value of BETA for which first (K=1) generalized Morse 
%   wavelet has time-bandcenter product equal to P.
%
%   For the first generalized Morse wavelet, the time bandcenter product
%   has the simple expression P=2*SQRT(BETA * GAMMA); 
%
%   The time-bandcenter product P describes the number of oscillations 
%   in the wavelet.  Specifically, if T is the wavelet half-width in time,
%   and F is the wavelet central (cyclic) frequency, then P = 2TF.  The
%   time half-width is considered to be the standard deviation of the
%   demodulated wavelet.
%
%   Usage:  be=morsebeta(ga,p);
%
%   'morsebeta --t' runs a test.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(ga, '--t')
    morsebeta_test,return
end
 
be=frac(p.^2,4*ga);

function[]=morsebeta_test
 
be=(1:1:10);
ga=(2:1:10);
[ga,be]=meshgrid(ga,be);

p=2*sqrt(be.*ga);

be2=morsebeta(ga,p);

tol=1e-10;
reporttest('MORSEBETA',aresame(be2,be,tol))
