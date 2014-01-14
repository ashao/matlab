function[f]=fourier(N)
%FOURIER  The Fourier frequencies for a given length time series.
%
%   F=FOURIER(N) returns the positive Fourier frequencies for a time series
%   of length N.  The units of F are cyclic so that the Nyquist is at 1/2.
%
%   Usage: f=fourier(N);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2011 J.M. Lilly --- type 'help jlab_license' for details
 
if iseven(N)
    f=(0:1./N:1/2-1./N)';
elseif isodd(N)
    f=(0:1./N:1/2)';
end
