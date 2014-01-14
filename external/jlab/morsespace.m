function[fs]=morsespace(ga,be,high,low,D)
%MORSESPACE  Logarithmically-spaced frequencies for generalized Morse wavelets.
%
%   F=MORSESPACE(GAMMA,BETA,HIGH,LOW,D) generates a logarithmically-spaced
%   frequency array for a generalized Morse wavelet.
%
%   The generalized Morse wavelet is specified by parameters GAMMA and BETA.
%   
%   The frequency has frequencies arranged in decending order.  The first
%   (largest) value of F is just larger than HIGH and the smallest is just 
%   larger than LOW.
%
%   HIGH, LOW, and F are all *radian* frequencies.
%
%   The array spacing interval is controlled by the density D.  When D=1, 
%   the peak of one wavelet is approximately located at the half-power 
%   point of the next highest and next lowest wavlet.
%
%   Higher values of D mean more overlap in the frequency domain.  A choice
%   of D=4 is reasonable for general usage.
%   __________________________________________________________________
%
%   Amplitude cutoff
% 
%   F=MORSESPACE(GAMMA,BETA,{ALPHA,HIGH},LOW,D) optionally enforces an
%   amplitude cutoff to determine the highest frequency wavelet. Here
%   ALPHA and HIGH are grouped together in a cell array.
%
%   In this format, the amplitude of the highest frequency wavelet at the
%   Nyquist frequency, divided by its maximum value at the peak frequency, 
%   will not exceed ALPHA.  
%
%   ALPHA is between zero and one.  As ALPHA decreses, the amplitude of
%   the highest-frequency wavelet at the Nyquist frequency decrease.
%
%   The highest frequency in F will then be the minimum of HIGH and the 
%   frequency determined by the cutoff value ALPHA.
%   __________________________________________________________________
%
%   See also MORSEWAVE, WAVETRANS.
%  
%   Usage: f=morsespace(ga,be,high,low,D);
%          f=morsespace(ga,be,{alpha,high},low,D);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009 J.M. Lilly --- type 'help jlab_license' for details

if iscell(high)
    %Recall pi is Nyquist
    high=min(high{2},pi*frac(morsefreq(ga,be),morsehigh(ga,be,high{1})));
end
r=1+frac(1,D*morseprops(ga,be));
N=floor(frac(log(frac(high,low)),log(r)));
fs=high*ones(N+1,1)./r.^[0:N]';


