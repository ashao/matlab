%
% JSIGNAL  Signal processing, wavelet and spectral analysis
%
%   See also: 
%   CONTENTS, JARRAY, JMATH, JPOLY, JGRAPH, JSTRINGS, JSTATS, *JSIGNAL*,
%   JELLIPSE, JCELL, VTOOLS, JOCEANS, JSPHERE, JSATFUN, JTRIADS, JPAPERS
%   _________________________________________________________________
%
% Multitaper analysis
%   mspec      - Multitaper power spectrum.
%   mtrans     - Multitaper "eigentransform" computation.
%   sleptap    - Slepian tapers.
%   hermfun    - Orthonormal Hermite functions. 
%   hermeig    - Eigenvalues of orthonormal Hermite functions. 
%
% Wavelet analysis
%   wavetrans  - Wavelet transform.
%   morsewave  - Generalized Morse wavelets. [See also below]
%   morlwave   - Morlet wavelet.  [See also below]
%
% Multi-transform polarization analysis
%   msvd       - Singular value decomposition for polarization analysis.
%
% Wavelet ridge analysis
%   ridgewalk  - Extract wavelet transform ridges, including bias estimates.
%   ridgemap   - Map wavelet ridge properties onto original time series.
%
% Frequency, bandwidth,and related quantities
%   instfreq   - Instantaneous frequency, bandwidth, and their generalizations.
%   jointfreq  - Joint instantaneous frequency, bandwidth, and curvature.
%   ellband    - Bandwidth of "elliptical" or bivariate analytic signals.
%   periodindex- Returns time index in increments of instantaneous period.
%
% Assorted other transforms   
%   wigdist    - Wigner distribution (alias-free algorithm).
%   slidetrans - Sliding-window ('moving-window') Fourier transform. 
%   hiltrans   - Hilbert transform.
%   anatrans   - Analytic part of signal.
%
% Morlet wavelet details
%   morlwave   - Morlet wavelet. 
%   morlfreq   - Compute Morlet wavelet carrier frequency given peak frequency.
%
% Generalized Morse details
%   morsewave  - Generalized Morse wavelets of Olhede and Walden (2002). 
%   morseprops - Properties of the demodulated generalized Morse wavelets.
%   morsespace - Logarithmically-spaced frequencies for generalized Morse wavelets.
%   morsefreq  - Minimum and maximum frequencies of Morse wavelets.
%   morsebox   - Heisenberg time-frequency box for generalized Morse wavelets.
%   morsearea  - Time-frequency concentration area of Morse wavelets.
%   morsehigh  - High-frequency cutoff of the generalized Morse wavelets.
%   morsemom   - Frequency-domain moments of generalized Morse wavelets.
%   morsederiv - Frequency-domain derivatives of generalized Morse wavelets.
%   morsexpand - Generalized Morse wavelets via a time-domain Taylor series.
%   morsecfun  - Morse wavelet "C"-function.
%
% Plotting tools
%   wavespecplot - Plot of wavelet spectra together with time series.
%   edgeplot   - Draws limits of edge-effect region on wavelet transform.
%   timelabel  - Put month, day, or hour labels on a time axes.
%
% Low-level functions
%   fourier    -  The Fourier frequencies for a given length time series.
%   ridgequantity - The ``ridge quantity'' associated with a wavelet transform.
%   ridgeinterp - Interpolate quantity values onto ridge locations.
%   powermean   - Power-weighted mean along a specified dimension.
%   morseproj   - Projection coefficient for two generalized Morse wavelets.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2011 J.M. Lilly --- type 'help jlab_license' for details      

%Very low-level functions
%   morsea      - Returns the generalized Morse wavelet amplitude "a".
%   timeseries_boundary - Apply boundary conditions to data before transform.

help Jsignal


if 0 
             %Multitaper analysis
             mspec      %- Multitaper power spectrum.
             mtrans     %- Multitaper "eigentransform" computation.
             sleptap    %- Slepian tapers.
             hermfun    %- Orthonormal Hermite functions. 
             hermeig    %- Eigenvalues of orthonormal Hermite functions. 
          
           %Wavelet analysis
             wavetrans  %- Wavelet transform.
             morsewave  %- Generalized Morse wavelets. [See also below]
             morlwave   %- Morlet wavelet.  [See also below]
          
           %Multi%-transform polarization analysis
             msvd       %- Singular value decomposition for polarization analysis.
          
           %Wavelet ridge analysis
             ridgewalk  %- Extract wavelet transform ridges, including bias estimates.
             ridgemap   %- Map wavelet ridge properties onto original time series.
          
           %Bandwidth and stability
             instfreq   %- Instantaneous frequency, bandwidth, and their generalizations.
             jointfreq  %- Joint instantaneous frequency, bandwidth, and curvature.
             ellband    %- Bandwidth of "elliptical" or bivariate analytic signals.
          
           %Assorted other transforms   
             wigdist    %- Wigner distribution (alias%-free algorithm).
             slidetrans %- Sliding%-window ('moving%-window') Fourier transform. 
             hiltrans   %- Hilbert transform.
             anatrans   %- Analytic part of signal.
          
           %Morlet wavelet details
             morlwave   %- Morlet wavelet. 
             morlfreq   %- Compute Morlet wavelet carrier frequency given peak frequency.
          
           %Generalized Morse details
             morsewave  %- Generalized Morse wavelets of Olhede and Walden (2002). 
             morseprops %- Properties of the demodulated generalized Morse wavelets.
             morsespace %- Logarithmically%-spaced frequencies for generalized Morse wavelets.
             morsefreq  %- Minimum and maximum frequencies of Morse wavelets.
             morsebox   %- Heisenberg time%-frequency box for generalized Morse wavelets.
             morsearea  %- Time%-frequency concentration area of Morse wavelets.
             morsehigh  %- High%-frequency cutoff of the generalized Morse wavelets.
             morsemom   %- Frequency%-domain moments of generalized Morse wavelets.
             morsederiv %- Frequency%-domain derivatives of generalized Morse wavelets.
             morsexpand %- Generalized Morse wavelets via a time%-domain Taylor series.
             morsecfun  %- Morse wavelet "C"%-function.
          
           %Plotting tools
             wavespecplot %- Plot of wavelet spectra together with time series.
             edgeplot   %- Draws limits of edge%-effect region on wavelet transform.
             timelabel  %- Put month, day, or hour labels on a time axes.
          
           %Low-level functions
             fourier    %-  The Fourier frequencies for a given length time series.
             ridgelen      %- Wavelet ridge length expressed as number of full cycles.
             ridgequantity %- The ``ridge quantity'' associated with a wavelet transform.
             ridgeinterp %- Interpolate quantity values onto ridge locations.
             powermean   %- Power%-weighted mean along a specified dimension.
             morseproj   %- Projection coefficient for two generalized Morse wavelets.
end