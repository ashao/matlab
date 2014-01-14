%
% JTRIADS   Gravity-capillary wave triad interactions.
%
%   See also: 
%   CONTENTS, JARRAY, JMATH, JPOLY, JGRAPH, JSTRINGS, JSTATS, JSIGNAL,
%   JELLIPSE, JCELL, VTOOLS, JOCEANS, JSPHERE, JSATFUN, *JTRIADS*, JPAPERS
%   _________________________________________________________________
%
% Gravity-capillary wave basics  
%   om         - Dispersion relation for gravity-capillary waves.
%   kmin       - Wavenumber of minimum phase speed for gravity-capillary waves.
%   gc_params  - Parameters for gravity-capillary waves.
%
% Resonant triads
%   triadres   - Solves the triad resonance condition given two waves.
%   vtriadres  - Returns resonant wave triads given a sum wavenumber.
%   isres      - Test whether input wavenumbers form a resonant triad.
%   rescoeff   - Interaction coefficients for a resonant wave triad.
%   triadevolve  - McGoldrick's evolution functions for a resonant wave triad.
%
% Triad interaction functions
%   kfun       - Combined interaction functions for wave triads.
%   hfun       - Sea surface height interaction functions for wave triads.
%   dfun       - Velocity potential interaction functions for wave triads.
%   i2ss       - Convert from i=1:4 notation to (s1,s2) notation.
%   ss2i       - Convert from (s1,s2) notation to i=1:4 notation.
%
% Spectra of sets of interacting discrete waves
%   dmspec     - Computes spectrum and bispectrum for discrete wave modes.
%   dmstd      - Computes standard deviation from discrete-mode formulation.
%   dmskew     - Computes skewness from discrete-mode formulation.
%   dmasym     - Computes asymmetry from discrete-mode formulation.
%
% Wavenumber gridding tools
%   wavegrid      - Makes a complex-valued valued (x+iy) grid.
%   iswavegrid    - Tests whether a matrix is in WAVEGRID format. 
%   k2sub         - Convert a complex-valued wavenumber into an I,J subscript pair.
%   sub2k         - Convert an I,J subscript pair into a wavenumber.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2011 J.M. Lilly --- type 'help jlab_license' for details   

help Jtriads

if 0            
             %Gravity-capillary wave basics  
               om         %- Dispersion relation for gravity%-capillary waves.
               kmin       %- Wavenumber of minimum phase speed for gravity%-capillary waves.
               gc_params  %- Parameters for gravity%-capillary waves.
            
             %Resonant triads
               triadres   %- Solves the triad resonance condition given two waves.
               vtriadres  %- Returns resonant wave triads given a sum wavenumber.
               isres      %- Test whether input wavenumbers form a resonant triad.
               rescoeff   %- Interaction coefficients for a resonant wave triad.
               triadevolve  %- McGoldrick's evolution functions for a resonant wave triad.
            
             %Triad interaction functions
               kfun       %- Combined interaction functions for wave triads.
               hfun       %- Sea surface height interaction functions for wave triads.
               dfun       %- Velocity potential interaction functions for wave triads.
               i2ss       %- Convert from i=1:4 notation to (s1,s2) notation.
               ss2i       %- Convert from (s1,s2) notation to i=1:4 notation.
            
             %Spectra of sets of interacting discrete waves
               dmspec     %- Computes spectrum and bispectrum for discrete wave modes.
               dmstd      %- Computes standard deviation from discrete%-mode formulation.
               dmskew     %- Computes skewness from discrete%-mode formulation.
               dmasym     %- Computes asymmetry from discrete%-mode formulation.
            
             %Wavenumber gridding tools
               wavegrid      %- Makes a complex%-valued valued (x+iy) grid.
               iswavegrid    %- Tests whether a matrix is in WAVEGRID format. 
               k2sub         %- Convert a complex%-valued wavenumber into an I,J subscript pair.
               sub2k         %- Convert an I,J subscript pair into a wavenumber.
end