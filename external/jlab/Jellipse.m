%
% JELLIPSE  Bivariate and trivariate signals analysis tools
%
%   See also: 
%   CONTENTS, JARRAY, JMATH, JPOLY, JGRAPH, JSTRINGS, JSTATS, JSIGNAL,
%   *JELLIPSE*, JCELL, VTOOLS, JOCEANS, JSPHERE, JSATFUN, JTRIADS, JPAPERS
%   _________________________________________________________________
%
%     -- see Lilly and Gascard (2006), Nonlin. Proc. Geophys.
%            Lilly and Olhede (2010), IEEE Trans. Sig. Proc.
%            Lilly (2011), Submitted to IEEE Trans. Sig. Proc.
%
% Ellipse analysis of bivariate or trivariate time series
%   ellparams  - Ellipse parameters of a modulated bivariate or trivariate oscillation.
%   ellsig     - Creates a modulated elliptical signal in two or three dimensions.
%   normvect   - Unit normal vector to the ellipse plane in three dimensions.
%
% Bivariate spectral matrices
%   specdiag   - Diagonalize a 2 x 2 spectral matrix.
%   polparams   - Spectral matrix polarization parameters.
%
% Operations on bivariate elliptical time series
%   elldiff    - Ellipse differentiation.
%   ellrad     - Average and instantaneous ellipse radius.
%   ellvel     - Average and instantaneous ellipse velocities.
%   ellband    - Bandwidth of "elliptical" or bivariate analytic signals.
%   ellipseplot  - Plot ellipses.
%
% Other conversions 
%   ecconv     - Converts between eccentricity measures.
%   kl2ab      - Converts ellipse parameters Kappa and Lambda to A and B.
%   ab2kl      - Converts A and B to ellipse parameters Kappa and Lambda.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2011 J.M. Lilly --- type 'help jlab_license' for details      

help Jellipse

if 0
         % Ellipse analysis of bivariate or trivariate time series
             ellparams  %- Ellipse parameters of a modulated bivariate or trivariate oscillation.
             elsig      %- Creates a modulated elliptical signal in two or three dimensions.
             normvect   %- Unit normal vector to the ellipse plane in three dimensions.
          
           %Bivariate spectral matrices
             specdiag   %- Diagonalize a 2 x 2 spectral matrix.
             polparam   %- Spectral matrix polarization parameters.
          
           %Operations on bivariate elliptical time series
             elldiff    %- Ellipse differentiation.
             ellrad     %- Average and instantaneous ellipse radius.
             ellvel     %- Average and instantaneous ellipse velocities.
             ellband    %- Bandwidth of "elliptical" or bivariate analytic signals.
             ellipseplot  %- Plot ellipses.
          
           %Other conversions 
             ecconv     %- Converts between eccentricity measures.
             kl2ab      %- Converts ellipse parameters Kappa and Lambda to A and B.
             ab2kl      %- Converts A and B to ellipse parameters Kappa and Lambda.      
end