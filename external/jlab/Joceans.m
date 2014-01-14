%
% JOCEANS   Oceanography-specific functions
%
%   See also: 
%   CONTENTS, JARRAY, JMATH, JPOLY, JGRAPH, JSTRINGS, JSTATS, JSIGNAL,
%   JELLIPSE, JCELL, VTOOLS, *JOCEANS*, JSPHERE, JSATFUN, JTRIADS, JPAPERS
%   _________________________________________________________________
%
% Date conversion 
%   monthfrac  - Convert date from 'datenum' format to 'month.fraction'.
%   yearfrac   - Convert date from 'datenum' format to 'year.fraction'.
%   yf2num     - Converts date in 'year.fraction' format to  'datenum' format.
%   mjd2num    - Converts Modified Julian Dates to 'datenum' format.
%
% Latitude and longitude
%   corfreq    - Coriolis frequency in cycles per hour.
%   spheredist - Computes great circle distances on a sphere.
%   inregion   - Tests whether lat/lon points lie within a specified box.
%   latlon2xy  - Converts latitude and longitude into local Cartesian coordinates.
%   xy2latlon  - Converts local Cartesian coordinates into latitude and longitude.
%   latlon2uv  - Converts latitude and logitude to velocity.
%   deg360     - Converts degrees to the range [0, 360].
%   deg180     - Converts degrees to the range [-180, 180].
%   lonshift   - Shifts longitude origin for plotting purposes.
%   [-- see also JSPHERE]
%
% Plotting functions
%   stickvect  - Plots "stick vectors" for multicomponent velocity time series.
%   denscont   - Density contour overlay for oceanographic T/S plots.
%   uvplot     - Plot u and v components of velocity on the same axis.
%   provec     - Generate progressive vector diagrams (simple and fancy).
%   hodograph  - Generate hodograph plots (simple and fancy).
%   latratio   - Set plot aspect ratio for latitude / longitude plot.
%   timelabel  - Put month, day, or hour labels on a time axes.
%
% Eddy modelling
%   rankineeddy  - Velocity and streamfunction for a Rankine vortex.
%   gaussianeddy - Velocity and streamfunction for a Gaussian vortex.
%   twolayereddy - Velocity and streamfunction for a 2-layer Rankine vortex.
%
% Water mass analysis
%   heatstorage  - Water column heat storage from 1-D mixing.
%
% Miscellaneous functions
%   radearth     - The radius of the earth in kilometers. 
%   tidefreq     - Frequencies of the eight major tidal compenents.
%   heat2evap    - Transform latent heat loss into units of evaporation.
%
% Reading various datasets
%   read_sand  - Read topography data from the Smith and Sandwell Database.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2011 J.M. Lilly --- type 'help jlab_license' for details        
 
help Joceans
if 0
          %Date conversion 
             monthfrac  %- Convert date from 'datenum' format to 'month.fraction'.
             yearfrac   %- Convert date from 'datenum' format to 'year.fraction'.
             yf2num     %- Converts date in 'year.fraction' format to  'datenum' format.
             mjd2num    %- Converts Modified Julian Dates to 'datenum' format.
          
           %Latitude and longitude
             corfreq    %- Coriolis frequency in cycles per hour.
             spheredist %- Computes great circle distances on a sphere.
             inregion   %- Tests whether lat/lon points lie within a specified box.
             latlon2xy  %- Converts latitude and longitude into local Cartesian coordinates.
             xy2latlon  %- Converts local Cartesian coordinates into latitude and longitude.
             latlon2uv  %- Converts latitude and logitude to velocity.
             deg360     %- Converts degrees to the range [0, 360].
             deg180     %- Converts degrees to the range [%-180, 180].
             lonshift   %- Shifts longitude origin for plotting purposes.
          
           %Plotting functions
             stickvect  %- Plots "stick vectors" for multicomponent velocity time series.
             denscont   %- Density contour overlay for oceanographic T/S plots.
             uvplot     %- Plot u and v components of velocity on the same axis.
             provec     %- Generate progressive vector diagrams (simple and fancy).
             hodograph  %- Generate hodograph plots (simple and fancy).
             latratio   %- Set plot aspect ratio for latitude / longitude plot.
             timelabel  %- Put month, day, or hour labels on a time axes.
          
           %Eddy modelling
             rankineeddy  %- Velocity and streamfunction for a Rankine vortex.
             gaussianeddy %- Velocity and streamfunction for a Gaussian vortex.
             twolayereddy %- Velocity and streamfunction for a 2%-layer Rankine vortex.
          
           %Water mass analysis
             heatstorage  %- Water column heat storage from 1%-D mixing.
          
           %Miscellaneous functions
             radearth     %- The radius of the earth in kilometers. 
             tidefreq     %- Frequencies of the eight major tidal compenents.
             heat2evap    %- Transform latent heat loss into units of evaporation.
          
          %Reading various datasets
             read_sand  %- Read topography data from the Smith and Sandwell Database. 
end

