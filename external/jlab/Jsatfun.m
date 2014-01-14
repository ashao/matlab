%
% JSATFUN  Satellite data treatment and design.
% 
%   See also: 
%   CONTENTS, JARRAY, JMATH, JPOLY, JGRAPH, JSTRINGS, JSTATS, JSIGNAL,
%   JELLIPSE, JCELL, VTOOLS, JOCEANS, JSPHERE, *JSATFUN*, JTRIADS, JPAPERS
%   _________________________________________________________________
%
% Alongtrack satellite data
%   turningpoint  - Find turning points, i.e. local extrema, in time series.
%   orbitbreaks   - Separate orbit into passes based on turning points.
%   ssh2eke       - Converts alongtrack sea surface height to eddy kinetic energy.
%
% Basic satellite geometry                        
%   ze2dist    - Converts beam zenith angle into distance to surface.          
%   ze2inc     - Converts beam zenith angle into incidence angle.              
%   latlon2zeaz  - Compute zenith and azimuth angles for satellite beam.
%   zeaz2latlon  - Compute latitude and longitude viewed by satellite beam.
%
% Aquarius satellite functions
%   aquaplot   - Plot Aquarius satellite radiometer footprint.
%   aquaprint  - Compute Aquarius satellite radiometer footprints.
%   aquasal    - Aquarius salinity change with brightness temperature.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2011 J.M. Lilly --- type 'help jlab_license' for details   

help Jsatfun

if 0
	 % Loading satellite altimetry data
	   pf_extract %- Extract alongtrack Pathfinder data from specified region.   
             pf_params  %- Load satellite parameters from Pathfinder format file.
          
           %Alongtrack satellite data
             trackfill     %- Despiking and filling for alongtrack satellite data.   
             track2grid    %- Interpolate alongtrack satellite data onto a grid.            
             turningpoint  %- Find turning points, i.e. local extrema, in time series.
             orbitbreaks   %- Separate orbit into passes based on turning points.
             ssh2eke       %- Converts alongtrack sea surface height to eddy kinetic energy.
          
           %Basic satellite geometry                        
             ze2dist    %- Converts beam zenith angle into distance to surface.          
             ze2inc     %- Converts beam zenith angle into incidence angle.              
             latlon2zeaz  %- Compute zenith and azimuth angles for satellite beam.
             zeaz2latlon  %- Compute latitude and longitude viewed by satellite beam.
          
           %Aquarius satellite functions
             aquaplot   %- Plot Aquarius satellite radiometer footprint.
             aquaprint  %- Compute Aquarius satellite radiometer footprints.
             aquasal    %- Aquarius salinity change with brightness temperature.
end

