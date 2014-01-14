function [ track ] = aviso_track_timeseries( datapath, restype, tracknumber )
% AVISO_TRACK_TIMESERIES: Extracts AVISO along-track altimetry by tracknumber
% into a format that allows for timeseries analysis
%
% [ track ] = aviso_track_timeseries(datapath, tracknumber)(
%
% Inputs:
%   datapath (cell): Contains all the directories containing altimetry data
%      Useful if altimetry files are separated by mission
%      E.g.: datapath={'/path/to/j1_cf/'; '/path/to/tp_cf'}
%   restype (string): Either 'vxxc' or 'vfec' if unfiltered or filtered
%       data, respectively, will be extracted
%   tracknumber (integer): Which track number to extract (1-254)
% Output:
%   track (structure)
%       tracknum: Number of the track
%       sla (nsla,nfiles): sea level anomaly (SLA)
%       lat (nsla): latitudes of SLA measurement
%       lon (nsla): longitude of SLA measruement
%       time (nsla,nfiles): timestamp of measurement
%
% DEPENDENCIES:
%   SNCTOOLS/MEXCDF: http://mexcdf.sourceforge.net/
% Author: Andrew E. Shao
% Email: ashao@uw.edu
% Date: 13 July 2013

npaths = length(datapath);

%% Loop over the first directory to reconstruct the satellite track
datafiles = dir([datapath{1} filesep '*' restype '*.nc']);
nfiles = length(datafiles);

lat = [];
lon = [];

for i=1:nfiles
    if mod(i,50)==0
        fprintf('Determining satellite track from file %d/%d\n',i,nfiles);
    end
    file = [datapath{1} filesep datafiles(i).name];
    alllat = nc_varget(file,'latitude');
    alllon = nc_varget(file,'longitude');
    tracks = nc_varget(file,'track');
    
    [lat unidx] = unique([lat ; alllat(tracks == tracknumber)]);
    lon = [lon ; alllon(tracks == tracknumber)];
    lon = lon(unidx);
    
    
end

track.lat = lat;
track.lon = lon;
nsla = length(track.lat);

%% Make the fill matrices to hold all the data
% Total number of files
ntotfiles = 0;
for i=1:npaths
    
    ntotfiles = ntotfiles + ...
        length(dir([datapath{i} filesep '*' restype '*.nc']));
    
end

fillmat2d = nan(nsla,ntotfiles);
track.sla = fillmat2d;
track.time = fillmat2d;

%% Extract the data
counter=0;
for i = 1:npaths
    files = dir([datapath{i} filesep '*' restype '*.nc']);
    nfiles = length(files);
    for j=1:nfiles
        
        counter = counter+1;
        if mod(counter,50)==0
            fprintf('Track %d: Extracting from file %d/%d\n',tracknumber, ...
                counter,ntotfiles);
        end
        % Extrac the data from netcdf file
        lat = nc_varget( [datapath{i} filesep files(j).name],'latitude');
        lon = nc_varget( [datapath{i} filesep files(j).name],'longitude');
        time = nc_varget( [datapath{i} filesep files(j).name],'time');
        sla = nc_varget( [datapath{i} filesep files(j).name],'SLA');
        tracks = nc_varget( [datapath{i} filesep files(j).name],'track');
        
        % Truncate the data for only the desired track
        trackidx = tracks == tracknumber;
        lat = lat(trackidx);
        lon = lon(trackidx);
        time = time(trackidx);
        sla = sla(trackidx);
        
        [null ia_lat ib_lat] = intersect(lat,track.lat);
        [null ia_lon ib_lon] = intersect(lon,track.lon);
        if length(ia_lat) ~= length(ia_lon)
            error('Length of matchinig lat/lon indices not equal!');
        end
        track.sla(ia_lat,counter)=sla(ia_lon);
        track.time(ia_lat,counter)=time(ia_lon);
        
        
    end
end
%% Add the tracknumber to outarray
track.tracknum = tracknumber;

