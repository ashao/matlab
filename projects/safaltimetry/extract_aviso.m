function [ array ] = extract_aviso( datapath )

tpfiles = dir([datapath filesep 'tp_cf' filesep '*vxxc*.nc]');
j1files = dir([datapath filesep 'j1_cf' filesep '*vxxc*.nc]');
j2files = dir([datapath filesep 'j2_cf' filesep '*vxxc*.nc]');

%% Get track lats/lons

fillmat = zeros(3129,254)*NaN;
trackpaths.lat = fillmat; % Assumes 3129 points along track (254 tracks)
trackpaths.lon = fillmat;

for fileidx = 1:length(tpfiles)
    
    datafile=[datapath filesep 'tp_cf' filesep tpfiles(fileidx).name];
    tracks = 
    

end