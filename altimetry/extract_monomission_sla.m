function track = extract_monomission_sla(datapath,outpath,tracknum)
%% extract_monomission_sla Concatenate along-track SLA for AVISO CF format
% Input:
%   datapath (char): Path where the datafiles are stored
%   tracknum (int #1-254): The number of the theoretical track
%  Output:
%   track (struct):
%       time [nfiles nSLA]: Time of the measurement in the Matlab serial day
%           NOTE: Reference to 0000-01-01 00:00:00
%                 NOT 1950-01-01 00:00 (AVISO's reference)
%       lat [nSLA]: latitude of measurement
%       lon [nSLA]; longitude of measurement
%       SLA [nfiles nSLA]: Sea level anomaly

% Split the data so they're organized according to mission order
tic;
tpfiles = dir([datapath filesep '*tp*.nc']);
j1files = dir([datapath filesep '*j1*.nc']);
j2files = dir([datapath filesep '*j2*.nc']);
allfiles = [tpfiles ; j1files ; j2files];
%%
ntp = length(tpfiles);
nj1 = length(j1files);
nj2 = length(j2files);
nfiles = length(allfiles);

% Following deltaT based on e-mail from AVISO
deltaT.tp = 1.07858;
deltaT.j1 = 1.02;
detlaT.j2 = 1.01941747572816;

lat=[];
lon=[];
time = [];

% Get the ground track based on the TP mission
fprintf('Getting ground track for track %d\n',tracknum)
for fidx = 1:ntp    
    fname = [datapath filesep tpfiles(fidx).name];
    templat = nc_varget(fname,'latitude');
    templon = nc_varget(fname,'longitude');
    temptrack = nc_varget(fname,'track');
    temptime = nc_varget(fname,'time');
    trackidx = temptrack == tracknum;
    
    [lat ilat]=unique([lat ; templat(trackidx)],'first');
    lon=unique([lon ; templon(trackidx)],'first');
    
end
%%

% Preallocate arrays and store the ground track
nSLA = length(lat);
fillmat2d = nan(nfiles,nSLA);

track.number=tracknum;
track.time = fillmat2d;
track.SLA = fillmat2d;
track.distflag = isnan(fillmat2d);
track.lat = lat;
track.lon = lon;

%%
fprintf('Extracting SLA (%d files): ',nfiles);
for fidx = 1:nfiles
    
    if mod(fidx,100)==0
        fprintf('%d ',fidx);
    end
    
    fname = [datapath filesep allfiles(fidx).name];
    templat = nc_varget(fname,'latitude');
    templon = nc_varget(fname,'longitude');
    temptrack = nc_varget(fname,'track');
    temptime = nc_varget(fname,'time');
    tempSLA = nc_varget(fname,'SLA');
    trackidx = temptrack == tracknum;
    templat=templat(trackidx);
    templon=templon(trackidx);
    temptime = temptime(trackidx);
    
    ndata = sum(trackidx);
 
    for datidx = 1:ndata;
    
     [dist idx] = min( (templat(datidx)-track.lat).^2 + (templon(datidx)-track.lon).^2);
     
     track.SLA(fidx,idx)=tempSLA(datidx);
     track.time(fidx,idx)=temptime(datidx);
     
     if dist > 0.1
         
        track.distflag(fidx,idx)=false;
     else
        track.distflag(fidx,idx) = true;
     end
    
    end

end
fprintf('\nTime Elapsed: %f\n',toc);
track.time = track.time + datenum(1950,1,1);

if sum(track.distflag(:))<length(track.distflag(:))
       fprintf('Distance is greater than 0.1 in Track: %d\n',tracknum);
end
        filename = sprintf('t%03d.mat',tracknum);
        save([outpath filename],'track');

