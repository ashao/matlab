tpfile = 'C:\Users\ashao\Documents\GitHub\matlab\projects\safaltimetry\TraceTheorique_TP.nc';
fields = {'DataIndexes','Latitudes','Longitudes','NbPoints','Tracks'};
nfields = length(fields);
for ifield = 1:nfields
    fprintf('Extracting %s\n',fields{ifield})
    tpdata.(fields{ifield}) = nc_varget(tpfile,fields{ifield});
end

ntracks = length(tpdata.Tracks);


sidx = 1;
for tidx=1:ntracks
    eidx = sidx + tpdata.NbPoints(tidx) - 1;
    groundpath(tidx).number = tidx;
    groundpath(tidx).lat = tpdata.Latitudes(sidx:eidx);
    groundpath(tidx).lon = tpdata.Longitudes(sidx:eidx);
    sidx = eidx+1;
    fprintf('Working on Track %03d: %d Points\n',tidx,tpdata.NbPoints(tidx));
end

%% Verify groundpaths
clf
axesm('eckert4','MapLatLimit',[-90 90],'MapLonLimit',[0 360],'Frame','On','Grid','On')
load coast
plotm(lat, long,'k')
for tidx = 1:ntracks
    plotm(groundpath(tidx).lat,groundpath(tidx).lon,'LineWidth',2)
    pause(0.1)
end
%% Text extraction on an altimetry file
altfile = 'C:\Users\ashao\Data\dt_ref_global_j1_sla_vxxc_20041020_20041027_20100503.nc';
%altfile = 'C:\Users\ashao\Data\dt_ref_global_tp_sla_vxxc_19940406_19940413_20100503.nc';
fields = {'SLA','longitude','track','time','latitude','cycle'};
nfields = length(fields);
for ifield = 1:nfields
    
    fprintf('Extracting %s\n',fields{ifield})
    altdata.(fields{ifield}) = nc_varget(altfile,fields{ifield});
    
end

%%
ntracks = length(unique(altdata.track));

fill2d = nan(254,3129);
track.lat = fill2d;
track.lon = fill2d;
track.SLA = fill2d;


for itrack = 1:1
    
    tracknum = 1;
    idx = altdata.track == tracknum;
    lats = altdata.latitude(idx);
    lons = altdata.latitude(idx);
    sla = altdata.SLA(idx);
    
    [C ia ib] = intersect(lats,groundpath(tracknum).lat);
    
    track.lat(tracknum,ib) = lats(ia);
    track.lon(tracknum,ib) = lons(ia);
    track.SLA(tracknum,ib) = sla(ia);
    
    
    
    
    
    
end