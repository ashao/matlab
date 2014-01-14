esm2gpath = '/ltraid2/ashao/uw-apl/models/esm2g/historical/';
core2path = '/ltraid1/ashao/COREv2/clim/';

%% CORE 2 Forcing
core.lat = nc_varget([core2path 'ncar_rad.15JUNE2009.nc'],'LAT');
core.lon = nc_varget([core2path 'ncar_rad.15JUNE2009.nc'],'LON');
lwdn = nc_varget([core2path 'ncar_rad.15JUNE2009.nc'],'LWDN_MOD');
swdn = nc_varget([core2path 'ncar_rad.15JUNE2009.nc'],'SWDN_MOD');

[nday nlat nlon] = size(lwdn);
eidx = cumsum(eomday(2003,1:12));
core.lwdn = zeros(12,nlat,nlon);
core.swdn = zeros(12,nlat,nlon);
sidx = 1;

for mon = 1:12        
    avgidx = sidx:eidx(mon)
    core.lwdn(mon,:,:) = mean(lwdn(avgidx,:,:));
    core.swdn(mon,:,:) = mean(swdn(avgidx,:,:));
    sidx = eidx(mon)+1;
end

%% ESM2G

load metrics;

esm2g.lon = metrics.geolon.data;
esm2g.lat = metrics.geolat.data;
files.lwdn = dir([esm2gpath 'rlds*.nc']);
files.swdn = dir([esm2gpath 'rsntds*.nc']);
fileidx = 19:length(files.lwdn);
[nlat nlon] = size(esm2g.lon);
ntime = length(fileidx)*60;

lwdn=zeros(ntime,nlat,nlon);
swdn=zeros(ntime,nlat,nlon);

counter = 0;

sidx = 1;
for file = 1:length(fileidx)
    
    eidx = sidx+60-1;
    lwdn(sidx:eidx,:,:) = nc_varget([esm2gpath files.lwdn(file).name],'rlds');
    swdn(sidx:eidx,:,:) = nc_varget([esm2gpath files.swdn(file).name],'rsntds');
    sidx = eidx+1;
    
end

esm2g.swdn = zeros(12,nlat,nlon);
esm2g.lwdn = zeros(12,nlat,nlon);

for mon = 1:12
    
    esm2g.swdn(mon,:,:) = mean(swdn(mon:12:end,:,:));
    esm2g.lwdn(mon,:,:) = mean(lwdn(mon:12:end,:,:));    
    
end

%% Plot Patterns
m_proj('Stereographic','lon',0,'lat',-90)
mopts = 'm_coast(''patch'',[0 0 0]+0.5); m_grid;';
colormap(othercolor('BuDRd_12'))
clevels=-1000:5:1000;
crange = [-1000 1000];
for mon=1:12
   
   subplot(1,2,1) 
   m_contourf(core.lon,core.lat,core.lwdn(mon,:,:),clevels,'LineColor','None')
   
   eval(mopts);
   cax=colorbar('SouthOutside');
   xlabel(cax,'Shortwave Heat Flux')
   title(sprintf('Month %d CORE2',mon))
   caxis(crange);
   
   subplot(1,2,2)
   m_contourf(esm2g.lon,esm2g.lat,esm2g.lwdn(mon,:,:),clevels,'LineColor','None')
   
   caxis(crange);
   eval(mopts)
   cax=colorbar('SouthOutside');
   xlabel(cax,'Shortwave Heat Flux')
   title(sprintf('Month %d ESM2G',mon))
   caxis(crange)
   
   filename = sprintf('swdn_compare_mon_%02d.eps',mon);
   exportfig(gcf,filename,'color','cmyk')
   
end
