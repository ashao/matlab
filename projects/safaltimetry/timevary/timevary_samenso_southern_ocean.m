sstfile = 'C:\Users\ashao\Data\sst.mnmean.nc';
maskfile = 'C:\Users\ashao\Data\lsmask.nc';

sst.lat = double(nc_varget(sstfile,'lat'));
sst.lon = double(nc_varget(sstfile,'lon'));
sst.data = nc_varget(sstfile,'sst');
sst.time = nc_varget(sstfile,'time') + datenum(1800,1,1);
timeidx = sst.time > datenum(1992,1,1) & sst.time < datenum(2013,1,1);

mask = flipud(logical(nc_varget(maskfile,'mask')));
%%

enso = read_psd('http://www.esrl.noaa.gov/psd/gcos_wgsp/Timeseries/Data/nino34.long.data');
%%
ensoidx = interp1(enso.time,enso.idx,sst.time(timeidx));
ensoidx = (ensoidx-mean(ensoidx))/std(ensoidx);
solatidx = sst.lat < -30;
nsolat = sum(solatidx);

[ntime nlat nlon] = size(sst.data);
sst.cov.enso = zeros(nlat,nlon);
sst.so.data = sst.data(timeidx,solatidx,:);

%% Remove annual cycle


ntime = sum(timeidx);
sst.so.filtseasonal = zeros(ntime,nsolat,nlon);
for monidx = 1:12
    
   nyears = length(monidx:12:ntime);
   sst.so.filtseasonal(monidx:12:end,:,:) = squeeze(sst.so.data(monidx:12:end,:,:)) - ...
       repmat(mean(sst.so.data(monidx:12:end,:,:)),[nyears,1,1]);
    
    
end
%%
[ntime nlat nlon] = size(sst.so.data);
sst.so.maxcorr = zeros(nlat,nlon);
sst.so.maxlag = zeros(nlat,nlon);

for latidx = 1:sum(solatidx);
    for lonidx = 1:nlon
        
        dataanom = squeeze(sst.so.filtseasonal(:,latidx,lonidx));
        dataanom = (dataanom-mean(dataanom)/std(dataanom));
        [C lags] = xcorr(dataanom,ensoidx,160,'biased');                
        [sst.so.maxcorr(latidx,lonidx) maxidx] = max(abs(C));
        sst.so.maxlag(latidx,lonidx) = lags(maxidx);
            
    end
end

%%
sst.so.maxcorr(mask(latidx,:)) = NaN;
sst.so.maxlag(mask(latidx,:)) = NaN;
clf
cmap = 'BuDRd_12';

minlat = -90;
maxlat = -30;
mapprops = {'stereo','MapLonLimit',[-180 180], ...
    'MapLatLimit',[minlat maxlat],'Frame','On','Grid', 'off','MeridianLabel','Off', ...
    'ParallelLabel','Off','PLineLocation',[-60 -30],'MLineLocation',-180:60:120, ...
    'MLabelParallel',maxlat-2,'LabelRotation','on'};
axesm(mapprops{:})
contourfm(sst.lat(solatidx),sst.lon,double(sst.so.maxcorr),-0.5:0.025:0.5,'LineStyle','None')

% [cs v] = contourm(sst.lat(solatidx),sst.lon,double(sst.so.maxlag)/12,-10:10:10,'LineColor','black');

caxis([0 0.5])
colormap(othercolor(cmap))
cax= colorbar('SouthOutside');
xlabel(cax,'Correlation \rho_{OISST, NINO34}')

tightmap
geoshow('landareas.shp', 'FaceColor', [1 1 1]-0.5)