%% Annually averaged SAM
mbfile = '/ltraid4/ecmwf/era-interim/monthly/700mb.geopotential.1979.2013.nc';
fields = {'longitude','latitude','z','time'};
for fidx = 1:length(fields)
    
    ecmwf.(fields{fidx}) = double(nc_varget(mbfile,fields{fidx}));
    
end
ecmwf.time = ecmwf.time/24 + datenum(1900,1,1);
%% Truncate time and latitude in accordance with SAM definition
latrange = [-90 -20];
lonrange = [ 0 360 ];
timerange = [datenum(1979,1,1) datenum(2001,1,1)-1];
latidx = findrange(ecmwf.latitude,latrange);
lonidx = findrange(ecmwf.longitude,lonrange);
timeidx = findrange(ecmwf.time,timerange);

time = ecmwf.time(timeidx);
lat = ecmwf.latitude(latidx);
lon = ecmwf.longitude(lonidx);
[longrid latgrid] = meshgrid(lon,lat);

wts = cosd(latgrid);
wts = wts./sum(wts(:));
wtgrid = [];
wtgrid(1,:,:) = wts;
wtgrid = repmat(wtgrid,[length(time),1,1]);

zhght_wt = ecmwf.z(timeidx,latidx,lonidx).*wtgrid;
zhght_wt = zhght_wt(:,:);
zhght_wt = detrend(zhght_wt')';
[ntime npts ] = size(zhght_wt);
nyears = 2000-1979;
%%
for mon = 1:12
   
    nlen = length(mon:12:ntime);
    zhght_wt(mon:12:end,:)=zhght_wt(mon:12:end,:) - ...
        repmat(mean(zhght_wt(mon:12:end,:)),[nlen,1,1]);;
    
end
%%
[U S V] = svds(zhght_wt,1);
clf
subplot(2,1,1)
sam.eof = zeros(size(latgrid));
sam.eof(:) = V;
time = time;
idx = U*S;
worldmap([-90 -20],[0 360]);
pcolorm(latgrid,longrid,sam.eof)
plotm(coast.lat,coast.long,'LineWidth',2)
subplot(2,1,2)
plot(sam.time,sam.idx)
datetick; grid on;
%% Annual Averaged SAM index
% Project all data onto the SAM mode
wtgrid = [];
wtgrid(1,:,:) = wts;
wtgrid = repmat(wtgrid,[length(ecmwf.time),1,1]);
zhght_wt = ecmwf.z(:,latidx,lonidx).*wtgrid;
zhght_wt = zhght_wt(:,:);
zhght_wt = detrend(zhght_wt')';
for mon = 1:12
   
    nlen = length(mon:12:length(ecmwf.time));
    zhght_wt(mon:12:end,:)=zhght_wt(mon:12:end,:) - ...
        repmat(mean(zhght_wt(mon:12:end,:)),[nlen,1,1]);;
    
end
sam.time = ecmwf.time;
sam.idx = zhght_wt*V;
sam.idx = sam.idx./std(sam.idx);
plot(sam.time,sam.idx,'LineWidth',2')
%%
sam.annual.time = 1979:2013 + 0.5;
for i = 1:length(sam.annual.time);
   
    timerange = [datenum(1979+i,1,1) datenum(1979+i+1,1,1)-1];
    timeidx = findrange(sam.time,timerange);
    sam.annual.idx(i) = nanmean(sam.idx(timeidx));
    
    
end
plot(sam.annual.time,sam.annual.idx)