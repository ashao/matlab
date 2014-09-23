load /ltraid3/ashao/uw-apl/projects/saf_altimetry/ecmwf.taucurl.mat
load ~/uw-apl/projects/saf_altimetry/longterm_saf.mat
load ~/uw-apl/projects/saf_altimetry/timevary.mat
%% Plot variance by track
clf
worldmap([-90 -30],[0 360])
timevary.mean.std = std(timevary.mean.lat,0,2);
timevary.mean.avglat = mean(timevary.mean.lat,2);
timevary.mean.avglon = mean(timevary.mean.lon,2);

%% Plot against bathymetry
subplot(1,2,1)
worldmap([-90 -30],[0 360])
[null sortidx] = sort(timevary.track.north.avglon);
scatterm(timevary.mean.avglat,timevary.mean.avglon,30,timevary.mean.std,'filled')
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5])
colorbar

subplot(1,2,2)
highidx = timevary.mean.std > mean(timevary.mean.std)*1.1;
worldmap([-90 -30],[0 360])
pcolorm(topo.lat,topo.lon,topo.z);
scatterm(timevary.mean.avglat(highidx),timevary.mean.avglon(highidx),'filled')
% contourm(ecmwf.latitude,ecmwf.longitude,squeeze(ecmwf.taucurl(1,:,:))*1e7,[-1 1],'LineColor','Black')
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5])
%% ACC Forcing
latrange = [-55 -40];
lonrange = [20 60];
timerange = [datenum(1992,1,1) datenum(2013,1,1)];

tidx = ecmwf.time >= min(timerange) & ecmwf.time < max(timerange);
latidx = ecmwf.latitude >= min(latrange) & ecmwf.latitude <= max(latrange);
lonidx = ecmwf.longitude >= min(lonrange) & ecmwf.longitude <= max(lonrange);

acc.taucurl = ecmwf.taucurl(tidx,latidx,lonidx);
acc.time = ecmwf.time(tidx);
acc.latitude = ecmwf.latitude(latidx);
acc.longitude = ecmwf.longitude(lonidx);
%% Calculate monthly anomalies
    acc.taucurlanom = zeros(size(acc.taucurl));
for mon = 1:12
   
    idx = mon:12:length(acc.time);
    monmean = repmat(mean(acc.taucurl(idx,:,:)),[length(idx) 1 1]);
    acc.taucurlanom(idx,:,:) = acc.taucurl(idx,:,:) - monmean;    
    
end

%% EOF
notnan = ~isnan(acc.taucurlanom(1,:));
[U S V] = svds(acc.taucurlanom(1:end,notnan));
pcs = U*S;
%%

nmode = 1;
clf;
subplot(2,1,1)
worldmap(latrange + [-10 10],lonrange +[ -10 10 ] )
mode = nan(length(acc.latitude),length(acc.longitude));
mode(notnan) = V(:,nmode);

pcolorm(acc.latitude,acc.longitude,mode)
colorbar
caxis([-1 1]*1e-1)
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5])
subplot(2,1,2)
area(acc.time,pcs(:,nmode))
grid on;
datetick

sam = sam_index;
tidx = sam.time >= min(timerange) & sam.time <= max(timerange);
corrcoef(sam.idx(tidx),pcs(:,nmode))
