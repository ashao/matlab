load /ltraid3/ashao/uw-apl/projects/saf_altimetry/ecmwf.taucurl.mat
%%
    ecmwf.taucurlanom = zeros(size(ecmwf.taucurl));
for mon = 1:12
   
    idx = mon:12:length(ecmwf.time);
    monmean = repmat(mean(ecmwf.taucurl(idx,:,:)),[length(idx) 1 1]);
    ecmwf.taucurlanom(idx,:,:) = ecmwf.taucurl(idx,:,:) - monmean;    
    
end

%%
[ntime nlat nlon] = size(ecmwf.taucurl);
ecmwf.taucurl_ann = zeros(ntime/12,nlat,nlon);
sidx = 1;
eidx = sidx+11;
counter = 0;
while eidx <= length(ecmwf.time);
   idx = sidx:eidx;
   counter = counter+1;
   ecmwf.taucurl_ann(counter,:,:) = mean(ecmwf.taucurl(idx,:,:));
   time_avg(counter) = mean(ecmwf.time(idx));
   sidx = sidx+12;
   eidx = sidx+11;
    
end


%%
for t = 1:ntime/12
    
   pcolor(ecmwf.longitude,ecmwf.latitude,squeeze(ecmwf.taucurl_ann(t,:,:)))
   shading flat;caxis([-1 1]*3e-7)
   drawnow
    
end

%%
latrange = [-60 -30];
lonrange = [0 360];

latidx = ecmwf.latitude >= min(latrange) & ecmwf.latitude <= max(latrange);
lonidx = ecmwf.longitude >= min(lonrange) & ecmwf.longitude <= max(lonrange);

acc.taucurl_ann = ecmwf.taucurl_ann(:,latidx,lonidx);
acc.time = time_avg;
acc.latitude = ecmwf.latitude(latidx);
acc.longitude = ecmwf.longitude(lonidx);
%%
notnan = ~isnan(acc.taucurl_ann(1,:,:));
data = detrend(acc.taucurl_ann(1:end,notnan));

[U S V] = svds(data);
ecmwf.time = acc.time;
ecmwf.pcs = U*S;
%%

% subplot(2,1,1)

nmode = 1;clf

subplot(2,1,1)
mode = zeros(length(acc.latitude),length(acc.longitude));
mode(notnan) = V(:,nmode);
worldmap([-90 max(latrange)],[0 360])
pcolorm(acc.latitude,acc.longitude,mode)
caxis([-1 1]*1e-2)
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5])
 subplot(2,1,2)
 hold on;
 plot(acc.time,pcs(:,nmode)./std(pcs(:,nmode)))
datetick
 xlim([datenum(1993,1,1) datenum(2014,1,1)])

 grid on