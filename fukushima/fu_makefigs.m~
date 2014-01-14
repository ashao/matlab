
datapath='/scratch/data/offtrac/input/';
ncfile='/scratch/models/offtrac/runs/fukushima_inert/FUKUSHIMA_INERT.0722.nc';
geolat=nc_varget([datapath 'metrics.nc'],'geolat');
geolon=nc_varget([datapath 'metrics.nc'],'geolon');

fu_max_inventory;
m_proj('Mercator','lon',[-240 -110],'lat',[0 60]);
figure
for i=1:12
subplot(4,3,i)
m_contourf(geolon,geolat,squeeze(cs137_inv_time(12*i,:,:)))
colorbar
m_grid
m_coast;
title(sprintf('%d Months since April 2011',12*i))
end




sew_latidx=165;
sew_lonidx=131;
sew_cs137=nc_varget(ncfile,'mn_cs137',[0 0 sew_latidx sew_lonidx],[-1 -1 1 1]);
sew_depth=cumsum(nc_varget(ncfile,'mn_h',[0 0 sew_latidx sew_lonidx],[-1 -1 1 1]),2);
figure
hold on
plot(sew_cs137(60,:),-sew_depth(60,:))
plot(sew_cs137(120,:),-sew_depth(60,:),'-o')
plot(sew_cs137(180,:),-sew_depth(60,:),'-x')
legend('2016','2021','2026')
ylabel('Depth (m)')
xlabel('Dilution Factor')
title('Seward, AK 57.5\circN 149.5\circW') 


hi_latidx=126;
hi_lonidx=127;
figure
hold on

hi_cs137=nc_varget(ncfile,'mn_cs137',[0 0 hi_latidx hi_lonidx],[-1 -1 1 1]);
hi_depth=cumsum(nc_varget(ncfile,'mn_h',[0 0 hi_latidx hi_lonidx],[-1 -1 1 1]),2);
plot(hi_cs137(60,:),-hi_depth(60,:))
plot(hi_cs137(120,:),-hi_depth(60,:),'-o')
plot(hi_cs137(180,:),-hi_depth(60,:),'-x')
legend('2016','2021','2026')
ylabel('Depth (m)')
xlabel('Dilution Factor')
title('Honolulu, HI 18.5\circN 153.5\circW');

figure
subplot(3,1,1)
timeidx=60;
cs137=squeeze(nc_varget(ncfile,'mn_cs137',[timeidx 0 0 0],[1 -1 -1 -1]));
depth=-cumsum(squeeze(nc_varget(ncfile,'mn_h',[120 0 0 0],[1 -1 -1 -1])));
contourf(repmat(geolat(:,110)',[49 1]),squeeze(depth(:,:,110)),squeeze(cs137(:,:,110)))
axis([0 60 -500 0]) 
title('2016 along 170\circW')
colorbar
xlabel('Latitude (\circN)')
ylabel('Depth (m)')
timeidx=120;
subplot(3,1,2)
cs137=squeeze(nc_varget(ncfile,'mn_cs137',[timeidx 0 0 0],[1 -1 -1 -1]));
depth=-cumsum(squeeze(nc_varget(ncfile,'mn_h',[120 0 0 0],[1 -1 -1 -1])));
contourf(repmat(geolat(:,110)',[49 1]),squeeze(depth(:,:,110)),squeeze(cs137(:,:,110)))
axis([0 60 -500 0]) 
title('2021 along 170\circW')
colorbar
xlabel('Latitude (\circN)')
ylabel('Depth (m)')
timeidx=180;
subplot(3,1,3)
cs137=squeeze(nc_varget(ncfile,'mn_cs137',[timeidx 0 0 0],[1 -1 -1 -1]));
depth=-cumsum(squeeze(nc_varget(ncfile,'mn_h',[120 0 0 0],[1 -1 -1 -1])));
contourf(repmat(geolat(:,110)',[49 1]),squeeze(depth(:,:,110)),squeeze(cs137(:,:,110)))
axis([0 60 -500 0])
title('2026 along 170\circW')
colorbar
xlabel('Latitude (\circN)')
ylabel('Depth (m)')


figure
subplot(3,1,1)
timeidx=60;
cs137=squeeze(nc_varget(ncfile,'mn_cs137',[timeidx 0 0 0],[1 -1 -1 -1]));
depth=-cumsum(squeeze(nc_varget(ncfile,'mn_h',[120 0 0 0],[1 -1 -1 -1])));
contourf(repmat(geolon(141,:),[49 1]),squeeze(depth(:,141,:)),squeeze(cs137(:,141,:)))
title('2016 along 33.5\circN')
axis([-220 -120 -500 0])
colorbar
xlabel('Longitude (\circW)')
ylabel('Depth (m)')
timeidx=120;
subplot(3,1,2)
cs137=squeeze(nc_varget(ncfile,'mn_cs137',[timeidx 0 0 0],[1 -1 -1 -1]));
depth=-cumsum(squeeze(nc_varget(ncfile,'mn_h',[120 0 0 0],[1 -1 -1 -1])));
contourf(repmat(geolon(141,:),[49 1]),squeeze(depth(:,141,:)),squeeze(cs137(:,141,:)))
title('2021 along 170\circW')
colorbar
xlabel('Longitude (\circW)')
ylabel('Depth (m)')
axis([-220 -120 -500 0])
timeidx=180;
subplot(3,1,3)
cs137=squeeze(nc_varget(ncfile,'mn_cs137',[timeidx 0 0 0],[1 -1 -1 -1]));
depth=-cumsum(squeeze(nc_varget(ncfile,'mn_h',[120 0 0 0],[1 -1 -1 -1])));
contourf(repmat(geolon(141,:),[49 1]),squeeze(depth(:,141,:)),squeeze(cs137(:,141,:)))
title('2026 along 170\circW')
colorbar
xlabel('Longitude (\circW)')
ylabel('Depth (m)')
axis([-220 -120 -500 0])
