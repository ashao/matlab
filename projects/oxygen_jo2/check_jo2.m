inpath = '/home/ashao/uw-apl/models/offtrac/active_runs/oxygen/';
jo2 = nc_varget([inpath 'oxygen.woa09.2410.nc'],'mn_jo2',[2400-12 0 0 0],[12 -11 -1 -1]);
load metrics;
colormap(othercolor('BuDRd_12'))
%%
pcolor(metrics.geolon.data,metrics.geolat.data,squeeze(mean(jo2(:,1,:,:)))); 
shading flat
hold on
contour(metrics.geolon.data,metrics.geolat.data,squeeze(mean(jo2(:,1,:,:))),[0 0],'LineColor','black');
caxis([-1 1]*0.2e-8)
colorbar
%%
hfile= '/home/ashao/uw-apl/data/offtrac/input/normalyear/H-clim.nc';
h = nc_varget(hfile,'h');
depth = cumsum(h,2) - h./2;
depth(:,1,:,:) = 0;
%%
clf
plotidx = find(abs(metrics.geolat.data-30.5)<0.2 & abs(metrics.geolon.data+152.5)<0.5);
monidx = 1;
plot((jo2(monidx,:,plotidx)),(depth(monidx,:,plotidx)))
set(gca,'ydir','reverse')
ylim([0 1000])
xlim([-1 1]*2e-9)
%%
idx = abs(metrics.geolon.data + 152.5)<0.5;
lats = metrics.geolat.data(idx);
