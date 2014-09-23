start4d = [0 1 0 0]; count4d = [12 1 -1 -1];
h.clim = double(nc_varget( ...
    '/ltraid3/ashao/uw-apl/data/offtrac/input/H-clim.nc','h',start4d,count4d));
h.hind = double(nc_varget( ...
    '/ltraid3/ashao/uw-apl/data/offtrac/input/normalyear/H-hind.nc','h',start4d,count4d));
%%
dhdt.clim = double(diff(h.clim));
dhdt.hind = double(diff(h.hind));

%%
load metrics
colormap(othercolor('Blues9'))
subplot(2,1,1)
worldmap([-80 80],[-280 80])
pcolorm(metrics.geolat.data,metrics.geolon.data,squeeze(max(dhdt.clim)))
colorbar
caxis([0 500])

subplot(2,1,2)
worldmap([-80 80],[-280 80])
pcolorm(metrics.geolat.data,metrics.geolon.data,squeeze(max(dhdt.hind)))
colorbar
caxis([0 500])
%%
figure
worldmap([-80 80],[-280 80])
pcolorm(metrics.geolat.data,metrics.geolon.data,max(h.clim)-max(h.hind));
colorbar
