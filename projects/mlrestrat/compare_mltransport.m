ocfile = '/ltraid1/ashao/HIM/hyak_store/COMBINE/month/ocean_month.nc';
start=[0 0 0 0];
count=[120 2 -1 -1];
uh=nc_varget(ocfile,'uh',start,count);
vh=nc_varget(ocfile,'vh',start,count);

uhml=nc_varget(ocfile,'uhml',start,count);
vhml=nc_varget(ocfile,'vhml',start,count);

uhgm=nc_varget(ocfile,'uhGM',start,count);
vhgm=nc_varget(ocfile,'vhGM',start,count);

uhtot = uh+uhml+uhgm;
vhtot = vh+vhml+vhgm;

magtot = sqrt(uhtot.^2+vhtot.^2);
magml = sqrt(uhml.^2+vhml.^2);

%%
[ntime ndepth nlat nlon] = size(magml);
load metrics
geolon = metrics.geolon.data;
geolat = metrics.geolat.data;
m_proj('Mercator','lon',[-280 80],'lat',[-60 60])

colormap(othercolor('Spectral9'))
subplot(2,1,1)
m_pcolor(geolon,geolat,mean(atan2(vhml(:,1,:,:),uhml(:,1,:,:)))*180/pi)
shading flat;
caxis([-1 1]*180)
cax=colorbar;
ylabel(cax,'Transport Angle')
m_grid;
m_coast('patch',[0 0 0]);
title('UH Restrat. ML')

subplot(2,1,2)
m_pcolor(geolon,geolat,mean(atan2(vhtot(:,1,:,:),uhtot(:,1,:,:)))*180/pi)
shading flat;
caxis([-1 1]*180)
cax=colorbar;
ylabel(cax,'Transport Angle')
title('Large Scale Flow')
m_grid;
m_coast('patch',[0 0 0]);
print -r300 -dpng transport.angle.png
%%
colormap(othercolor('Set36'))
subplot(2,1,1)
m_pcolor(geolon,geolat,log10(mean(magml(:,1,:,:))))
shading flat;
caxis([4 6])
cax=colorbar;
ylabel(cax,'log10(Transport Magnitude)')
m_grid;
m_coast('patch',[0 0 0]);
title('UH Restrat. ML')

subplot(2,1,2)
m_pcolor(geolon,geolat,log10(mean(magtot(:,1,:,:))))
shading flat;
caxis([4 6])
cax=colorbar;
ylabel(cax,'log10(Transport Magnitude)')
m_grid;
m_coast('patch',[0 0 0]);
title('UH Restrat. ML')
print -r300 -dpng transport.magnitude.png
% packboth(2,2)