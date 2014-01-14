m_proj('Mercator','lon',[-240 -80],'lat',[0 60]);
for i=1:240
    clf
subplot(4,2,1)
lonidx=30;
contourf(repmat(array.geolat(:,lonidx)',[30 1]),-squeeze(array.depth(i,:,:,lonidx)),double(squeeze(array.mn_inert(i,:,:,lonidx))),20)
xlim([0 60])
ylim([-600 0])
caxis([0 5e-5])
title(sprintf('Longitude %3.1f Month: %3.3d',array.lonh(lonidx),i))
shading flat
% colorbar
subplot(4,2,2)
lonidx=45;
contourf(repmat(array.geolat(:,lonidx)',[30 1]),-squeeze(array.depth(i,:,:,lonidx)),double(squeeze(array.mn_inert(i,:,:,lonidx))),20)
xlim([0 60])
ylim([-600 0])
caxis([0 5e-5])
title(sprintf('Longitude %3.1f Month: %3.3d',array.lonh(lonidx),i))
shading flat
% colorbar
subplot(4,2,3)
lonidx=60;
contourf(repmat(array.geolat(:,lonidx)',[30 1]),-squeeze(array.depth(i,:,:,lonidx)),double(squeeze(array.mn_inert(i,:,:,lonidx))),20)
xlim([0 60])
ylim([-600 0])
caxis([0 5e-5])
title(sprintf('Longitude %3.1f Month: %3.3d',array.lonh(lonidx),i))
shading flat

subplot(4,2,4)
lonidx=75;
contourf(repmat(array.geolat(:,lonidx)',[30 1]),-squeeze(array.depth(i,:,:,lonidx)),double(squeeze(array.mn_inert(i,:,:,lonidx))),20)
xlim([0 60])
ylim([-600 0])
caxis([0 5e-5])
title(sprintf('Longitude %3.1f Month: %3.3d',array.lonh(lonidx),i))
shading flat
% colorbar
subplot(4,2,5)
lonidx=90;
contourf(repmat(array.geolat(:,lonidx)',[30 1]),-squeeze(array.depth(i,:,:,lonidx)),double(squeeze(array.mn_inert(i,:,:,lonidx))),20)
xlim([0 60])
ylim([-600 0])
caxis([0 5e-5])
title(sprintf('Longitude %3.1f Month: %3.3d',array.lonh(lonidx),i))
shading flat
% colorbar
subplot(4,2,6)
lonidx=105;
contourf(repmat(array.geolat(:,lonidx)',[30 1]),-squeeze(array.depth(i,:,:,lonidx)),double(squeeze(array.mn_inert(i,:,:,lonidx))),20)
xlim([0 60])
ylim([-600 0])
caxis([0 5e-5])
title(sprintf('Longitude %3.1f Month: %3.3d',array.lonh(lonidx),i))
shading flat


% colorbar
% title(sprintf('Latitude 45.5 Month: %d',i))
subplot(4,2,[7 8])
m_pcolor(array.geolon,array.geolat,double(squeeze(array.mn_cs137(i,1,:,:))))
colorbar
shading flat
caxis([0 5e-5])
m_grid
drawnow;
M(i)=getframe(gcf);
end
movie2avi(M, 'fukushima_lon_xsect.avi');