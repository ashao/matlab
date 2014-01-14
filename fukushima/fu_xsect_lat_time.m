m_proj('Mercator','lon',[-240 -80],'lat',[0 60]);
for i=1:240
    clf    
subplot(4,2,1)
latidx=lats(1);
contourf(repmat(array.geolon(latidx,:),[30 1]),-squeeze(array.depth(i,:,latidx,:)),double(squeeze(array.mn_inert(i,:,latidx,:))),20)
xlim([-240 -80])
ylim([-600 0])
caxis([0 2e-5])
title(sprintf('Latitude %3.1f Month: %3.3d',array.lath(latidx),i))
shading flat
% colorbar
subplot(4,2,2)
latidx=lats(2);
contourf(repmat(array.geolon(latidx,:),[30 1]),-squeeze(array.depth(i,:,latidx,:)),double(squeeze(array.mn_inert(i,:,latidx,:))),20);
xlim([-240 -80])
ylim([-600 0])
caxis([0 2e-5])
title(sprintf('Latitude %3.1f Month: %3.3d',array.lath(latidx),i))
shading flat
% colorbar
% colorbar
subplot(4,2,3)
latidx=lats(3);
contourf(repmat(array.geolon(latidx,:),[30 1]),-squeeze(array.depth(i,:,latidx,:)),double(squeeze(array.mn_inert(i,:,latidx,:))),20);
xlim([-240 -80])
ylim([-600 0])
caxis([0 2e-5])
title(sprintf('Latitude %3.1f Month: %3.3d',array.lath(latidx),i))
shading flat
% colorbar

subplot(4,2,4)
latidx=lats(4);
contourf(repmat(array.geolon(latidx,:),[30 1]),-squeeze(array.depth(i,:,latidx,:)),double(squeeze(array.mn_inert(i,:,latidx,:))),20);
xlim([-240 -80])
ylim([-600 0])
caxis([0 2e-5])
title(sprintf('Latitude %3.1f Month: %3.3d',array.lath(latidx),i))
shading flat
% colorbar
subplot(4,2,5)
latidx=lats(5);
contourf(repmat(array.geolon(latidx,:),[30 1]),-squeeze(array.depth(i,:,latidx,:)),double(squeeze(array.mn_inert(i,:,latidx,:))),20);
xlim([-240 -80])
ylim([-600 0])
caxis([0 2e-5])
title(sprintf('Latitude %3.1f Month: %3.3d',array.lath(latidx),i))
shading flat
% colorbar
subplot(4,2,6)
latidx=lats(6);
contourf(repmat(array.geolon(latidx,:),[30 1]),-squeeze(array.depth(i,:,latidx,:)),double(squeeze(array.mn_inert(i,:,latidx,:))),20);
xlim([-240 -80])
ylim([-600 0])
caxis([0 2e-5])
title(sprintf('Latitude %3.1f Month: %3.3d',array.lath(latidx),i))
shading flat
% colorbar


% colorbar
% title(sprintf('Latitude 45.5 Month: %d',i))
subplot(4,2,[7 8])
m_pcolor(array.geolon,array.geolat,double(squeeze(array.mn_cs137(i,1,:,:))))
colorbar
shading flat
m_grid
% pause(0.1)

drawnow
M(i)=getframe(gcf);
end
movie2avi(M, 'fukushima_lat_xsect.avi');
