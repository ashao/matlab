infiles.u = '/ltraid4/ncep/NCEPV1/surf/uwnd.mon.mean.nc';
infiles.v = '/ltraid4/ncep/NCEPV1/surf/vwnd.mon.mean.nc';
infiles.sam = '/ltraid3/ashao//uw-apl/data/ncep/sam_from_ncep.mat';
outpath = '~/uw-apl/figs/sam_analysis/';

%%
load(infiles.sam);
ncep.time = nc_varget(infiles.u,'time')/24+datenum(1,1,1);
ncep.lat = nc_varget(infiles.u,'lat');
ncep.lon = nc_varget(infiles.u,'lon');
ncep.u = nc_varget(infiles.u,'uwnd');
ncep.v = nc_varget(infiles.v,'vwnd');

%%
sam.pos.idx = sam_ncep.pc > 0.5;
sam.neg.idx = sam_ncep.pc < 0.5;
sam.pos.u = squeeze(mean(ncep.u(sam.pos.idx,:,:)));
sam.pos.v = squeeze(mean(ncep.v(sam.pos.idx,:,:)));
sam.neg.u = squeeze(mean(ncep.u(sam.neg.idx,:,:)));
sam.neg.v = squeeze(mean(ncep.v(sam.neg.idx,:,:)));
sam.pos.utot = sqrt(sam.pos.u.^2+sam.pos.v.^2);
sam.neg.utot = sqrt(sam.neg.u.^2+sam.neg.v.^2);

%%
colormap(othercolor('Purples9'))
m_proj('Stereographic','lon',0,'lat',-90,'radius',70);
subplot(1,2,1)
m_contourf(ncep.lon,ncep.lat,sam.pos.utot,-10:2:10)
m_grid;
m_coast('line','Color','Black','Linewidth',2);
caxis([0 10])
cax = colorbar;
title('SAM > -0.5')
ylabel(cax,'Wind Speed (m/s)');

subplot(1,2,2)
m_contourf(ncep.lon,ncep.lat,sam.neg.utot,-10:2:10)
% shading flat
m_grid;
m_coast('line','Color','Black','Linewidth',2);
caxis([0 10])
cax=colorbar;
title('SAM < -0.5')
ylabel(cax,'Wind Speed (m/s)');

saveas(gcf,[outpath 'utot.pos.neg.ncep.eps'],'epsc')
%%
clf
colormap(flipud(othercolor('RdYlBu10')))
m_proj('Stereographic','lon',0,'lat',-90,'radius',70);
m_contourf(ncep.lon,ncep.lat,sam.pos.utot-sam.neg.utot,-5:1:10)
m_grid;

caxis([-3 3])
m_coast('line','Color','Black','Linewidth',2);
cax=colorbar;
ylabel(cax,'Wind Speed (m/s)')
title('Positive SAM - Negative SAM')
saveas(gcf,[outpath 'utot.diff.ncep.eps'],'epsc')