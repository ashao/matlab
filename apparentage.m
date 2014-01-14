layeridx=16;
midx=(1990-1936)*12;

tsfile='/ltraid2/ashao/uw-apl/data/offtrac/input/normalyear/ts.nc';
offtrac.temp=nc_varget(tsfile,'temp',[0 layeridx 0 0],[-1 1 -1 -1]);
offtrac.salt=nc_varget(tsfile,'salt',[0 layeridx 0 0],[-1 1 -1 -1]);

runpath='/ltraid2/ashao/uw-apl/models/offtrac/run_archive/cfcsf6.normalyear/';
satrun=[runpath 'cfcs.sf6.saturation.0899.nc'];
gasexrun=[runpath 'cfcs.sf6.control.0899.nc'];

mn_cfc11.conc.sat=1.0*nc_varget(satrun,'mn_cfc11',[midx layeridx 0 0],[12 1 -1 -1]);
mn_cfc11.conc.gasex=1.0*nc_varget(gasexrun,'mn_cfc11',[midx layeridx 0 0],[12 1 -1 -1]);

load metrics
wet=logical(metrics.wet.data);
mn_cfc11.sat(:,~wet)=NaN;
mn_cfc11.gasex(:,~wet)=NaN;


% Load files with solubility information
load /ltraid2/ashao/data/tracers/cfc11.mat

% Rereference temperature to surface for solubility calculations
P=ones(size(offtrac.salt))*2000;
PR=zeros(size(offtrac.salt));
offtrac.temp0=sw_ptmp(offtrac.salt,offtrac.temp,P,PR);

% Calculate solubility
sol=trac_calcsol(min(offtrac.temp0+274.15),nanmean(offtrac.salt),cfc11.sol);

% Calculate apparent age
mn_cfc11.atmcon.sat=nanmean(mn_cfc11.conc.sat)./sol;
mn_cfc11.atmcon.gasex=nanmean(mn_cfc11.conc.gasex)./sol;
yidx=cfc11.time<1992 & cfc11.time>1950;
mn_cfc11.age.sat=1990-interp1(cfc11.Nval(yidx),cfc11.time(yidx),mn_cfc11.atmcon.sat);
mn_cfc11.age.gasex=1990-interp1(cfc11.Nval(yidx),cfc11.time(yidx),mn_cfc11.atmcon.gasex);

% Plot difference
% subplot(3,1,1)
colormap(flipud(othercolor('BuDRd_12')))
lonrange=[-260 -80];
latrange=[20 65];
m_proj('Mercator','lon',lonrange,'lat',latrange)
m_contourf(metrics.geolon.data,metrics.geolat.data,...
    nanmean(mn_cfc11.conc.gasex-mn_cfc11.conc.sat),[-0.5:0.05:0],'LineColor','None');
% shading flat
title(sprintf('Average Difference in Concentation Layer %d',layeridx))
caxis([-0.5 0.5])
m_grid;
colorbar
% saveas(gcf,sprintf('~/figs/cfc11diff.layer%02d.eps',layeridx),'epsc');


subplot(2,1,1)

colormap((othercolor('BuDRd_12')))
m_proj('Mercator','lon',lonrange,'lat',latrange)
m_contourf(metrics.geolon.data,metrics.geolat.data,mn_cfc11.atmcon.sat,0:10:300);
m_coast('patch',[0.5 0.5 0.5]);
title(sprintf('Saturated Layer %d',layeridx))
caxis([0 270])
m_grid;
colorbar


subplot(2,1,2)

colormap((othercolor('BuDRd_12')))
m_proj('Mercator','lon',lonrange,'lat',latrange)
data=squeeze((mn_cfc11.age.gasex-mn_cfc11.age.sat)./mn_cfc11.age.sat*100);
m_contourf(metrics.geolon.data,metrics.geolat.data,data,0:10:80);
m_coast('patch',[0.5 0.5 0.5]);
title(sprintf('Saturated Layer %d',layeridx))
caxis([0 80])
m_grid;
colorbar
%%
latidx=(metrics.lath.data<55 & metrics.lath.data> 30);
lonidx=(metrics.lonh.data>-180 & metrics.lonh.data < -130);

npacmask=zeros(210,360);
npacmask(latidx,lonidx)=1;
npacmask=npacmask & wet;

meanerr=nansum(makevec(data.*npacmask.*metrics.Ah.data))./nansum(makevec(npacmask.*metrics.Ah.data));
meanagesat=nansum(makevec(squeeze(mn_cfc11.age.sat).*metrics.Ah.data.*npacmask))./ ...
    nansum(makevec(npacmask.*metrics.Ah.data));
meanagegasex=nansum(makevec(squeeze(mn_cfc11.age.gasex).*metrics.Ah.data.*npacmask))./ ...
    nansum(makevec(npacmask.*metrics.Ah.data));
fprintf('Layer %02d\n',layeridx+1);
fprintf('Mean Age Saturation: %f\n',meanagesat);
fprintf('Mean Age Gasex: %f\n',meanagegasex);
fprintf('Mean Percent Error %f\n\n',meanerr);

