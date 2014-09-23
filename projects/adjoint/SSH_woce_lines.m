sshfile = '/ltraid4/ashao/HIM/hyak_store/NORMALYEAR/month/ice_month.nc';
ssh = squeeze(mean(nc_varget(sshfile,'SSH')));
%%
load metrics
geolat = metrics.geolat.data;
geolon = metrics.geolon.data;
iswet = logical(metrics.wet.data);
geolat = [geolat geolat(:,1)];
geolon = [geolon geolon(:,1)+360];
ssh(~iswet) = NaN;

%%
clf
data = [ssh ssh(:,1)]-nanmean(ssh(:));
axesm('gstereo','MapLatLimit',[-80 80],'MapLonLimit',[-60 300])
contourfm(geolat,geolon,data,-2:0.1:2,'LineStyle','none')
colormap(flipud(lbmap(256,'redblue')))
tightmap;plabel;mlabel;
linem([-70:60],-152*ones(size(-70:60)),'LineWidth',2,'Color','black');
linem([-65:20],95*ones(size(-65:20)),'LineWidth',2,'Color','black');
linem(30*ones(size(-230:-120)),-230:-120,'LineWidth',2,'Color','black');
linem(-30*ones(size(35:110)),35:110,'LineWidth',2,'Color','black');

linem([-70:60],-152*(ones(size(-70:60))),'LineWidth',2,'Color','black');
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]+0.2);
cax = colorbar;
caxis([-1.25 1.25])
set(cax,'ytick',-1.25:0.25:1.25)
ylabel(cax,'Sea Surface Height [m]')