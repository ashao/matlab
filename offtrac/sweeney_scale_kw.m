%%
corepath = '/ltraid4/ashao/COREv2/clim/';
windfiles = {'u_10.15JUNE2009.nc','v_10.15JUNE2009.nc'};
core.lat=nc_varget([corepath windfiles{1}],'LAT');
core.lon=nc_varget([corepath windfiles{1}],'LON');
core.u=nc_varget([corepath windfiles{1}],'U_10_MOD');
core.v=nc_varget([corepath windfiles{2}],'V_10_MOD');
[ntime nlat nlon] = size(core.u);
core.utot = reshape(sqrt(core.u.^2+core.v.^2),[ntime nlat*nlon]);
[longrid latgrid]=meshgrid(core.lon,core.lat);
%% Get wet mask for HIM/GOLD
load metrics
[shiftlon shiftwet]=shiftlonrange(metrics.geolon.data+280,...
    metrics.wet.data,360);
wetmask = griddata(shiftlon,metrics.geolat.data,shiftwet, ...
    longrid,latgrid);
wetmask(wetmask < 0.5)=false;
wetmask(wetmask >= 0.5)=true;
%%
wts = cosd(latgrid).*wetmask;
wts = wts/nansum(wts(:));
wts(isnan(wts))=0;
wts(wts<=0)=0;
wts = wts(:);
glob.utotavg=core.utot*wts;

meanu10=sqrt((mean(glob.utotavg)^2+var(glob.utotavg)))
sweeney_scale = 15.4;
alpha = 15.7/(mean(glob.utotavg)^2+var(glob.utotavg))
% alpha = 14.6/(mean(glob.utotavg)^2+var(glob.utotavg,wts))
alpha_m_s = alpha * 2.7777777e-6