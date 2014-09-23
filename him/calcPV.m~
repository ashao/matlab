clear all

hfile = '/ltraid2/ashao/uw-apl/data/offtrac/input/normalyear/H-clim.nc';
h=squeeze(mean(nc_varget(hfile,'h')));
load metrics.mat
latgrid(1,:,:) = metrics.geolat.data;
latgrid = repmat(latgrid,[49 1 1]);
PV = sw_f(latgrid)./h;
dry = ~logical(metrics.wet.data);
PV(:,dry)=NaN;
%%
layerfile='/ltraid1/ashao/HIM/himw/him_sis/INPUT/Global_HIM_IC.nc';
rho=nc_varget(layerfile,'Layer');
interface=nc_varget(layerfile,'Interface');
deltarho = diff(interface);
rho = repmat(rho,[1 210 360]);
deltarho = repmat(deltarho,[1 210 360]);

PV = PV.*deltarho./rho;

%%
geolon= metrics.geolon.data;
geolat= metrics.geolat.data;
m_proj('Mercator','lon',[-260 -100],'lat',[10 60])
m_pcolor(geolon,geolat,PV(15,:,:))
colorbar
caxis([1e-10 1e-9])