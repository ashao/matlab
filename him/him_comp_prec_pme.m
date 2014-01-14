mainpath='/models/HIM/';

precfile=[mainpath 'ncar_precip_clim.nc'];
pmefile=[mainpath 'PmE.100yr.nc'];
runofffile=[mainpath 'runoff.nc'];
metricfile=[mainpath 'faces.360.210.nc'];
wet_grid=nc_varget(metricfile,'wet');
area_grid=nc_varget(metricfile,'Ah');

%% CORE Forcing
rho0=1035.0;
precarray=him_diag_ncar_prec(precfile);
precarray=structfun( @(x) (x./rho0), precarray, 'UniformOutput', false);
precarray.lat=precarray.lat.*rho0;
precarray.lon=precarray.lon.*rho0;

%% PmE HIM Output
pmearray=him_diag_pme(pmefile,metricfile,'PmE');
pmearray.runoff=nc_varget(runofffile,'runoff')./rho0;
pmearray.runoff(~wet_grid)=NaN;

%% Regrid runoff data
[ncarlon ncarlat]=meshgrid(precarray.lon,precarray.lat);
ncarlon=him_convert_to_grid_lon(ncarlon);
runoff_ncar=griddata(pmearray.geolon,pmearray.geolat,pmearray.runoff, ...
    ncarlon, ncarlat);


%% Account for runoff in annual precipitation and apply landmask
precarray.net_annual_prec=precarray.net_annual_prec+runoff_ncar;
precarray.mean_annual_prec=nanmean(makevec(precarray.net_annual_prec));

net_annual_pme_ncar=griddata(pmearray.geolon,pmearray.geolat, ...
    squeeze(nanmean(pmearray.net_annual_pme,1)), ncarlon, ncarlat);


