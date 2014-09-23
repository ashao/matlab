tsfile = '/home/ashao/uw-apl/data/offtrac/input/normalyear/ts.nc';
hfile = '/home/ashao/uw-apl/data/offtrac/input/normalyear/H-clim.nc';
cfcfile = '/home/ashao/uw-apl/models/offtrac/active_runs/cfcsf6.normalyear/cfcs.sf6.saturation.0898.nc';

start4d = [0 0 0 0];
count4d = [1 -1 -1 -1];

offtrac.temp = nc_varget(tsfile,'temp',start4d,count4d);
offtrac.salt = nc_varget(tsfile,'salt',start4d,count4d);
offtrac.h = nc_varget(hfile,'h',start4d,count4d);

monidx = (2009-1936)*12;
start4d = [ (2009-1936)*12 0 0 0];
time = nc_varget(cfcfile,'Time',monidx,1);
cfc11.conc = nc_varget(cfcfile,'mn_cfc11',start4d,count4d);
cfc12.conc = nc_varget(cfcfile,'mn_cfc12',start4d,count4d);
sf6.conc = nc_varget(cfcfile,'mn_sf6',start4d,count4d);

trac = tracer_properties;
cfc11.sol = trac_calcsol(offtrac.temp+273.15,offtrac.salt, ...
    trac.cfc11.sol_coeffs.vol);
cfc12.sol = trac_calcsol(offtrac.temp+273.15,offtrac.salt, ...
    trac.cfc12.sol_coeffs.vol);
sf6.sol = trac_calcsol(offtrac.temp+273.15,offtrac.salt, ...
    trac.sf6.sol_coeffs.vol);

cfc11.atmc = cfc11.conc./cfc11.sol;
cfc12.atmc = cfc12.conc./cfc12.sol;
sf6.atmc = sf6.conc./sf6.sol;

%%
load metrics
clevels = 0:20:300;
hmask = offtrac.h > 1e-2;
for i=1:49
contourf(metrics.geolon.data,metrics.geolat.data,squeeze(cfc11.atmc(i,:,:).*hmask(i,:,:)))
colorbar
title(sprintf('Layer %d',i))
pause
end