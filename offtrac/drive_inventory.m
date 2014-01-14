
%% Extract Data
% ncfile='/models/offtrac/gasex/gasex_nonlinear_atmospheric.0863.nc';
 ncfile='/scratch/models/offtrac/kw_parameterizations/wanninkhof.0863.nc';
%  ncfile='/models/offtrac/gasex/clim_test.0863.nc';
metricfile='/scratch/data/offtrac/input/metrics.nc';

offinv=off_comp_inventory(ncfile,metricfile);
% Correct model time
offinv.time=offinv.time+datenum(1936,1,1);

%% Plot the inventories of each tracer

h1=figure(1)
clf
subplot(2,1,1)
hold on
plot(offinv.time,cumsum(offinv.cfc11_mon_global_flux)./offinv.cfc11_global_inv,'k-','LineWidth',2);
plot(offinv.time,cumsum(offinv.cfc12_mon_global_flux)./offinv.cfc12_global_inv,'b-','LineWidth',2);
grid on
datetick
title('Ratio of (Cumulative Flux)/(Global Inventory)')
xlabel('Model Year')
legend('CFC-11','CFC-12','Location','Best')
subplot(2,1,2)
title('Global Surface flux')
hold on
plot(offinv.time,offinv.cfc11_mon_global_flux,'k-','LineWidth',2)
plot(offinv.time,offinv.cfc12_mon_global_flux,'b-','LineWidth',2)
grid on
datetick
saveas(h1,'/scratch/figs/sabine_meeting/cfcs_inventory_flux.eps','epsc')
