infile = '/ltraid3/ashao/uw-apl/models/offtrac/active_runs/oxygen/oxygen.12011.nc';
load metrics;
oxy = nc_varget(infile,'mn_oxygen');

%%
load metrics
for i = 1:49
   cla 
   data = squeeze(oxy(end,i,:,:));
   nneg = find(data(:)<0 & logical(metrics.wet.data(:)))
   metrics.geolon.data(nneg)
   metrics.geolat.data(nneg)
% axesm('MapProjection','Robinson','MapLatLimit',[-90 90],'MapLonLimit',[-280 80]);
% pcolorm(metrics.geolat.data,metrics.geolon.data,data)
% caxis([-0.1 0])
% colorbar
fprintf('Layer %d\n',i)
pause
end