load metrics
infile = 'C:\Users\ashao\data\phos_init.mat';
load(infile);

%%
clf
axesm('MapProjection','Robinson','MapLatLimit',[-90 90],'MapLonLimit',[-280 80]);
plotdata = squeeze(phos_out.data(1,1,:,:));
plotdata(~logical(metrics.wet.data))=NaN;
% plotdata = double(plotdata==0 & logical(metrics.wet.data));
cla
contourfm(metrics.geolat.data,metrics.geolon.data,plotdata);
% shading flat;
colorbar
caxis([0 2])