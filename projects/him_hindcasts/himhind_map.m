infile.mimoc = '/home/ashao/uw-apl/data/mimoc/MIMOC_ML_v2.2_PT_S_MLP_month08.nc';
mimoc.lat = nc_varget(infile.mimoc,'LATITUDE');
mimoc.long = nc_varget(infile.mimoc,'LONGITUDE');
mimoc.salt = nc_varget(infile.mimoc,'ABSOLUTE_SALINITY_MIXED_LAYER');
mimoc.temp = nc_varget(infile.mimoc,'CONSERVATIVE_TEMPERATURE_MIXED_LAYER');
mimoc.mldepth = nc_varget(infile.mimoc,'DEPTH_MIXED_LAYER');
mimoc.den = sw_pden(mimoc.salt,mimoc.temp,0,0);

mimoc = structfun(@double,mimoc,'UniformOutput',false);
%%
figure(1)
clf;
latidx = mimoc.lat < -20;
colormap(othercolor('BuOrR_14'))
hold on;

subplot(1,2,2)
% worldmap([-90 -20],[0 360]);
map_proj = { 'MapProjection','ortho', ...
    'MapLatLimit',[-90 -20],'Frame','on','Grid', 'On','MeridianLabel','On', ...
    'ParallelLabel','On','MLineLocation',-180:60:179, ...
    'LabelRotation','On','MLabelParallel',-35,'MlineException',0, ...
    'PLabelMeridian','prime', 'FontName','MyriadPro-Regular','Origin',[-90 60 0]} ;
axesm(map_proj{:})
contourfm(mimoc.lat(latidx),mimoc.long,mimoc.mldepth(latidx,:),0:10:400,'LineStyle','none')
[cs v]=contourm(mimoc.lat(latidx),mimoc.long,mimoc.den(latidx,:)-1000,[26.4 26.8 27.2],'LineColor','black');
clabelm(cs,v);
geoshow('landareas.shp','FaceColor',[0 0 0]+0.8)
cax = colorbar('SouthOutside');
xlabel(cax,'Mixed Layer Depth (m)')
%%
clf
% subplot(1,2,1); hold on
latidx = mimoc.lat >= -70 & mimoc.lat <= 20;
lonidx = mimoc.long >= 20 & mimoc.long <= 160;
worldmap([-70 -20],[20 160])
    setm(gca,'Frame','on','Grid', 'On','MeridianLabel','On', ...
    'ParallelLabel','On')
contourfm(mimoc.lat(latidx),mimoc.long(lonidx),mimoc.mldepth(latidx,lonidx),0:10:400,'LineStyle','none')
[cs v]=contourm(mimoc.lat(latidx),mimoc.long(lonidx),mimoc.den(latidx,lonidx)-1000,[26.4 26.8 27.2],'LineColor','black');
clabelm(cs,v);
i5lons = 30:115;
linem(ones(size(i5lons))*-32,i5lons,'LineWidth',3,'color','black')
linem(ones(size(i5lons))*-32,i5lons,'LineWidth',1,'color','red')
linem([-60 -20],[90 90],'LineWidth',3,'color','black')
linem([-60 -20],[90 90],'LineWidth',1,'color',[0.2 0.2 1])

geoshow('landareas.shp','FaceColor',[0 0 0]+0.8)
cax = colorbar('SouthOutside');
xlabel(cax,'Mixed Layer Depth (m)')
%% SAM
coast = load('coast.mat');
load /ltraid3/ashao/uw-apl/data/ncep/sam_from_ncep.mat
sam = structfun(@double,sam_ncep,'UniformOutput',false);
clf
% subplot(2,1,1)
map_proj = { 'MapProjection','ortho', ...
    'MapLatLimit',[-90 -20],'Frame','on','Grid', 'Off','MeridianLabel','On', ...
    'ParallelLabel','On','MLineLocation',-180:60:179, ...
    'LabelRotation','On','MLabelParallel',-35,'MlineException',0, ...
    'PLabelMeridian','prime', 'FontName','MyriadPro-Regular','Origin',[-90 60 0]} ;
axesm(map_proj{:})

contourfm(sam.lat,sam.lon,sam.eof*1000,-50:5:20,'LineStyle','None');
plotm(coast.lat,coast.long,'LineWidth',2,'Color','black')
cax=colorbar;
ylabel(cax,'700mb geopotential height anomaly (m)')
gridm on
box off
% caxis([-50 20])
% subplot(2,1,2)
% area(sam.time,smooth(sam.pc./std(sam.pc),12,'mean'))
% daspect([1e4 10 1])
% xlim([datenum(1948,1,1) datenum(2015,1,1)])
% datetick('keeplimits','keepticks')
% ylim([-1.5 1.5])
saveas(gcf,'/ltraid3/ashao/uw-apl/figs/him_hindcasts/sam_ncep.eps','epsc')
% grid on