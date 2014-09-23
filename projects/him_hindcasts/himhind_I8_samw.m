%% Load the data from himhind_PV.m
load ~/uw-apl/models/HIM/hindcast/hindcast.PV.mat
load metrics
%%
inpath.sigmatheta = '/ltraid3/ashao/uw-apl/models/HIM/hindcast/';
denfiles = dir([inpath.sigmatheta 'sigmatheta*.mat']);


%% Meridional Transects with salinity
figure
infile.salt = '/ltraid3/ashao/uw-apl/data/offtrac/input/normalyear/salt.hind.nc'
lineidx = 5;
wocenames = {'S1','A23','S2','I6','I8','I9','S3','P11','P14','P15', ...
    'P16','P17','P18','P19'};
wocelons = [-68.5 -30.5 0.5 30.5 -270.5 -245.5 -215.5 -205.5 -190.5 ...
    -170.5 -150.5 -135.5 -105.5 -90.5];
lon = wocelons(lineidx);
lonidx = find( abs(lon - metrics.lonh.data)<0.1 );
latgrid = repmat(metrics.geolat.data(:,lonidx)',[49 1]);
levels = [1e-20 linspace(1e-11,1e-9,30)];
colormap(flipud(othercolor('Blues5')))
outpath = '~/uw-apl/figs/him_hindcasts/';

layermask = ones(49,210);
% layermask(19:23,:)=0;
%%

for i=1:60
    clf
    depth = squeeze(hindcast.depth(i,:,:,lonidx));
    PV = abs(double(squeeze(hindcast.PV(i,:,:,lonidx))));
    start4d = [(i-1)*12 0 0 0];
    count4d = [12 -1 -1 -1];
    salt = squeeze(mean(nc_varget(infile.salt,'salt',start4d,count4d)));
    salt = squeeze(salt(:,:,lonidx));
    %     dens = squeeze(mean(pden(:,:,:,lonidx)));
    subplot(1,2,1)
    hold on
    contourf(latgrid,depth,log10(PV).*layermask,-11:0.1:-9,'LineColor','None')
%     pcolor(latgrid,depth,log10(PV).*layermask);caxis([-11 -9])
    load([inpath.sigmatheta denfiles(i).name]);
    
    dens = squeeze(mean(pden(:,:,:,lonidx)));
    [cs v] =contour(latgrid,depth,dens-1000,25:.2:28,'LineColor','Black');
    clabel(cs,v);
    set(gca,'ydir','reverse')
    xlim([-60 20])
    ylim([0 1000])
    caxis(log10([1e-2 1]*5e-10))
    caxis([-10.2 -9.5])
    
    cax=  colorbar('SouthOutside');
    xlabel('log10(PV)')
    title(sprintf('%d %s: %f',i+1947,wocenames{lineidx},lon));
    
    subplot(1,2,2)
    hold on
    contourf(latgrid,depth,salt,33:0.05:35,'LineColor','None')
%     pcolor(latgrid,depth,salt);caxis([33 35])
    
    load([inpath.sigmatheta denfiles(i).name]);
    
    dens = squeeze(mean(pden(:,:,:,lonidx)));
    [cs v] =contour(latgrid,depth,dens-1000,25:.2:28,'LineColor','Black');
    clabel(cs,v);
    set(gca,'ydir','reverse')
    xlim([-60 20])
    ylim([0 1000])
%     caxis(log10([1e-2 1]*5e-10))
%     caxis([-10.2 -9.5])
    
    cax=  colorbar('SouthOutside');
    xlabel('Salinity')
    title(sprintf('%d %s: %f',i+1947,wocenames{lineidx},lon));
    caxis([34 35])
    pause
end
%% PLot SAMW at three different points in time
lonidx = find(abs(-270.5-metrics.lonh.data)<0.1);

figure
yrange = [0 1500];
startyear = 1947;
latrange = [-60 -10];
pvrange = [-10.5 -9.5];
subplot(3,1,1) % 1948
i = 1955 - startyear;
denlevels = 26.4:0.2:27.4;
depth = squeeze(hindcast.depth(i,:,:,lonidx));
PV = abs(double(squeeze(hindcast.PV(i,:,:,lonidx))));
start4d = [(i-1)*12 0 0 0];
count4d = [12 -1 -1 -1];
% salt = squeeze(mean(nc_varget(infile.salt,'salt',start4d,count4d)));
% salt = squeeze(salt(:,:,lonidx));
%     dens = squeeze(mean(pden(:,:,:,lonidx)));
hold on
contourf(latgrid,depth,log10(PV).*layermask,-11:0.1:-9,'LineColor','None')
%     pcolor(latgrid,depth,log10(PV).*layermask);caxis([-11 -9])
load([inpath.sigmatheta denfiles(i).name]);
dens = squeeze(mean(pden(:,:,:,lonidx)));
[cs v] =contour(latgrid,depth,dens-1000,denlevels,'LineColor','Black');
clabel(cs,v);
set(gca,'ydir','reverse')
xlim(latrange)
ylim(yrange)
% caxis([34 35])
% caxis(log10([1e-2 1]*5e-10))
caxis(pvrange)
cax=  colorbar('EastOutside');
ylabel(cax,'log10(PV)')
title('I9 1955')

subplot(3,1,2) % 1948
i = 1975 - startyear;
depth = squeeze(hindcast.depth(i,:,:,lonidx));
PV = abs(double(squeeze(hindcast.PV(i,:,:,lonidx))));
start4d = [(i-1)*12 0 0 0];
count4d = [12 -1 -1 -1];
% salt = squeeze(mean(nc_varget(infile.salt,'salt',start4d,count4d)));
% salt = squeeze(salt(:,:,lonidx));
%     dens = squeeze(mean(pden(:,:,:,lonidx)));
hold on
contourf(latgrid,depth,log10(PV).*layermask,-11:0.1:-9,'LineColor','None')
%     pcolor(latgrid,depth,log10(PV).*layermask);caxis([-11 -9])
load([inpath.sigmatheta denfiles(i).name]);
dens = squeeze(mean(pden(:,:,:,lonidx)));
[cs v] =contour(latgrid,depth,dens-1000,denlevels,'LineColor','Black');
clabel(cs,v);
set(gca,'ydir','reverse')
xlim(latrange)
ylim(yrange)
% caxis([34 35])
caxis(pvrange)
cax=  colorbar('EastOutside');
ylabel(cax,'log10(PV)')
title('I9 1975')

subplot(3,1,3) % 1948
i = 1995 - startyear;
depth = squeeze(hindcast.depth(i,:,:,lonidx));
PV = abs(double(squeeze(hindcast.PV(i,:,:,lonidx))));
start4d = [(i-1)*12 0 0 0];
count4d = [12 -1 -1 -1];
% salt = squeeze(mean(nc_varget(infile.salt,'salt',start4d,count4d)));
% salt = squeeze(salt(:,:,lonidx));
%     dens = squeeze(mean(pden(:,:,:,lonidx)));
hold on
contourf(latgrid,depth,log10(PV).*layermask,-11:0.1:-9,'LineColor','None')
%     pcolor(latgrid,depth,log10(PV).*layermask);caxis([-11 -9])
load([inpath.sigmatheta denfiles(i).name]);
dens = squeeze(mean(pden(:,:,:,lonidx)));
[cs v] =contour(latgrid,depth,dens-1000,denlevels,'LineColor','Black');
clabel(cs,v);
set(gca,'ydir','reverse')
xlim(latrange)
ylim(yrange)
% caxis([34 35])
caxis(pvrange)
cax=  colorbar('EastOutside');
ylabel(cax,'log10(PV)')
title('I9 1995')
colormap(flipud(othercolor('Blues9')))
saveas(gcf,'/ltraid3/ashao/uw-apl/figs/him_hindcasts/I9_samw_1955_1975_1995.eps','epsc')
%% PLot SAMW at three different points in time
lonidx = find(abs(-270.5-metrics.lonh.data)<0.1);

figure
yrange = [0 1500];
startyear = 1947;
latrange = [-60 -10];
pvrange = [-10.5 -9.5];
subplot(3,1,1) % 1948
i = 1955 - startyear;
denlevels = 26.4:0.2:27.4;
depth = squeeze(hindcast.depth(i,:,:,lonidx));
PV = abs(double(squeeze(hindcast.PV(i,:,:,lonidx))));
start4d = [(i-1)*12 0 0 0];
count4d = [12 -1 -1 -1];
% salt = squeeze(mean(nc_varget(infile.salt,'salt',start4d,count4d)));
% salt = squeeze(salt(:,:,lonidx));
%     dens = squeeze(mean(pden(:,:,:,lonidx)));
hold on
contourf(latgrid,depth,log10(PV).*layermask,-11:0.1:-9,'LineColor','None')
%     pcolor(latgrid,depth,log10(PV).*layermask);caxis([-11 -9])
load([inpath.sigmatheta denfiles(i).name]);
dens = squeeze(mean(pden(:,:,:,lonidx)));
[cs v] =contour(latgrid,depth,dens-1000,denlevels,'LineColor','Black');
clabel(cs,v);
set(gca,'ydir','reverse')
xlim(latrange)
ylim(yrange)
% caxis([34 35])
% caxis(log10([1e-2 1]*5e-10))
caxis(pvrange)
cax=  colorbar('EastOutside');
ylabel(cax,'log10(PV)')
title('I9 1955')

subplot(3,1,2) % 1948
i = 1975 - startyear;
depth = squeeze(hindcast.depth(i,:,:,lonidx));
PV = abs(double(squeeze(hindcast.PV(i,:,:,lonidx))));
start4d = [(i-1)*12 0 0 0];
count4d = [12 -1 -1 -1];
% salt = squeeze(mean(nc_varget(infile.salt,'salt',start4d,count4d)));
% salt = squeeze(salt(:,:,lonidx));
%     dens = squeeze(mean(pden(:,:,:,lonidx)));
hold on
contourf(latgrid,depth,log10(PV).*layermask,-11:0.1:-9,'LineColor','None')
%     pcolor(latgrid,depth,log10(PV).*layermask);caxis([-11 -9])
load([inpath.sigmatheta denfiles(i).name]);
dens = squeeze(mean(pden(:,:,:,lonidx)));
[cs v] =contour(latgrid,depth,dens-1000,denlevels,'LineColor','Black');
clabel(cs,v);
set(gca,'ydir','reverse')
xlim(latrange)
ylim(yrange)
% caxis([34 35])
caxis(pvrange)
cax=  colorbar('EastOutside');
ylabel(cax,'log10(PV)')
title('I9 1975')

subplot(3,1,3) % 1948
i = 1995 - startyear;
depth = squeeze(hindcast.depth(i,:,:,lonidx));
PV = abs(double(squeeze(hindcast.PV(i,:,:,lonidx))));
start4d = [(i-1)*12 0 0 0];
count4d = [12 -1 -1 -1];
% salt = squeeze(mean(nc_varget(infile.salt,'salt',start4d,count4d)));
% salt = squeeze(salt(:,:,lonidx));
%     dens = squeeze(mean(pden(:,:,:,lonidx)));
hold on
contourf(latgrid,depth,log10(PV).*layermask,-11:0.1:-9,'LineColor','None')
%     pcolor(latgrid,depth,log10(PV).*layermask);caxis([-11 -9])
load([inpath.sigmatheta denfiles(i).name]);
dens = squeeze(mean(pden(:,:,:,lonidx)));
[cs v] =contour(latgrid,depth,dens-1000,denlevels,'LineColor','Black');
clabel(cs,v);
set(gca,'ydir','reverse')
xlim(latrange)
ylim(yrange)
% caxis([34 35])
caxis(pvrange)
cax=  colorbar('EastOutside');
ylabel(cax,'log10(PV)')
title('I9 1995')
colormap(flipud(othercolor('Blues9')))
saveas(gcf,'/ltraid3/ashao/uw-apl/figs/him_hindcasts/I9_samw_1955_1975_1995.eps','epsc')
%% PLot AAIW at three different points in time along I5
figure
latidx = find(abs(-32.5-metrics.lath.data)<0.1);
yrange = [0 1500];
startyear = 1947;
latrange = [-280 80];
salrange = [34.2 35.2];
denlevels = 26.4:0.2:27.4;
latgrid = repmat(metrics.geolon.data(latidx,:),[49 1]);

subplot(3,1,1) % 1948
i = 1955 - startyear;
depth = squeeze(hindcast.depth(i,:,latidx,:));
PV = abs(double(squeeze(hindcast.PV(i,:,latidx,:))));
start4d = [(i-1)*12 0 0 0];
count4d = [12 -1 -1 -1];
salt = squeeze(mean(nc_varget(infile.salt,'salt',start4d,count4d)));
salt = squeeze(salt(:,latidx,:));
%     dens = squeeze(mean(pden(:,:,:,lonidx)));
hold on
contourf(latgrid,depth,salt,33:0.05:35,'LineColor','None')
%     pcolor(latgrid,depth,log10(PV).*layermask);caxis([-11 -9])
load([inpath.sigmatheta denfiles(i).name]);
dens = squeeze(mean(pden(:,:,latidx,:)));
[cs v] =contour(latgrid,depth,dens-1000,denlevels,'LineColor','Black');
clabel(cs,v);
set(gca,'ydir','reverse')
xlim(latrange)
ylim(yrange)
caxis(salrange)
% caxis(log10([1e-2 1]*5e-10))
% caxis([-10.2 -9.5])
cax=  colorbar('EastOutside');
ylabel(cax,'Salinity')
title('I9 1955')

subplot(3,1,2) % 1948
i = 1975 - startyear;
depth = squeeze(hindcast.depth(i,:,latidx,:));
PV = abs(double(squeeze(hindcast.PV(i,:,latidx,:))));
start4d = [(i-1)*12 0 0 0];
count4d = [12 -1 -1 -1];
salt = squeeze(mean(nc_varget(infile.salt,'salt',start4d,count4d)));
salt = squeeze(salt(:,latidx,:));
%     dens = squeeze(mean(pden(:,:,:,lonidx)));
hold on
contourf(latgrid,depth,salt,33:0.05:35,'LineColor','None')
%     pcolor(latgrid,depth,log10(PV).*layermask);caxis([-11 -9])
load([inpath.sigmatheta denfiles(i).name]);
dens = squeeze(mean(pden(:,:,latidx,:)));
[cs v] =contour(latgrid,depth,dens-1000,denlevels,'LineColor','Black');
clabel(cs,v);
set(gca,'ydir','reverse')
xlim(latrange)
ylim(yrange)
caxis(salrange)
cax=  colorbar('EastOutside');
ylabel(cax,'Salinity')
title('I9 1975')

subplot(3,1,3) % 1948
i = 1995 - startyear;
depth = squeeze(hindcast.depth(i,:,latidx,:));
PV = abs(double(squeeze(hindcast.PV(i,:,latidx,:))));
start4d = [(i-1)*12 0 0 0];
count4d = [12 -1 -1 -1];
salt = squeeze(mean(nc_varget(infile.salt,'salt',start4d,count4d)));
salt = squeeze(salt(:,latidx,:));
%     dens = squeeze(mean(pden(:,:,:,lonidx)));
hold on
contourf(latgrid,depth,salt,33:0.05:35,'LineColor','None')
%     pcolor(latgrid,depth,log10(PV).*layermask);caxis([-11 -9])
load([inpath.sigmatheta denfiles(i).name]);
dens = squeeze(mean(pden(:,:,latidx,:)));
[cs v] =contour(latgrid,depth,dens-1000,denlevels,'LineColor','Black');
clabel(cs,v);
set(gca,'ydir','reverse')
xlim(latrange)
ylim(yrange)
caxis(salrange)
cax=  colorbar('EastOutside');
ylabel(cax,'Salinity')
title('I9 1995')
colormap(flipud(othercolor('Blues9')))
% saveas(gcf,'/ltraid3/ashao/uw-apl/figs/him_hindcasts/I9_aaiw_1955_1975_1995.eps','epsc')
%% Hovmoeller diagram of Indian SAMW

colormap(flipud(othercolor('Blues6')))
lineidx = 5;
wocenames = {'S1','A23','S2','I6','I8','I9','S3','P11','P14','P15', ...
    'P16','P17','P18','P19'};
wocelons = [-68.5 -30.5 0.5 30.5 -270.5 -245.5 -215.5 -205.5 -190.5 ...
    -170.5 -150.5 -135.5 -105.5 -90.5];
lon = wocelons(lineidx);
lat = -35.5;
lonidx = find( abs(lon - metrics.lonh.data)<0.1 );
latidx = find( abs(lat - metrics.lath.data)<0.1 );
layidx = 21;

hov.PV.meridional = squeeze(hindcast.PV(:,layidx,:,lonidx));
hov.PV.zonal = squeeze(hindcast.PV(:,layidx,latidx,:));
time = 1948:2007;
clims = [-10.3 -9.5];
options.contour = {-10.6:0.15:-9,'Linecolor','black','LineWidth',1.5};
options.contourf = {-10.6:0.05:-9,'LineColor','None'};
clf
subplot(2,1,1); hold on;
set(gcf,'Renderer','painters')
contourf(time,metrics.geolat.data(:,lonidx),log10(hov.PV.meridional'),options.contourf{:});
contour(time,metrics.geolat.data(:,lonidx),log10(hov.PV.meridional'),options.contour{:})
xlabel('Year'); ylabel('Latitude')
caxis(clims)
ylim([-40 -10])
cax = colorbar;
ylabel(cax, 'log10(PV)')
grid on; box on;
title(sprintf('Meridional: %3.1f',lon))
subplot(2,1,2); hold on;
contourf(time,metrics.geolon.data(latidx,:),log10(hov.PV.zonal'),options.contourf{:});
contour(time,metrics.geolon.data(latidx,:),log10(hov.PV.zonal'),options.contour{:});
xlabel('Year'); ylabel('Longitude')
ylim([-279.5 -245.5])
caxis([-10.4 -9.5])
cax = colorbar;
ylabel(cax, 'log10(PV)')
grid on; box on;
title(sprintf('Zonal: %3.1f',lat))

%% SAMW properties

avg.pv = zeros(60,1);
avg.dens = zeros(60,1);
avg.vol = zeros(60,1);

layers = 19:23;
geomask = [];
geomask(1,:,:) = metrics.geolon.data > -250 & metrics.geolat.data < -15 & metrics.geolat.data > -40;
geomask = repmat(geomask,[length(layers),1,1]);
Ah = [];
Ah(1,:,:) = metrics.Ah.data;
Ah = repmat(Ah,[length(layers),1,1]);
for i = 1:60
   fprintf('%d/%d...',i,60);
    load([inpath.sigmatheta denfiles(i).name]);
    PV = abs(squeeze(hindcast.PV(i,layers,:,:)));        
    avgmask = geomask & PV < 1.5e-10;
    
    h = squeeze(hindcast.h(i,layers,:,:));
    pden = squeeze(mean(pden(:,layers,:,:)));
    vol = h.*Ah.*avgmask;
    wts = vol./nansum(makevec(vol));
    
    avg.pv(i) = nansum(makevec(PV.*wts));
    avg.vol(i) = nansum(makevec(vol));
    avg.dens(i) = nansum(makevec(pden.*wts));
    fprintf('Done!\n')
    
    
end
%% Plot SAMW Properties and compare to SAM
years= 1948:2007;
load /home/ashao/uw-apl/data/ncep/sam_from_ncep.mat

for t = 1:length(years)
   
    sam.year(t) = years(t);
    tidx = sam_ncep.time >= datenum(years(t),1,1) & sam_ncep.time < datenum(years(t)+1,1,1);
    
    sam.idx(t) = mean(sam_ncep.pc(tidx));
    
end


subplot(4,1,1)
plot(years,avg.pv,'k-x')
xlim([1948 2007])
ylabel('Ertel PV')
%  title(sprintf('Mean: %2.1e STD: %2.1e',mean(avg.pv),std(avg.pv)))
[C lag]=maxcorr(avg.pv,sam.idx,10,'coeff');
title(sprintf('Correlation with SAM: %2.1f Lag: %d Years',C,lag))
grid on


subplot(4,1,2)
plot(years,avg.dens-1000,'k-x')
xlim([1948 2007])
ylabel('Density (kg m^{-3})')
%  title(sprintf('Mean: %3.1 STD: %2.1e',mean(avg.dens-1000),std(avg.dens-1000)))
[C lag]=maxcorr(avg.dens,sam.idx,10,'coeff');
title(sprintf('Correlation with SAM: %2.1f Lag: %d Years',C,lag))
grid on

subplot(4,1,3)
plot(years,avg.vol,'k-x')
xlim([1948 2007])
ylabel('Volume (m^3)')
[C lag]=maxcorr(avg.vol*1e-15,sam.idx,10,'coeff');
title(sprintf('Correlation with SAM: %2.1f Lag: %d Years',C,lag))

%  title(sprintf('Mean: %2.1e STD: %2.1e',mean(avg.vol),std(avg.vol)))
grid on

subplot(4,1,4)
plot(sam.year,sam.idx/std(sam.idx),'k-x')
xlim([1948 2007])
ylabel('SAM Index')
ylim([-2.5 2.5])
% xlabel(sprintf('Mean: %2.1e STD: %2.1e',mean(avg.vol),std(avg.vol)))
grid on

saveas(gcf,'/ltraid3/ashao/uw-apl/figs/him_hindcasts/indian_samw_props.eps')
