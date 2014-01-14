%% Load the data from himhind_PV.m
load ~/uw-apl/models/HIM/hindcast/hindcast.PV.mat
load metrics
%%
inpath.sigmatheta = '/ltraid3/ashao/uw-apl/models/HIM/hindcast/';
denfiles = dir([inpath.sigmatheta 'sigmatheta*.mat']);

%% North Pacific Mode Water Surface plots
hindcast.PV(hindcast.h<1e-6) = NaN;
%
colormap(flipud(othercolor('RdYlBu11')))
m_proj('Robinson','lon',[-250 -100],'lat',[0 50])

for i=1:60
    layidx = 11:17;
    clf
    subplot(2,1,1)
    hold on
    m_contourf(metrics.geolon.data,metrics.geolat.data, ...
        min(abs(hindcast.PV(i,layidx,:,:)),[],2),linspace(1e-11,1e-9,40),'LineColor','None')
    m_contour(metrics.geolon.data,metrics.geolat.data, ...
        min(abs(hindcast.PV(i,layidx,:,:)),[],2),[2.5e-10 2.5e-10],'LineColor','Black')
    % shading flat
    m_coast('patch',[0 0 0]);
    caxis([1e-10 1e-9])
    m_grid;
    cax = colorbar;
    ylabel(cax,'PV')
    title(sprintf('Year %d',1948+i-1))
    subplot(2,1,2)ot
    
    % m_proj('Mercator','lon',[-240 -120],'lat',[10 50])
    m_contourf(metrics.geolon.data,metrics.geolat.data, ...
        sum(hindcast.h(i,layidx,:,:),2),0:10:500,'LineColor','none');
    % shading flat
    m_coast('patch',[0 0 0]);
    caxis([0 400])
    m_grid;
    cax = colorbar;
    ylabel(cax,'Layer Thickness')
    drawnow
end

%% North Pacific Mode Water Properties
layidx = 11:17;
hindcast.npac.PV_vol=zeros(60,1);
layer = [32.978 33.288 33.575 33.846 34.103 34.347 34.577];
load metrics
geoidx = metrics.geolon.data < -135 & ...
    metrics.geolon.data > -240 & ...
    metrics.geolat.data > 15 & ...
    metrics.geolat.data < 45;
clear temp
temp(1,:,:)=metrics.Ah.data;
Ahgrid = repmat(temp,[length(layidx),1,1]);
clear temp
temp(1:7,1,1) = layer;
pden = repmat(temp,[1 210 360]);

for i=1:60
    
    pvdata= squeeze(hindcast.PV(i,layidx,geoidx));
    voldata = squeeze(hindcast.h(i,layidx,geoidx)).* ...
        Ahgrid(:,geoidx);
    pdendata = pden(:,geoidx);
    pvidx = pvdata < 2.5e-10;
    voldata(~pvidx) = 0;
    pdendata(~pvidx) = 0;
    hindcast.npac.cmw.vol(i) = sum(voldata(:));
    hindcast.npac.cmw.meanPV(i) = nansum(nansum(voldata.*pvdata))./nansum(voldata(:));
    hindcast.npac.cmw.meanden(i) = nansum(nansum(voldata.*pdendata))./nansum(voldata(:));
end

subplot(3,1,1)
plot(1948:2007,hindcast.npac.cmw.vol-mean(hindcast.npac.cmw.vol),...
    'k-x','LineWidth',2);
xlim([1948 2007])
grid on
ylabel('Volume Anomaly')
title('North Pacific Central Mode Water')
subplot(3,1,2)
plot(1948:2007,hindcast.npac.cmw.meanPV-mean(hindcast.npac.cmw.meanPV), ...
    'k-x','LineWidth',2)
xlim([1948 2007])
ylabel('PV Anomaly')
grid on;
subplot(3,1,3)
plot(1948:2007,hindcast.npac.cmw.meanden-mean(hindcast.npac.cmw.meanden),'k-x','LineWidth',2)
xlim([1948 2007])
ylabel('Density Anomaly')
grid on;

%% North Pacific Transect

latidx = find(metrics.lath.data==20.5);
longrid = repmat(metrics.geolon.data(latidx,:),[49 1]);
colormap(flipud(othercolor('RdYlBu9')))
years = 1948:2007
outpath = '/ltraid3/ashao/uw-apl/figs/him_hindcasts/PV/';
% Set some plotting commands common to both subplots
toplim = 800;
labels.x = 'Longitude';
labels.y = 'Depth (m)';
labels.cax = 'PV';
plotcmds = ['caxis([5e-11 5e-10]);' ...
    'ylabel(labels.y);xlabel(labels.x);' ...
    'set(gca,''ydir'',''reverse'');' ...
    ll  'xlim([-240 -110]);'];

for i=1:1
    
    clf; hold on;
    load([inpath.sigmatheta denfiles(i).name]);
    % NOTE: pcolor shows the actual grid point by grid point temperature
    contourf(longrid, squeeze(hindcast.depth(i,:,latidx,:)), ...
        double(squeeze(hindcast.PV(i,:,latidx,:))),linspace(5e-11,5e-10,30),'LineColor','None')
    % shading flat
    [cs v] =contour(longrid,squeeze(hindcast.depth(i,:,latidx,:)), ...
        squeeze(mean(pden(:,:,latidx,:)))-1000,24.7:.2:27,'LineColor','Black');
    clabel(cs,v);
    
    
    ylim([0 toplim])
    eval(plotcmds)
    
    cax = colorbar('SouthOutside');
    xlabel(cax,labels.cax)
    title(sprintf('Year %d',years(i)))
    
    % saveas(gcf,[outpath sprintf('PV.%d.eps',years(i))],'epsc');
    
end

%% Find Pacific SAMW
hindcast.PV(hindcast.h<1e-6) = NaN;
%
colormap(flipud(othercolor('RdYlBu11')))
m_proj('Robinson','lon',[-220 -70],'lat',[-60 0])

% for layidx = 5:49
%     disp(layidx)
for i=1:60
    layidx = 14:21;
    clf
    subplot(2,1,1)
    hold on
    m_contourf(metrics.geolon.data,metrics.geolat.data, ...
        min(abs(hindcast.PV(i,layidx,:,:)),[],2),linspace(1e-11,1e-9,40),'LineColor','None')
    m_contour(metrics.geolon.data,metrics.geolat.data, ...
        min(abs(hindcast.PV(i,layidx,:,:)),[],2),[1.5e-10 1.5e-10],'LineColor','Black')
    % shading flat
    m_coast('patch',[0 0 0]);
    caxis([1e-10 1e-9])
    m_grid;
    cax = colorbar;
    ylabel(cax,'PV')
    title(sprintf('Year %d',1948+i-1))
    subplot(2,1,2)
    
    % m_proj('Mercator','lon',[-240 -120],'lat',[10 50])
    m_contourf(metrics.geolon.data,metrics.geolat.data, ...
        sum(hindcast.h(i,layidx,:,:),2),0:10:500,'LineColor','none');
    % shading flat
    m_coast('patch',[0 0 0]);
    caxis([0 500])
    m_grid;
    cax = colorbar;
    ylabel(cax,'Layer Thickness')
    drawnow
end

%% Zonal Plot transects of PV
latidx = find( abs(-30.5 - metrics.lath.data)<0.1 );
latavg = (latidx-10):(latidx+10);
longrid = repmat(metrics.geolon.data(latidx,:),[49 1]);
longrid = [longrid longrid+360];
levels = [1e-20 linspace(1e-11,1e-9,30)];

sectorlons = {[30 120],[-210 -60],[-60 30]};
sectortitle = {'Indian','Pacific','Atlantic'}
colormap(othercolor('BuDRd_12'))
for i=1:1
    clf
    depth = squeeze(hindcast.depth(i,:,latidx,:));
    PV = abs(double(squeeze(hindcast.PV(i,:,latidx,:))));
    depth = [depth depth];
    PV = [PV PV];
    dens = squeeze(mean(pden(:,:,latidx,:)));
    dens = [dens dens];
    for j = 1:length(sectorlons)
        
        lonidx = longrid(1,:) > min(sectorlons{j}) & ...
            longrid(1,:) < max(sectorlons{j});
        subplot(1,3,j)
        hold on
        
        %         PVavg = nanmean(abs(double(squeeze(hindcast.PV(i,:,latavg,:)))),2);
        %         PVavg = squeeze(PVavg);
        %         PVavg = [PVavg PVavg];
        %         size(PVavg)
        %         break
        %         PVplot = PV - repmat(nanmean(PV(:,lonidx),2),[1 720]);
        %         contourf(longrid(:,lonidx),depth(:,lonidx),PVplot(:,lonidx),levels*.1,'LineColor','None')
        pcolor(longrid(:,lonidx),depth(:,lonidx),PV(:,lonidx))
        shading flat;
        load([inpath.sigmatheta denfiles(i).name]);
        
        
        [cs v] =contour(longrid(:,lonidx),depth(:,lonidx),dens(:,lonidx)-1000,25:.2:28,'LineColor','Black');
        clabel(cs,v);
        set(gca,'ydir','reverse')
        xlim(sectorlons{j})
        ylim([0 1500])
        caxis([1e-11 1e-9]*0.5)
        xlabel(sectortitle{j})
        colorbar('SouthOutside')
    end
    %     drawnow; pause(1)
end

%% Meridional Transects with salinity
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


for i=1
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
    contourf(latgrid,depth,log10(PV),-11:0.1:-9,'LineColor','None')
    
    load([inpath.sigmatheta denfiles(i).name]);
    
    dens = squeeze(mean(pden(:,:,:,lonidx)));
    [cs v] =contour(latgrid,depth,dens-1000,25:.2:28,'LineColor','Black');
    clabel(cs,v);
    set(gca,'ydir','reverse')
    xlim([-60 -10])
    ylim([0 1500])
    caxis(log10([1e-2 1]*5e-10))
    caxis([-10.2 -9.5])
    
    cax=  colorbar('SouthOutside');
    xlabel('log10(PV)')
    title(sprintf('%d %s: %f',i+1947,wocenames{lineidx},lon));
    
    subplot(1,2,2)
    hold on
    contourf(latgrid,depth,salt,33:0.05:35,'LineColor','None')
    
    load([inpath.sigmatheta denfiles(i).name]);
    
    dens = squeeze(mean(pden(:,:,:,lonidx)));
    [cs v] =contour(latgrid,depth,dens-1000,25:.2:28,'LineColor','Black');
    clabel(cs,v);
    set(gca,'ydir','reverse')
    xlim([-60 -10])
    ylim([0 1500])
%     caxis(log10([1e-2 1]*5e-10))
%     caxis([-10.2 -9.5])
    
    cax=  colorbar('SouthOutside');
    xlabel('Salinity')
    title(sprintf('%d %s: %f',i+1947,wocenames{lineidx},lon));
    caxis([33.9 34.5])
end


