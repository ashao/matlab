climfile = '/home/ashao/uw-apl/development/Offtrac/branches/ttd_bp/ttd.normalyear.0730.nc';
hindfile = '/home/ashao/uw-apl/development/Offtrac/branches/ttd_bp/ttd.hindcast.0730.nc';
load metrics
%% Truncate to just the Indian Ocean
lonidx(1) = min( find(metrics.lonh.data>20) );
lonidx(2) = min( find(metrics.lonh.data>-230 & ...
    metrics.lonh.data < 0) );
inlat = metrics.lath.data < 20 & metrics.lath.data > -60;
latidx = min(find(inlat));
nlat = sum(inlat);
%%
nlon = length(lonidx(1):length(metrics.lonh.data)) + length(1:lonidx(2))-1;
ntime = 720;
nlay = 49;
lonrange1 = 1:(length(lonidx(1):360));
lonrange2 = (lonrange1(end)+1);
lonrange2 = lonrange2:(lonrange2+length(1:lonidx(2)-2));
%%
ttd.clim = zeros(ntime,nlay,nlat,nlon);
ttd.hind = zeros(ntime,nlay,nlat,nlon);

start1 = [ 0 0 latidx-1 lonidx(1)-1 ];
count1 = [ -1 -1 nlat inf ];
start2 = [ 0 0 latidx-1 0 ];
count2 = [-1 -1 nlat lonidx(2)-1 ];

ttd.clim(:,:,:,lonrange1) = nc_varget(climfile,'mn_ttd',start1,count1);
ttd.clim(:,:,:,lonrange2) = nc_varget(climfile,'mn_ttd',start2,count2);

ttd.hind(:,:,:,lonrange1) = nc_varget(hindfile,'mn_ttd',start1,count1);
ttd.hind(:,:,:,lonrange2) = nc_varget(hindfile,'mn_ttd',start2,count2);

npts = nlat*nlon;
ttd.similarity = zeros(nlay,nlat,nlon);
for lay = 1:nlay
    fprintf('Layer %d\n',lay)
    for pt = 1:npts
        R = corrcoef(squeeze(ttd.clim(:,lay,pt)),squeeze(ttd.hind(:,lay,pt)));
        ttd.similarity(lay,pt) = R(2,1);
        
    end
end

%%
save -v7.3 ~/uw-apl/projects/him_hindcasts/indian_ocean_ttd.mat ttd
%%
load ~/uw-apl/projects/him_hindcasts/indian_ocean_ttd.mat ttd
load /home/ashao/uw-apl/models/HIM/hindcast/hindcast.PV.mat
%%
ttd.lat = nc_varget(climfile,'lath',latidx-1,nlat);
ttd.lon(lonrange1) = nc_varget(climfile,'lonh',lonidx(1)-1,inf);
ttd.lon(lonrange2) = nc_varget(climfile,'lonh',0,lonidx(2)-1);
ttd.lon = wrapTo360(ttd.lon);
%% Plot along isopycnals
coast = load('coast.mat');
colormap(flipud(othercolor('Blues9')))
for lay = 1:nlay
    clf
    worldmap([min(ttd.lat) max(ttd.lat)],[min(ttd.lon) max(ttd.lon)]);
    contourfm(ttd.lat,ttd.lon,squeeze(ttd.similarity(lay,:,:)),0.1:0.05:1,'LineStyle','none')
    [cs v] = contourm(metrics.geolat.data,metrics.geolon.data,double(squeeze(mean(hindcast.h(:,lay,:,:)))),10:50:200,'LineColor','black','LineWidth',2);
    clabelm(cs,v)
    linem(coast.lat,coast.long,'LineWidth',2','Color','black')    
    title(sprintf('Layer %02d',lay))
    colorbar
    caxis([0.5 1])
    pause(0.1)
    
end
%% Plot a transect along I8
nlay = 49;
figure;
subplot(2,1,1); hold on;
lonplot = 60.5;
lonidx = metrics.lonh.data == lonplot;
depths = squeeze(mean(hindcast.depth(:,1:nlay,:,lonidx)));
latgrid = double(repmat(metrics.lath.data(:)',[nlay 1]));
contourf(latgrid,depths,squeeze(ttd.similarity(:,:,lonidx)),0:0.05:1,'LineStyle','none')
% contour(latgrid,depths,squeeze(log10(mean(hindcast.PV(:,1:nlay,inlat,lonidx_hind)))),[-10 10],'LineColor','black');
colorbar
caxis([0.5 1])
ylim([0 1500]);xlim([-65 20])
set(gca,'ydir','reverse')
title(sprintf('Transect along %2.1f',lonplot))
colormap(flipud(othercolor('Blues9')))

subplot(2,1,2)
lonplot = -270.5; hold on;
lonidx = metrics.lonh.data ==lonplot;
% lonidx_hind = find (abs(metrics.lonh.data - lonplot) < 0.1);
depths = squeeze(mean(hindcast.depth(:,1:nlay,:,lonidx)));
latgrid = double(repmat(metrics.lath.data(:)',[nlay 1]));
contourf(latgrid,depths,squeeze(ttd.similarity(:,:,lonidx)),0:0.05:1,'LineStyle','none')
% contour(latgrid,depths,squeeze(log10(mean(hindcast.PV(:,1:nlay,inlat,lonidx_hind)))),[-10 10],'LineColor','black');
colorbar
ylim([0 1500]);xlim([-65 20])
caxis([0.5 1])
set(gca,'ydir','reverse')
title(sprintf('Transect along %2.1f',lonplot))
colormap(flipud(othercolor('Blues9')))
