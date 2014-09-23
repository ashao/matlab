inpath = '/ltraid3/ashao/uw-apl/projects/saf_altimetry/longterm/curvefit/processed2/';
trackpath = '/ltraid3/ashao/uw-apl/projects/saf_altimetry/bootstrap_skewness/';
files = dir([inpath 't*.mat']);
nfiles = length(files)


%% Plot the centers as a scatter plot

temp.lat = nan(nfiles,1);
temp.lon = nan(nfiles,1);
temp.width = nan(nfiles,1);

acc.north = temp;
acc.mean = temp;
acc.south = temp;

for tidx = 1:nfiles
    
    load([inpath files(tidx).name])
    load([trackpath files(tidx).name])
    optpar = cell2mat(opt_track.optpar');
    R2 = cell2mat(opt_track.R2);
    skews = cellfun(@(x)(max(abs(x))),opt_track.skewness);
    
    nskews = length(skews);
    skewidx = zeros(nskews,1);
    
    for i = 1:nskews
        skewidx(i) = find(abs(track.skewness)==skews(i));
    end
    
    sigidx  = track.sigidx(skewidx);
    
    plotidx = find(R2 < 0.2 & optpar(:,4)' > 0.55 & sigidx');
    
    
    nplot = length(plotidx);
    fprintf('Plotting Track %d: %d centers\n',tidx,nplot)
    if ~isempty(plotidx)
        [acc.mean.lat(tidx) acc.mean.lon(tidx)] = meanm(opt_track.processed.lat(plotidx),opt_track.processed.lon(plotidx));
        acc.mean.width(tidx) = mean(opt_track.processed.width(plotidx));
        if length(plotidx) > 1
            % Northern extent
            [maxlat maxidx] = max(opt_track.processed.lat(plotidx));
            maxidx = plotidx(maxidx);                        
            acc.north.lat(tidx) = maxlat;
            acc.north.lon(tidx) = opt_track.processed.lon(maxidx);
            acc.north.width(tidx) = opt_track.processed.width(maxidx);
            % Average of all transitions
            
            % Southern extent
            [minlat minidx] = min(opt_track.processed.lat(plotidx));
            minidx = plotidx(minidx);
            minlon = opt_track.processed.lon(minidx);
            acc.south.lat(tidx) = minlat;
            acc.south.lon(tidx) = minlon;
            acc.south.width(tidx) = opt_track.processed.width(minidx);
        end
    end
end

nanidx = isnan(acc.south.lat) | isnan(acc.south.lon);
acc.south.lat(nanidx) = [];
acc.south.lon(nanidx) = [];
acc.mean.lon = wrapTo360(acc.mean.lon);
% scatterm(acc.north.lat,acc.north.lon,'bo','filled')
% scatterm(acc.mean.lat,acc.mean.lon,'ko','filled');
% scatterm(acc.south.lat,acc.south.lon,'ro','filled')
%%
coast=load('coast.mat');
clf
map_proj = { 'MapProjection','stereo','MapLonLimit',[-180 180], ...
    'MapLatLimit',[-90 -30],'Frame','off','Grid', 'off','MeridianLabel','off', ...
    'ParallelLabel','Off','MLineLocation',-180:60:179, ...
    'LabelRotation','On','MLabelParallel',-31,'MlineException',0, ...
    'PLabelMeridian','prime','PLineLocation',-70:10:-40} ;
axesm(map_proj{:})
colors = {'r','k','b'};
extents = {'north','mean','south'};
plotlon = 0:0.5:360;
cla

for i = 1:length(extents)
    centers = sortrows([acc.(extents{i}).lat acc.(extents{i}).lon],2);
    notnan = ~isnan(centers(:,1)) & ~isnan(centers(:,2));
    centers = centers(notnan,:);
    centers(end+1,:) = centers(1,:);
    centers(end,2) = centers(end,2)+360;
    pp = splinefit(centers(:,2),centers(:,1),30,'r');
    plotlat = ppval(pp,plotlon);
    acc.(extents{i}).spline.lat = plotlat;
    acc.(extents{i}).spline.lon = plotlon;
    plotm(plotlat,plotlon,colors{i},'LineWidth',2)
    scatterm(centers(:,1),centers(:,2),colors{i},'filled')
end
%
plotm(coast.lat,coast.long,'k','LineWidth',2);
framem; gridm; plabel; mlabel; tightmap;



%% Plot jet centers aound the Kerguelen plateau and Drake Passage
clf
% axesm('MapProjection','Stereo','MapLatLimit',[-65 -30],'MapLonLimit',[0 360],'Frame','On','Grid','On');
latlim = [-70 -35]; lonlim = [-90 -10];
worldmap(latlim,lonlim)
[latgrat longrat elev] = satbath(20,latlim,lonlim);
contourfm(latgrat,longrat,-elev,0:250:5000,'Color','none')
for tidx = 1:1:254
    fprintf('Track %d/%d\n',tidx,nfiles)
    load([inpath files(tidx).name])
    load([trackpath files(tidx).name])
    optpar = cell2mat(opt_track.optpar');
    R2 = cell2mat(opt_track.R2);
    lons = wrapTo180(opt_track.processed.lon);
    plotidx = R2 < 0.1 & optpar(:,4)' > 0.6 & lons > min(lonlim) & ...
        lons < max(lonlim);
    if sum(plotidx)>0
        sum(plotidx)
        % Northern extentR2 < 0.2 & optpar(:,4)' > 0.6 &
        plotm(opt_track.processed.lat(plotidx), ...
            lons(plotidx),'kx','LineWidth',2,'MarkerSize',10)
        drawnow
    end
    
end

load ACCfronts_mean
front_idx = 1:12
% colors = {'r','k','b'};
for i = 1:length(front_idx)
    plotm(Fronts.y{ front_idx(i) },Fronts.x{ front_idx(i) },'k','LineWidth',1.5)
end
cax = colorbar('SouthOutside');
xlabel(cax,'Bathymetry [m]')
geoshow('landareas.shp', 'FaceColor', [1 1 1]);

%% Plot jet centers aound the Kerguelen plateau
clf
% axesm('MapProjection','Stereo','MapLatLimit',[-65 -30],'MapLonLimit',[0 360],'Frame','On','Grid','On');
latlim = [-65 -40]; lonlim = [40 90];
worldmap(latlim,lonlim)
[latgrat longrat elev] = satbath(20,latlim,lonlim);
contourfm(latgrat,longrat,-elev,0:250:5000,'Color','none')
for tidx = 1:1:254
    fprintf('Track %d/%d\n',tidx,nfiles)
    load([inpath files(tidx).name])
    load([trackpath files(tidx).name])
    optpar = cell2mat(opt_track.optpar');
    R2 = cell2mat(opt_track.R2);
    lons = wrapTo360(opt_track.processed.lon);
    plotidx = R2 < 0.2 & optpar(:,4)' > 0.6 & lons > min(lonlim) & ...
        lons < max(lonlim);
    if sum(plotidx)>0
        sum(plotidx)
        % Northern extentR2 < 0.2 & optpar(:,4)' > 0.6 &
        plotm(opt_track.processed.lat(plotidx), ...
            lons(plotidx),'kx','LineWidth',2,'MarkerSize',10)
        drawnow
    end
    
end

load ACCfronts_mean
front_idx = 1:12
% colors = {'r','k','b'};
for i = 1:length(front_idx)
    plotm(Fronts.y{ front_idx(i) },Fronts.x{ front_idx(i) },'k','LineWidth',1.5)
end
cax = colorbar('SouthOutside');
xlabel(cax,'Bathymetry [m]')
geoshow('landareas.shp', 'FaceColor', [1 1 1]);

%%

%% Plot jet centers aound the Kerguelen plateau
clf
% axesm('MapProjection','Stereo','MapLatLimit',[-65 -30],'MapLonLimit',[0 360],'Frame','On','Grid','On');
latlim = [-90 -40]; lonlim = [0 360];
worldmap(latlim,lonlim)
% [latgrat longrat elev] = satbath(20,latlim,lonlim);
for tidx = 1:1:254
    fprintf('Track %d/%d\n',tidx,nfiles)
    load([inpath files(tidx).name])
    load([trackpath files(tidx).name])
    optpar = cell2mat(opt_track.optpar');
    R2 = cell2mat(opt_track.R2);
    lons = wrapTo360(opt_track.processed.lon);
    plotidx = R2 < 0.2 & optpar(:,4)' > 0.6 & lons > min(lonlim) & ...
        lons < max(lonlim);
    if sum(plotidx)>0
        sum(plotidx)
        % Northern extentR2 < 0.2 & optpar(:,4)' > 0.6 &
        plotm(opt_track.processed.lat(plotidx), ...
            lons(plotidx),'kx','LineWidth',2,'MarkerSize',10)
        drawnow
    end
    
end
% contourfm(latgrat,longrat,-elev,0:250:5000,'Color','none')
load ACCfronts_mean
front_idx = 1:12
% colors = {'r','k','b'};
for i = 1:length(front_idx)
    plotm(Fronts.y{ front_idx(i) },Fronts.x{ front_idx(i) },'k','LineWidth',2)
end
cax = colorbar('SouthOutside');
xlabel(cax,'Bathymetry [m]')
geoshow('landareas.shp', 'FaceColor', [1 1 1]);




%% Generate the countour for the northern boundary
coast=load('coast.mat');
clf
axesm('MapProjection','Stereo','MapLatLimit',[-90 -30],'MapLonLimit',[0 360],'Frame','On','Grid','On');
for tidx = 1:nfiles
    
    load([inpath files(tidx).name])
    load([trackpath files(tidx).name])
    
    % Northern extent
    [maxlat maxidx] = max(opt_track.processed.lat);
    scatterm(maxlat,opt_track.processed.lon(maxidx),'bo','filled')
end
plotm(coast.lat,coast.long,'k','LineWidth',2)
load ACCfronts_mean
front_idx = [3 5 8];
colors = {'r','k','b'};
for i = 1:length(front_idx)
    plotm(Fronts.y{ front_idx(i) },Fronts.x{ front_idx(i) },colors{i})
end

%%
[outlat outlon]= digitize_on_map;
saf.north.lat = [saf.north.lat outlat];
saf.north.lon = [saf.north.lon outlon]