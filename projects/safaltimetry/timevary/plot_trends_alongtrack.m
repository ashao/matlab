trackpath = 'C:\Users\ashao\Data\vxxc\';

%%
calculate_trend;
ntracks = length(saf.tracknum);

%%
clf
axesm('MapProjection','stereo','MapLonLimit',[0 360], ...
    'MapLatLimit',[-90 -40],'Frame','On','Grid', 'On','MeridianLabel','On', ...
    'ParallelLabel','On','PLineLocation',-90:10:-40,'MLineLocation',0:60:360, ...
    'LabelRotation','on','MLabelParallel',-50,'MlineException',0)
colormap(othercolor('BuDRd_12'))
whitebg(gcf,[0.9 0.9 0.9])
%%
% plot_fronts;
% hold on;
for t = 1:ntracks
    fprintf('Track %d/%d\n', t,ntracks)
    load([trackpath sprintf('t%03d.mat',saf.tracknum(t))]);
    
    saf.north.meanlat(t) = mean(saf.north.lat(:,t));
    saf.north.meanlon(t) = interp1(track.lat,track.lon,saf.north.meanlat(t));
    
    saf.center.meanlat(t) = mean(saf.center.lat(:,t));
    saf.center.meanlon(t) = interp1(track.lat,track.lon,saf.center.meanlat(t));
    
    saf.south.meanlat(t) = mean(saf.south.lat(:,t));
    saf.south.meanlon(t) = interp1(track.lat,track.lon,saf.south.meanlat(t));
%     scatterm(meanlat,meanlon,30,saf.north.trend(t),'filled')
    
%     plotm(track.lat,track.lon)
end
%%
clf
extent = {'north','center','south'};
for i=1:length(extent)
subplot(1,3,i)
axesm('MapProjection','stereo','MapLonLimit',[-180 180], ...
    'MapLatLimit',[-90 -30],'Frame','On','Grid', 'On','MeridianLabel','On', ...
    'ParallelLabel','On','PLineLocation',-90:15:-40,'MLineLocation',-180:60:180, ...
    'LabelRotation','on','MLabelParallel',-35,'MlineException',0)
colormap(othercolor('BuDRd_12'))
whitebg(gcf,[0.9 0.9 0.9])
scatterm(saf.(extent{i}).meanlat,saf.(extent{i}).meanlon,30,saf.(extent{i}).trend*110*365,'filled');
geoshow('landareas.shp', 'FaceColor', [0.15 0.5 0.15])
colorbar;
caxis([-1 1]*5e-4*110*365)
title(sprintf('%s Extent',extent{i}))
end
%%
clf
for i=1:length(extent)
[sortmeanlon sortidx] = sort(saf.(extent{i}).meanlon);
sortmeanlat = saf.(extent{i}).meanlat(sortidx);
% idx = saf.north.trend>0;

subplot(3,1,i)
area(sortmeanlon,smooth(saf.(extent{i}).trend,10,'mean')*110*365)
ylim([-10 10])
grid on
xlabel('Longitude')
ylabel('Trend (km/yr)')
end