trackpath = '/ltraid4/aviso/alongtrack/sla/tracks/';
clf
axesm('MapProjection','stereo','MapLonLimit',[0 360], ...
    'MapLatLimit',[-90 -40],'Frame','On','Grid', 'On','MeridianLabel','On', ...
    'ParallelLabel','On','PLineLocation',-90:10:-40,'MLineLocation',0:60:360, ...
    'LabelRotation','on','MLabelParallel',-50,'MlineException',0)
trackfiles = dir([trackpath '*.mat']);
ntracks = length(trackfiles);
%%
hold on;
for t=1:3:ntracks
   
    load([trackpath trackfiles(t).name])
    plotm(track.lat,track.lon)
    
end

