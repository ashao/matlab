%%

inpath.tracks = 'C:\Users\ashao\data\aviso\';
files = dir([inpath.tracks 't*.mat']);
ntracks = length(files);

%%
map_proj = { 'MapProjection','eqdazim','MapLonLimit',[0 360], ...
    'MapLatLimit',[-90 -30],'Frame','on','Grid', 'On','MeridianLabel','On', ...
    'ParallelLabel','On','MLineLocation',-180:60:179, ...
    'LabelRotation','On','MLabelParallel',-35,'MlineException',0, ...
    'PLabelMeridian','prime', 'FontName','MyriadPro-Regular'} ;
%%
% h=figure('Visible','Off');
figure;
cla;clf;
axesm(map_proj{:})
colormap(othercolor('BuOrR_14'))
for track = 1:10
    
    if mod(track-1,10)==0
       fprintf('Track %d/%d\n',track,ntracks) 
    end
    load([inpath.tracks files(track).name])
    skew = skewness(track.sla,1,2);
    try 
        scatterm(track.lat,track.lon,1,skew,'filled')
    catch
        fprintf('Error in track %d\n',track.tracknum);
    end
    caxis([-1 1])    
end
cax = colorbar('SouthOutside');
xlabel('Skewness')
geoshow('landareas.shp','FaceColor',[0 0 0]+1)
% gridm;plabel;mlabel;framem;