inpath = 'C:\Users\ashao\data\vxxc\';
outpath = 'C:\Users\ashao\altimetry\longterm\';
if ~exist(outpath)
    mkdir(outpath)
end

files = dir([inpath '*.mat']);
nfiles = length(files);

%%
minlat = -70;
maxlat = -30;
mapprops = {'MapProjection','eqdcyl','MapLonLimit',[-180 180], ...
        'MapLatLimit',[minlat maxlat],'Frame','On','Grid', 'off','MeridianLabel','Off', ...
        'ParallelLabel','Off','PLineLocation',[-90:10:maxlat],'MLineLocation',-180:60:120, ...
        'MLabelParallel',maxlat-2,'LabelRotation','on'};
startt = 1;    
for t = startt:nfiles
    clf;
    
    load([inpath files(t).name])
%     
%     
%     subplot(2,1,2); hold on;
%     skew = skewness(track.sla,0,2);
%     plot(track.lat,skew,'k');
%     plot(track.lat,smooth(skew,10,'mean'),'k','LineWidth',2);
%     xlim([-65 maxlat]);
%     ylim([-1 1])
%     grid on    
    
%     subplot(2,1,1)
    axesm(mapprops{:});
    geoshow('landareas.shp', 'FaceColor', [1 1 1]-0.1)
    plot_fronts;
    plotm(track.lat,track.lon,'k','LineWidth',2);
    plabel;mlabel;gridm;tightmap;    
    daspect([3 1 1])
    [track.latrange,null] = inputm(2);
    fprintf('Track %d Latrange: %f %f\n',t,track.latrange(1),track.latrange(2))
    save([outpath files(t).name],'track');
end