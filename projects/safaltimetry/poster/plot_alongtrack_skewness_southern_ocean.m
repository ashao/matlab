inpath = 'C:\Users\ashao\Data\vxxc_matlab';
files = dir([inpath filesep 't*.mat']);
nfiles = length(files);
ntracks = 1;
%%
clf
minlat = -90;
maxlat = -20;
mapprops = {'stereo','MapLonLimit',[-180 180], ...
    'MapLatLimit',[minlat maxlat],'Frame','On','Grid', 'off','MeridianLabel','Off', ...
    'ParallelLabel','Off','PLineLocation',[-60 -30],'MLineLocation',-180:60:120, ...
    'MLabelParallel',maxlat-2,'LabelRotation','on'};
axesm(mapprops{:})
coast = load('coast.mat');
set(gca,'Color',[1 1 1]*1);
caxis([-0.5 0.5])
% worldmap

for tidx = 1:3:254
    fprintf('Track %d\n',tidx)
    load([inpath filesep files(tidx).name])
    plotidx = track.lat <= maxlat-5;
    plotlat = track.lat(plotidx);
    
    plotlon = track.lon(plotidx);
    plotskew = skewness(track.sla(:,plotidx));
    nplotpts = length(plotlon);
    
    truncidx = 1:5:nplotpts;
    scatterm(plotlat(truncidx),plotlon(truncidx),2,plotskew(truncidx),'filled')      
end
gridm; plabel; mlabel;
geoshow('landareas.shp', 'FaceColor', [1 1 1]-0.1,'LineWidth',2)
set(gca,'Color',[1 1 1]*0.6)
%%
outpath = 'C:\Users\ashao\Documents\GitHub\matlab\projects\safaltimetry\poster\';
export_fig([outpath 'skewness_global.eps'],'-eps','-cmyk')