function [ ] = plot_fronts( lineprops, saf_path  );
if nargin<2
    saf_path = '/home/ashao//matlab/projects/safaltimetry/saf.mat';
end
if nargin < 1
    lineprops = {'k.','MarkerSize',0.5};
end
load /home/ashao/uw-apl/data/mimoc/mimoc.ctemp.salt.mldepth.mat
load(saf_path);
hold on;

padidx = 1;
axesm('MapProjection','mercator','MapLonLimit',[0 360], ...
    'MapLatLimit',[-65 -40],'Frame','On','Grid', 'On','MeridianLabel','On', ...
    'ParallelLabel','On','PLineLocation',-90:5:-40,'MLineLocation',0:60:360,'LabelRotation','on','MLabelParallel',-50,'MlineException',0)
% m_pcolor(mimoc.lon,mimoc.lat,max(mimoc.mldepth));shading flat;
plotm(saf.sokolov.sam_m.lat(padidx:end-padidx),saf.sokolov.sam_m.lon(padidx:end-padidx),lineprops{:});
plotm(saf.sokolov.sam_n.lat(padidx:end-padidx),saf.sokolov.sam_n.lon(padidx:end-padidx),lineprops{:});
plotm(saf.sokolov.sam_s.lat(padidx:end-padidx),saf.sokolov.sam_s.lon(padidx:end-padidx),lineprops{:});
plotm(saf.orsi.lon,saf.orsi.lat,lineprops{:})
tightmap;

