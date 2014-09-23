%%
climfile = '/home/ashao/uw-apl/development/Offtrac/branches/ttd_bp/ttd.normalyear.annual.0730.nc';
hindfile = '/home/ashao/uw-apl/development/Offtrac/branches/ttd_bp/ttd.hindcast.annual.0730.nc';
load metrics
%%
ttd.clim = nc_varget(climfile,'mn_ttd');
ttd.hind = nc_varget(hindfile,'mn_ttd');
%%
nlay = 49;
wetmask = logical(metrics.wet.data);
ptlist = find(wetmask(:));
ttd.similarity = nan([49 size(wetmask)]);
   
for lay = 1:nlay
    fprintf('Layer %d\n',lay)
    for pt = 1:length(ptlist)
        ptidx = ptlist(pt);
        R = corrcoef(squeeze(ttd.clim(:,lay,ptidx)), ...
            squeeze(ttd.hind(:,lay,ptidx)));
        ttd.similarity(lay,ptidx) = R(1,2);
    end
end
%%
outpath = '~/uw-apl/figs/him_hindcasts/ttd_comparison/layer/';
mkdir(outpath)
coast = load('coast.mat');
figure;
for lay = 26
    clf;hold on;
    worldmap([-60 30],[20 140])
    contourfm(metrics.geolat.data,metrics.geolon.data, ...
        squeeze(ttd.similarity(lay,:,:)),0:0.1:1,'LineStyle','None')
    plotm(coast.lat,coast.long,'k')
    caxis([0.5 1])
    colorbar
    title(sprintf('Layer %d',lay))
    pause(0.1)
%     saveas(gcf,[outpath sprintf('ttd_compare_layer_%02d.eps',lay)],'epsc')
end
%% 
figure
latidx = metrics.lath.data == -32.5;
lonidx = metrics.lonh.data ==-270.5;
subplot(1,2,1); hold on;
plot(ttd.clim(:,20,latidx,lonidx)/sum(ttd.clim(:,20,latidx,lonidx)),'k--')
plot(ttd.hind(:,20,latidx,lonidx)/sum(ttd.hind(:,20,latidx,lonidx)),'k')
xlabel('Year since BP')
title('SAMW 32S 90E')
% legend('Normal year','Hindcast')
axis square
subplot(1,2,2); hold on;
plot(1:60,ttd.clim(:,26,latidx,lonidx),'k--')
plot(ttd.hind(:,26,latidx,lonidx),'k')
legend('Normal year','Hindcast','Location','South')
axis square
xlabel('Year since BP')
title('AAIW 32S 90E')