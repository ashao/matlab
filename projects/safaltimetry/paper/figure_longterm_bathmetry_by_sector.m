inpath = '/ltraid3/ashao/uw-apl/projects/saf_altimetry/longterm/range/processed/';
files = dir([inpath 't*.mat']);
nfiles = length(files)
trackpath = '/ltraid4/aviso/alongtrack/sla/vxxc_matlab/groundpath.mat';
load(trackpath);
outpath = '/ltraid3/ashao/uw-apl/matlab/projects/safaltimetry/paper/figs/';
%%


centers = zeros(nfiles*4,4);

sidx = 1;
for tidx = 1:nfiles
    
    fprintf('Track %03d\n',tidx)
    load([inpath files(tidx).name])
    optpar = cell2mat(opt_track.optpar');
    
    mintracklat = min(groundpath(tidx).lat);
    goodidx = optpar(:,2) > mintracklat+2;    
    optpar = optpar(goodidx,:);
   
    [ntrans npar] = size(optpar);
    
    idx = sidx:(sidx+ntrans-1);
    centers(idx,1) = tidx;
    centers(idx,2) = optpar(:,2);
    centers(idx,3) = intrplon(groundpath(tidx).lat,groundpath(tidx).lon,optpar(:,2));
    centers(idx,4) = sqrt(optpar(:,1).^2 + 2*optpar(:,3).^2);
    centers(idx,5) = cell2mat(opt_track.R2(goodidx));
    sidx = sidx+ntrans;
    
end
centers(sidx:end,:) = [];
%%
addpath('/ltraid4/topography/')
[latgrat longrat Z] = satbath(30,[-90 -30],[0 359.95]);
Z(Z>0) = 0;
% [slope null gradN, gradE] = gradientm(latgrat,longrat,sw_f(latgrat)./Z);
coast = load('coast.mat');
lonranges = [-70 20 ; 20 110 ; 110 200; 200 290]
latranges = [-65 -40 ; -65 -40; -70 -40 ; -70 -45]

%%
binwidth = 5;
lonbins = 0.5:binwidth:359.5;
nbins = length(lonbins);

longterm.north.lon=zeros(nbins,1);
longterm.north.lat=zeros(nbins,1);
longterm.north.width=zeros(nbins,1);

longterm.south.lon=zeros(nbins,1);
longterm.south.lat=zeros(nbins,1);
longterm.south.width=zeros(nbins,1);

longterm.mean.lon=zeros(nbins,1);
longterm.mean.lat=zeros(nbins,1);
longterm.mean.width=zeros(nbins,1);

for bin = 1:length(lonbins)
    
    lonidx = centers(:,3) > (lonbins(bin)-binwidth/2) & ...
        centers(:,3) < (lonbins(bin) + binwidth/2) & ...
        centers(:,5) > 0.8;
    if any(lonidx)
        
        truncenter = centers(lonidx,:);
        
        [longterm.north.lat(bin) nidx] = max(truncenter(:,2));
        longterm.north.lon(bin) = truncenter(nidx,3);
        longterm.north.width(bin) = truncenter(nidx,4);
        
        [longterm.south.lat(bin) sidx] = min(truncenter(:,2));
        longterm.south.lon(bin) = truncenter(sidx,3);
        longterm.south.width(bin) = truncenter(sidx,4);
        
        [longterm.mean.lat(bin) longterm.mean.lon(bin)] = ...
            meanm(truncenter(:,2),truncenter(:,3));        
        longterm.mean.width(bin) = mean(truncenter(:,4));
        
    else
        longterm.mean.lat(bin) = NaN;
        longterm.mean.lon(bin) = NaN;
        longterm.north.lat(bin) = NaN;
        longterm.north.lon(bin) = NaN;
        longterm.south.lat(bin) = NaN;
        longterm.south.lon(bin) = NaN;
    end
    
end
        longterm.mean.lon = wrapTo360(longterm.mean.lon);
        longterm.north.lon = wrapTo360(longterm.north.lon);
        longterm.south.lon = wrapTo360(longterm.south.lon);

        longterm.mean.lat(bin+1) = longterm.mean.lat(1);
        longterm.mean.lon(bin+1) = longterm.mean.lon(1);
        longterm.mean.width(bin+1) = longterm.mean.width(1);
        longterm.north.lat(bin+1) = longterm.north.lat(1);
        longterm.north.lon(bin+1) = longterm.north.lon(1);
        longterm.north.width(bin+1) = longterm.north.width(1);
        longterm.south.lat(bin+1) = longterm.south.lat(1);
        longterm.south.lon(bin+1) = longterm.south.lon(1);
        longterm.south.width(bin+1) = longterm.south.width(1);

        %% Width
        figure; hold on;
plot(longterm.north.lon(1:end-1),smooth(longterm.north.width(1:end-1)*110,5,'mean'),'r-x','LineWidth',2)
% plot(longterm.mean.lon(1:end-1),longterm.mean.width(1:end-1)*110,'k','LineWidth',2)
plot(longterm.south.lon(1:end-1),smooth(longterm.south.width(1:end-1)*110,5,'mean'),'b-x','LineWidth',2)
grid on; legend('North','South');
xlim([0 360])
nanmean(longterm.north.width(1:end-1))*110
nanmean(longterm.south.width(1:end-1))*110
xlabel('Longitude (\circE)');
ylabel('Jet width (m)')
saveas(gcf,[outpath 'fig_longterm_width.eps'],'epsc')
%%
figure
worldmap([-90 -30],[0 360])
extents = {'north','mean','south'};
colors = {'r','k','b'};
load ~/matlab/projects/safaltimetry/sokolov/sokolov.mat
for i=1:3
    
    
    extent = extents{i};
    scatterm(longterm.(extent).lat,longterm.(extent).lon,30,colors{i},'filled')
%     plotm(longterm.(extent).lat,longterm.(extent).lon,colors{i},'LineWidth',2)
%     [longterm.(extent).curve.lat longterm.(extent).curve.lon] = digitize_on_map;

    
end
        plotm(sokolov.pfm(:,2),sokolov.pfm(:,1),'k--','LineWidth',1)
    plotm(sokolov.safm(:,2),sokolov.safm(:,1),'k--','LineWidth',1)
    contourm(latgrat,longrat,Z,[-2000 -2000],'LineColor','black');
%     colormap(othercolor('Blues5'))
geoshow('landareas.shp', 'FaceColor',  [0.5 0.7 0.5])
saveas(gcf,[outpath 'fig_extents_global.eps'],'epsc')
% print(gcf,[outpath filesep sprintf('fig_longterm_bathymetry_southern_ocean',i)],'-dpng','-r200')

%% Basin by basin

for i= 1:4
    subplot(2,2,i)
    load ~/matlab/projects/safaltimetry/sokolov/sokolov.mat
%     figure(i);clf
    notnan = ~isnan(centers(:,3)) & centers(:,5) > 0.8;
    worldmap(latranges(i,:),lonranges(i,:))
    contourfm(latgrat,longrat,-Z,0:250:6000,'LineColor','none')
%     shading flat
    caxis([0 6000])

   
    scatterm(centers(notnan,2),centers(notnan,3),30,'y','filled','LineWidth',1)
    plotm(sokolov.pfm(:,2),sokolov.pfm(:,1),'k','LineWidth',2)
    plotm(sokolov.safm(:,2),sokolov.safm(:,1),'k','LineWidth',2)
        for t=[1 3]
    
    
        extent = extents{t};
        scatterm(longterm.(extent).lat,longterm.(extent).lon,30, 'kx')
%     plotm(longterm.(extent).lat,longterm.(extent).lon,colors{i},'LineWidth',2)
%     [longterm.(extent).curve.lat longterm.(extent).curve.lon] = digitize_on_map;

    
    end
    tightmap;
    colormap(othercolor('Blues5'))
    h=colorbar('SouthOutside');
    xlabel(h,'Ocean Basin Depth (m)')
    geoshow('landareas.shp', 'FaceColor',  [0.5 0.7 0.5])
    saveas(gcf,[outpath 'fig_extents_basin.eps'],'epsc')
%     print(gcf,[outpath filesep sprintf('fig_longterm_bathymetry_%d',i)],'-dpng','-r200')
end