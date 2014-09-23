inpath = '/ltraid3/ashao/uw-apl/projects/saf_altimetry/longterm/range/processed/';
files = dir([inpath 't*.mat']);
nfiles = length(files)
trackpath = '/ltraid4/aviso/alongtrack/sla/vxxc_matlab/groundpath.mat';
load(trackpath);
%%


centers = zeros(nfiles*4,4);

sidx = 1;
for tidx = 1:nfiles
    
    fprintf('Track %03d\n',tidx)
    load([inpath files(tidx).name])
    optpar = cell2mat(opt_track.optpar');
    [ntrans npar] = size(optpar);
    idx = sidx:(sidx+ntrans-1);
    centers(idx,1) = tidx;
    centers(idx,2) = optpar(:,2);
    centers(idx,3) = interp1(groundpath(tidx).lat,groundpath(tidx).lon,optpar(:,2));
    centers(idx,4) = sqrt(optpar(1).^2 + 2*optpar(3).^2);
    centers(idx,5) = cell2mat(opt_track.R2);
    sidx = sidx+ntrans;
    
end
centers(sidx:end,:) = [];
%%
addpath('/ltraid4/topography/')
[latgrat longrat Z] = satbath(1,[-90 -30],[0 359.9]);
Z(Z>0) = 0;
[slope null gradN, gradE] = gradientm(latgrat,longrat,sw_f(latgrat)./Z);
coast = load('coast.mat');
%%
load ~/matlab/projects/safaltimetry/sokolov/sokolov.mat
clf
notnan = ~isnan(centers(:,3)) & centers(:,5) > 0.5;
% axesm('MapProjection','Stereo','MapLatLimit',[-90 -30],'MapLonLimit',[0 360],'Frame','On','Grid','On');
worldmap([-65 -40],[60 90])
pcolorm(latgrat,longrat,Z)
caxis([-5000 -500])
scatterm(centers(notnan,2),centers(notnan,3),30,'k','filled','LineWidth',1)
% plotm(coast.lat,coast.long,'k','LineWidth',1)
plotm(sokolov.pfm(:,2),sokolov.pfm(:,1),'k','LineWidth',2)
plotm(sokolov.safm(:,2),sokolov.safm(:,1),'k','LineWidth',2)
% caxis([0 1]*100)
tightmap
colorbar
%% Kerguelen
load ~/matlab/projects/safaltimetry/sokolov/sokolov.mat
figure(1);clf
notnan = ~isnan(centers(:,3)) & centers(:,5) > 0.7;
% axesm('MapProjection','Stereo','MapLatLimit',[-90 -30],'MapLonLimit',[0 360],'Frame','On','Grid','On');
worldmap([-65 -40],[60 90])
pcolorm(latgrat,longrat,Z)
caxis([-5000 -500])
scatterm(centers(notnan,2),centers(notnan,3),30,'k','filled','LineWidth',1)
% plotm(coast.lat,coast.long,'k','LineWidth',1)
plotm(sokolov.pfm(:,2),sokolov.pfm(:,1),'k','LineWidth',2)
plotm(sokolov.safm(:,2),sokolov.safm(:,1),'k','LineWidth',2)
% caxis([0 1]*100)
tightmap
colorbar
%%
load ~/matlab/projects/safaltimetry/sokolov/sokolov.mat
figure(2);clf
notnan = ~isnan(centers(:,3)) & centers(:,5) > 0.8;
% axesm('MapProjection','Stereo','MapLatLimit',[-90 -30],'MapLonLimit',[0 360],'Frame','On','Grid','On');
worldmap([-70 -40],[0 90])
pcolorm(latgrat,longrat,Z)
caxis([-5000 -500])
scatterm(centers(notnan,2),centers(notnan,3),30,'k','filled','LineWidth',1)
% plotm(coast.lat,coast.long,'k','LineWidth',1)
plotm(sokolov.pfm(:,2),sokolov.pfm(:,1),'k','LineWidth',2)
plotm(sokolov.safm(:,2),sokolov.safm(:,1),'k','LineWidth',2)
% caxis([0 1]*100)
tightmap
colorbar
%%
load ~/matlab/projects/safaltimetry/sokolov/sokolov.mat
figure(3);clf
notnan = ~isnan(centers(:,3)) & centers(:,5) > 0.7;
% axesm('MapProjection','Stereo','MapLatLimit',[-90 -30],'MapLonLimit',[0 360],'Frame','On','Grid','On');
worldmap([-70 -40],[180 270])
pcolorm(latgrat,longrat,Z)
caxis([-5000 -500])
scatterm(centers(notnan,2),centers(notnan,3),30,'k','filled','LineWidth',1)
% plotm(coast.lat,coast.long,'k','LineWidth',1)
plotm(sokolov.pfm(:,2),sokolov.pfm(:,1),'k','LineWidth',2)
plotm(sokolov.safm(:,2),sokolov.safm(:,1),'k','LineWidth',2)
% caxis([0 1]*100)
tightmap
colorbar
%%
load ~/matlab/projects/safaltimetry/sokolov/sokolov.mat
figure(4);clf
notnan = ~isnan(centers(:,3)) & centers(:,5) > 0.8;
% axesm('MapProjection','Stereo','MapLatLimit',[-90 -30],'MapLonLimit',[0 360],'Frame','On','Grid','On');
worldmap([-70 -40],[270 360])
pcolorm(latgrat,longrat,Z)
caxis([-5000 -500])
scatterm(centers(notnan,2),centers(notnan,3),30,'k','filled','LineWidth',1)
plotm(coast.lat,coast.long,'k','LineWidth',1)
plotm(sokolov.pfm(:,2),sokolov.pfm(:,1),'k','LineWidth',2)
plotm(sokolov.safm(:,2),sokolov.safm(:,1),'k','LineWidth',2)
% caxis([0 1]*100)
tightmap
colorbar

%%
load ~/matlab/projects/safaltimetry/sokolov/sokolov.mat
figure(5);clf
notnan = ~isnan(centers(:,3)) & centers(:,5) > 0.7;
% axesm('MapProjection','Stereo','MapLatLimit',[-90 -30],'MapLonLimit',[0 360],'Frame','On','Grid','On');
worldmap([-70 -40],[90 180])
pcolorm(latgrat,longrat,Z)
caxis([-5000 -500])
scatterm(centers(notnan,2),centers(notnan,3),30,'k','filled','LineWidth',1)
plotm(coast.lat,coast.long,'k','LineWidth',1)
plotm(sokolov.pfm(:,2),sokolov.pfm(:,1),'k','LineWidth',2)
plotm(sokolov.safm(:,2),sokolov.safm(:,1),'k','LineWidth',2)
% caxis([0 1]*100)
tightmap
colorbar