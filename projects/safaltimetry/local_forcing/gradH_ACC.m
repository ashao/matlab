addpath('/ltraid4/topography/')
[latgrat longrat Z] = satbath(20,[-90 -30],[0 359.9]);
Z(Z>0) = 0;
[slope null gradN, gradE] = gradientm(latgrat,longrat,sw_f(latgrat)./Z);
% load ~/uw-apl/projects/saf_altimetry/timevary.mat
%%
worldmap([-90 -30],[0 360])
pcolorm(latgrat,longrat,abs(gradN));clf

% 
% %%
% extents = {'north','mean','south'};
% %%
% for i = [1 3]
%     
%     extent = extents{i};
%     scatterm(nanmean(timevary.(extent).lat,2),nanmean(timevary.(extent).lon,2),80,'black','filled')
%     scatterm(nanmean(timevary.(extent).lat,2),nanmean(timevary.(extent).lon,2),20,'yellow','filled')
%     
% end