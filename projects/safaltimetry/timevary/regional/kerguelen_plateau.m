load('C:\Users\ashao\Data\allcenters.mat')
clf
latrange = [-70 -40];
lonrange = [55 95];
worldmap(latrange,lonrange)
[latgrat longrat Z] = satbath(1,latrange,lonrange);
% contourm(latgrat,longrat,Z,'LineColor','black','LineWidth',2);
pcolorm(latgrat,longrat,Z);
notnan = ~isnan(centers(:,3));
scatterm(centers(notnan,2),centers(notnan,3))