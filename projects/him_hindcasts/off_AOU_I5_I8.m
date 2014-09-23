infile.oxygen = '/home/ashao/uw-apl/development/Offtrac/branches/o2_hindcast/oxygen.woa09.5542.nc';
infile.hfile = '/home/ashao/uw-apl/data/offtrac/input/normalyear/H-hind.nc';
infile.saltfile = '/home/ashao/uw-apl/data/offtrac/input/normalyear/TS-hind.nc';
load metrics
% Indices for end of June in each year (compare to Alvarez et al.)
years = [1987 1995 2002];
monindices = (years-1947)*12+1+5;
hindcast.oxygen = nan([length(monindices) 49 210 360]);
hindcast.o2sat = nan([length(monindices) 49 210 360]);
hindcast.depth = nan([length(monindices) 49 210 360]);
hindcast.salt = nan([length(monindices) 49 210 360]);

for t = 1:length(years)
    start4d = [monindices(t)-1 0 0 0];
    count4d = [1 -1 -1 -1];
    hindcast.oxygen(t,:,:,:) = (nc_varget(infile.oxygen,'mn_oxygen',start4d,count4d));
    hindcast.o2sat(t,:,:,:) = (nc_varget(infile.oxygen,'mn_o2sat',start4d,count4d));
    hindcast.h(t,:,:,:) = (nc_varget(infile.hfile,'h',start4d,count4d));
    hindcast.salt(t,:,:,:) = (nc_varget(infile.saltfile,'salt',start4d,count4d));
    hindcast.depth(t,:,:,:) = (cumsum(nc_varget(infile.hfile,'h',start4d,count4d)));
end
wetmask = logical(metrics.wet.data);
hindcast.oxygen(:,:,~wetmask)=NaN;
hindcast.o2sat(:,:,~wetmask)=NaN;
hindcast.depth(:,:,~wetmask)=NaN;

i5.latidx = find(metrics.lath.data == -32.5);
i5.lon = repmat(metrics.geolon.data(i5.latidx,:),[49,1]);
i5.oxygen = squeeze(hindcast.oxygen(:,:,i5.latidx,:));
i5.o2sat = squeeze(hindcast.o2sat(:,:,i5.latidx,:));
i5.aou = i5.o2sat-i5.oxygen;
i5.depth = squeeze(hindcast.depth(:,:,i5.latidx,:));
i5.h = squeeze(hindcast.h(:,:,i5.latidx,:));
i5.salt = squeeze(hindcast.salt(:,:,i5.latidx,:));

i8.lonidx = find(metrics.lonh.data ==-270.5);
i8.lat = repmat(metrics.geolat.data(:,i8.lonidx)',[49 1]);
i8.oxygen = squeeze(hindcast.oxygen(:,:,:,i8.lonidx));
i8.o2sat = squeeze(hindcast.o2sat(:,:,:,i8.lonidx));
i8.aou = i8.o2sat-i8.oxygen;
i8.depth = squeeze(hindcast.depth(:,:,:,i8.lonidx));
i8.h = squeeze(hindcast.h(:,:,:,i8.lonidx));
%%
colormap(flipud(othercolor('RdYlBu9')))
clf
subplot(3,1,1)
plotlon= [i5.lon i5.lon+360];
plotdepth = repmat(squeeze(mean(i5.depth(2:3,:,:))),[1 2]);
plotaou = repmat(squeeze(i5.aou(3,:,:)-i5.aou(2,:,:)),[1,2]);

% contourf(plotlon,plotdepth,plotaou,(-1:0.05:1)*0.02)
pcolor(plotlon,plotdepth,plotaou)
set(gca,'ydir','reverse')
ylim([0 1200]);xlim([30 120]); caxis([-1 1]*0.01); colorbar;
title('2002-1995')
shading flat;

subplot(3,1,2)
plotlon= [i5.lon i5.lon+360];
plotdepth = repmat(squeeze(mean(i5.depth(1:2,:,:))),[1 2]);
plotaou = repmat(squeeze(i5.aou(2,:,:)-i5.aou(1,:,:)),[1,2]);
% contourf(plotlon,plotdepth,plotaou,(-1:0.05:1)*0.02)
pcolor(plotlon,plotdepth,plotaou)
set(gca,'ydir','reverse')
ylim([0 1200]);xlim([30 120]); caxis([-1 1]*0.01)
colorbar
; colorbar;
title('1995-1987')
shading flat;

subplot(3,1,3); hold on;
hold on
plotlon= [i5.lon i5.lon+360];
plotdepth = repmat(squeeze(mean(i5.depth([1 3],:,:))),[1 2]);
plotaou = repmat(squeeze(i5.aou(3,:,:)-i5.aou(1,:,:)),[1,2])./1025*1e6;
% plotaou(19,:,:)=NaN;
% contourf(plotlon,plotdepth,plotaou,(-1:0.05:1)*0.02)
contourf(plotlon,plotdepth,plotaou,-20:1:20,'LineStyle','none')
% shading flat;
% contour(plotlon,plotdepth,plotaou,[-2 -2],'LineColor','black')
set(gca,'ydir','reverse')
ylim([0 1200]);xlim([30 120]); caxis([-10 10])
cax= colorbar;
ylabel(cax,'AOU micromol/kg')
title('2002-1987')
%%
clf
% subplot(3,1,3); hold on;
hold on
plotlon= [i5.lon i5.lon+360];
plotdepth = repmat(squeeze(i5.depth(3,:,:)),[1 2]);
plotaou = repmat(squeeze(i5.aou(3,:,:)),[1,2])./1025*1e6;
plotoxygen = repmat(squeeze(i5.oxygen(3,:,:)),[1 2])*1e3;
pcolor(plotlon,plotdepth,plotaou);
% plotaou(19,:,:)=NaN;
% contourf(plotlon,plotdepth,plotaou,(-1:0.05:1)*0.02)
% contourf(plotlon,plotdepth,plotoxygen,180:5:250,'LineStyle','none')

% shading flat;
% contour(plotlon,plotdepth,plotaou,[-2 -2],'LineColor','black')
set(gca,'ydir','reverse')
ylim([0 2000]);xlim([30 120]);
cax= colorbar;
ylabel(cax,'AOU micromol/kg')
title('2002-1987')
%% HovMoeller SAMW
start4d = [5 18 0 0];
count4d = [-1 1 -1 -1];
stride4d = [12 1 1 1];
samw.oxygen= nc_varget(infile.oxygen,'mn_oxygen',start4d,count4d,stride4d);
samw.o2sat= nc_varget(infile.oxygen,'mn_o2sat',start4d,count4d,stride4d);
samw.h = nc_varget(infile.hfile,'h',start4d,count4d,stride4d);
samw.salt = nc_varget(infile.saltfile,'salt',start4d,count4d,stride4d);
%%

clf; subplot(1,3,3);
plotlon= [i5.lon i5.lon+360];
plotaou = samw.o2sat-samw.oxygen;
plotaou = plotaou-repmat(mean(plotaou),[60 1]);
plotaou = repmat(squeeze(plotaou(:,i5.latidx,:)),[1 2]);
plotaou = plotaou./1025*1e6;
ploth = repmat(squeeze(samw.h(1:60,i5.latidx,:)),[1 2]);
plotsalt = repmat(squeeze(samw.salt(1:60,i5.latidx,:)),[1 2]);

contourf(plotlon(1,:),1947:2006,plotaou,[-1:0.05:1]*30,'LineStyle','none');
grid on;
cax = colorbar;
xlabel('Longitude')
ylabel(cax,'AOU-mean(AOU) [mol kg^{-1}]')
xlim([30 120]); caxis([-10 10])
daspect([20 10 1])

subplot(1,3,1)
ploth=ploth-repmat(mean(ploth),[60,1]);
contourf(plotlon(1,:),1947:2006,ploth,-120:10:120,'Linestyle','none')
grid on;
cax=colorbar;
caxis([-100 100])
xlim([30 120]);
xlabel('Longitude')
ylabel(cax,'h-mean(h) [m]')
daspect([20 10 1])

subplot(1,3,2)
plotsalt=plotsalt-repmat(mean(plotsalt),[60,1]);
contourf(plotlon(1,:),1947:2006,plotsalt,-0.3:0.025:0.3,'Linestyle','none')
grid on;
caxis([-0.150001 0.15])
cax=colorbar;
xlim([30 120]);

xlabel('Longitude')
ylabel(cax,'S-mean(S) [PSU]')
daspect([20 10 1])
saveas(gcf,'/ltraid3/ashao/uw-apl/figs/o2_hindcast/i5.aou.h.salt.eps','epsc')
% print /ltraid3/ashao/uw-apl/figs/o2_hindcast/i5.aou.h.salt.png -r200 -dpng
%% AOU Plot averaged on isopycnals
i5.isopycnalavg.aou = zeros(3,49);
i5.isopycnalavg.depth = zeros(3,49);
i5.isopycnalavg.h = zeros(3,49);
lon = [i5.lon i5.lon+360];
lonidx = lon(1,:)>60 & lon(1,:)<80;
layers = nc_varget('/ltraid4/ashao/HIM/himw/him_sis/INPUT/Global_HIM_IC.nc','Layer')';
for t = 1:3
    aou = repmat(squeeze(i5.aou(t,:,:)),[1,2]);    
    h = repmat(squeeze(i5.h(t,:,:)),[1,2]);
%     depth = repmat(squeeze(i5.depth(t,:,:)),[1,2]);
    for lay = 1:49
        
        wts = h(lay,lonidx)./nansum(h(lay,lonidx));
        i5.isopycnalavg.aou(t,lay) = nansum(aou(lay,lonidx).*wts).*layers(lay);
        i5.isopycnalavg.h(t,lay) = nansum(h(lay,lonidx).*wts);
        
    end
    
    
end
i5.isopycnalavg.depth = mean(cumsum(i5.isopycnalavg.h,2));

clf; 
subplot(1,2,1); hold on;
plot(i5.isopycnalavg.aou(1,:),i5.isopycnalavg.depth,'r-x','LineWidth',2)
plot(i5.isopycnalavg.aou(2,:),i5.isopycnalavg.depth,'b-x','LineWidth',2)
plot(i5.isopycnalavg.aou(3,:),i5.isopycnalavg.depth,'k-x','LineWidth',2)
legend('1987','1995','2002')
ylim(i5.isopycnalavg.depth([10 30]))
set(gca,'ydir','reverse')
grid on;

subplot(1,2,2); hold on;
plot(i5.isopycnalavg.aou(3,:)-i5.isopycnalavg.aou(2,:),i5.isopycnalavg.depth,'b-x','LineWidth',2)
plot(i5.isopycnalavg.aou(3,:)-i5.isopycnalavg.aou(1,:),i5.isopycnalavg.depth,'r-x','LineWidth',2)
plot(i5.isopycnalavg.aou(2,:)-i5.isopycnalavg.aou(1,:),i5.isopycnalavg.depth,'g-x','LineWidth',2)

scatter(i5.isopycnalavg.aou(3,19:20)-i5.isopycnalavg.aou(2,:),i5.isopycnalavg.depth,'b-x','LineWidth',2)
scatter(i5.isopycnalavg.aou(3,19:20)-i5.isopycnalavg.aou(1,:),i5.isopycnalavg.depth,'r-x','LineWidth',2)
scatter(i5.isopycnalavg.aou(2,19:20)-i5.isopycnalavg.aou(1,:),i5.isopycnalavg.depth,'g-x','LineWidth',2)


line([0 0],[0 5000],'Color','black')
grid on;
ylim(i5.isopycnalavg.depth([10 30]))
legend('2002-1995','2002-1987','1995-1987')
set(gca,'ydir','reverse')