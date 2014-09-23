h = nc_varget('/ltraid3/ashao/uw-apl/data/offtrac/input/normalyear/H-clim.nc','h',[0 0 0 0],[12 1 -1 -1])*2;
load metrics
vol = zeros(12,210,360);

adjidx = 12;

for idx = 1:12
   
    vol(idx,:,:) = h(adjidx,:,:);
    adjidx = adjidx - 1;
end
iswet = logical(metrics.wet.data);
vol(:,~logical(metrics.wet.data))=NaN;
%%
adjointttd = nc_varget('/home/ashao/uw-apl/development/Offtrac/branches/adjoint/adjoint.normalyear.12010.nc','mn_adjoint');
[ntime nlat nlon] = size(adjointttd);
%%
% calcarray = adjointttd;
monidx = 1;

for t= 1:ntime
%     fprintf('Month %d/%d\n',t,ntime)
    adjointttd(t,:,:) = adjointttd(t,:,:).*vol(monidx,:,:);
    monidx = monidx+1;
    if monidx==12;
        monidx=1;
    end
end

truncttd = adjointttd(:,iswet);
%% Calculate normalization factors
t=(1:ntime)/12;
timenorm = trapz(t,truncttd);
globalnorm = sum(timenorm);
%% Calculate meanage
meanage = nan(210,360);
meanage(iswet) = trapz(t,bsxfun(@times,truncttd,t'))./timenorm;

%%
coast = load('coast.mat');
clear geolat
clear geolon
data = nan(210,360);
data(iswet) = timenorm/globalnorm;
% geolat= [metrics.geolat.data metrics.geolat.data];
% geolon= [metrics.geolon.data metrics.geolon.data+360];
% data = [data data];
geolat = metrics.geolat.data;
geolon = metrics.geolon.data;
geolat = [geolat geolat(:,1)];
geolon = [geolon geolon(:,1)+360];
data = [data data(:,1)];
%%
clf
axesm('gstereo','MapLatLimit',[-80 80],'MapLonLimit',[-60 300])

contourfm(geolat,geolon,log10(data),-13:-1);
% [c v] = contourfm(geolat,geolon,log10(data),-12:1:0,'LineColor','k','LineWidth',1);
% colormap(othercolor('Greens9'))
colormap(flipud(lbmap(20,'brownblue')))

% scatterm(metrics.geolat.data(56,29),metrics.geolon.data(56,29),'kx')
% clabelm(c,v,[-3 -6 -9])

% [c v]= contourm(geolat,geolon,log10(data),-10:1:0,'LineColor','k','LineWidth',1);
scatterm(metrics.geolat.data(56,29),metrics.geolon.data(56,29),200,'kx','LineWidth',4)
scatterm(metrics.geolat.data(56,29),metrics.geolon.data(56,29),100,'cx','LineWidth',2)
plotm(coast.lat,coast.long,'k')
% clabelm(c, v,-8:2:0,'LabelSpacing',1e9)
cax = contourcbar;
caxis([-13 1])
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]+0.2);
ylabel(cax,'log10(Volume Fraction)')
 tightmap;plabel;mlabel;
% print ~/uw-apl/figs/adjoint/adjoint_volume_fraction.png -dpng -r200

%% Mean Age

% meanage(~logical(metrics.wet.data)) = NaN;
geolat= [metrics.geolat.data metrics.geolat.data(:,1)];
geolon= [metrics.geolon.data metrics.geolon.data(:,1)+360];
data = [meanage meanage(:,1)];

clf;cla
ax = axesm('gstereo','MapLatLimit',[-80 80],'MapLonLimit',[-60 300]);
levels = [1:10 20:10:100 150:50:1000];
levels = [10:10:100 150:50:1000];
levels = [10 25 50:50:350 450 550 700];
contourfm(geolat,geolon,data,levels,'LineColor','k')

% [c h] = contourm(geolat,geolon,data,[0:20:100 200:200:1000],'LineColor','black');
% clabelm(c,h,'LabelSpacing',1e6)
scatterm(metrics.geolat.data(56,29),metrics.geolon.data(56,29),200,'rx','LineWidth',2)
% cax = cbarf([1 1000],[1:10 20:10:100 150:50:1000])
% cax = colorbar;
colormap(flipud(othercolor('RdYlBu10')))
cax = contourcbar('peer',ax);
set(get(cax,'YLabel'),'String','Adjoint Mean Age [yr)]');
% ylabel(cax,'Adjoint Mean Age [yr)]')
colormap(flipud(lbmap(24,'redblue')))
set(cax,'Ytick',levels)
% caxis([0 2.5])
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
tightmap;plabel;mlabel;

% rgb2cm
% print ~/uw-apl/figs/adjoint/adjoint_mean_age.png -dpng -r200
%% Mean Age Contribution
tprime = [];
tprime(1,1,1:ntime) = (1:ntime)/12;
tprime = permute(repmat(tprime,[210 360,1]),[3 1 2]);
norm2d = squeeze(sum(adjointttd));
clear wt
wt(1,:,:) = norm2d;
wt = repmat(wt,[ntime,1,1]);
%%
% meanage = squeeze(sum(tprime.*adjointttd./wt)).*normttd;
% meanage(~logical(metrics.wet.data)) = NaN;

geolat= [metrics.geolat.data metrics.geolat.data(:,1)];
geolon= [metrics.geolon.data metrics.geolon.data(:,1)+360];
data = nan(210,360);
data(iswet) = meanage(iswet).*timenorm'./globalnorm;
data = [data data(:,1)];
clf;cla
axesm('gstereo','MapLatLimit',[-80 80],'MapLonLimit',[-60 300])
contourfm(geolat,geolon,log10(data),-10:0.25:-1,'LineStyle','none')
[c h] = contourm(geolat,geolon,log10(data),-10:2:0,'LineColor','black');
clabelm(c,h)
scatterm(metrics.geolat.data(56,29),metrics.geolon.data(56,29),200,'rx','LineWidth',2)
cax = colorbar;
ylabel(cax,'Age contribution ( log10(yr) )')
% caxis([0 1000])
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
caxis([-10 -1])
tightmap;plabel;mlabel;
saveas(gcf,'/ltraid3/ashao/uw-apl/figs/adjoint/mean_age_contribution.eps','epsc')
% saveas([