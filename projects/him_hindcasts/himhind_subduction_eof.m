wdfile = '/ltraid4/ashao/HIM/hyak_store/HINDCAST/WD-hind.nc';
wd = nc_varget(wdfile,'wd',[0 1 0 0],[inf 1 inf inf]);
himhind.time = nc_varget(wdfile,'Time');
himhind.time = (himhind.time - himhind.time(1))/365+1948;
load metrics
[nmonths nlat nlon] = size(wd);
nyears = 60;
load('~/uw-apl/data/ncep/sam_from_ncep.mat');
%%
padidx = 270;
himhind.wd = zeros(720,210,360+padidx);
for mon = 1:nmonths
   
    himhind.wd(mon,:,:) = [squeeze(wd(mon,:,:)) squeeze(wd(mon,:,1:padidx))];
    
end
%%
geolat = [metrics.geolat.data metrics.geolat.data(:,1:padidx)];
geolon = [metrics.geolon.data wrapTo360(metrics.geolon.data(:,1:padidx))];
wet = logical([metrics.wet.data metrics.wet.data(:,1:padidx)]);

indian.latrange = [-65 -20];
indian.lonrange = [30 160];
indian.mask = geolat > min(indian.latrange) & geolat < max(indian.latrange) & ...
    geolon < max(indian.lonrange) & geolon > min(indian.lonrange);
%%
himhind.annual.wd = zeros([nyears size(geolat)]);

for year = 1:nyears
    sidx = (year-1)*12+1;
    eidx = sidx+11;
    himhind.annual.wd(year,:,:) = mean(himhind.wd(sidx:eidx,:,:));
    
end
%%
indian.geolat = geolat(indian.mask);
indian.geolon = geolon(indian.mask);
indian.wet = wet(indian.mask);
indian.wd = himhind.annual.wd(:,indian.mask & wet);
[nyears npts] = size(indian.wd);

% data = indian.wd(:,indian.wet);
data = data - repmat(mean(data),[nyears 1 1]);
data = detrend(data')';
[U S V] = svd(data);
eofs = S*V;
%%
lambda = diag(S);
mode = 1;
clf
eof = nan(size(geolat));
eof(indian.mask & wet) = sum(eofs(1:10,:));
subplot(3,2,[1 2 3 4])
worldmap(indian.latrange,indian.lonrange)
% contourfm(geolat,geolon,eof*1e6,-2:0.1:2,'LineStyle','none')
contourfm(geolat,geolon,-eof*1e5,-2:0.2:2,'LineStyle','none')
caxis([-1 1])
cax = colorbar;
% ylabel(cax,'Ekman pumping [10^-6 1/(m s)]')
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5])
pcs = U*S;

%%
subplot(3,2,[5 6]); cla; hold on

pc = sum(pcs(:,1:2),2);
area(datenum(1948:2007,6,15),-pc./std(pc))
plot(sam_ncep.time,smooth(sam_ncep.pc./std(sam_ncep.pc),12,'mean'),'k-','LineWidth',2)

ylim([-2 2]);
xlim(datenum([1948 2008],1,1))
datetick('KeepTicks','KeepLimits');
legend('PC','SAM')
%% Regress onto SAM
idx = sam_ncep.time > datenum(1948,1,1) & sam_ncep.time < datenum(2008,1,1);
temptime = sam_ncep.time(idx);
tempsam = sam_ncep.pc(idx);
for t=1:nyears
    sidx = (t-1)*12 + 1;
    eidx = sidx+11;
    annualsam.time(t) = mean(temptime(sidx:eidx));
    annualsam.idx(t) = mean(tempsam(sidx:eidx));
    
end

for pt = 1:npts
    
    [wdcorr(pt) wdlag(pt)] = maxcorr(annualsam.idx,indian.wd(:,pt),20,'coeff');
    
end
%%
clf
subplot(2,1,1)
worldmap(indian.latrange,indian.lonrange)
mapdata = nan(size(geolat));
mapdata(indian.mask & wet) = wdcorr;
% contourfm(geolat,geolon,mapdata,-1:0.1:1,'LineStyle','none')
pcolorm(geolat,geolon,mapdata)
colorbar
subplot(2,1,2)
worldmap(indian.latrange,indian.lonrange)
mapdata = nan(size(geolat));
mapdata(indian.mask & wet) = wdlag;
pcolorm(geolat,geolon,mapdata)
colorbar
% sam_ncep.annual.time = 