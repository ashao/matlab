iodtab = urlread('http://www.jamstec.go.jp/frcgc/research/d1/iod/DATA/dmi_HadISST.txt');
iodtab = str2num(iodtab);
iod.idx = iodtab(:,3);
iod.time = datenum(iodtab(:,1),iodtab(:,2),15);

%% Hindcasts
wdfile = '/ltraid4/ashao/HIM/hyak_store/HINDCAST/WD-hind.nc';
rmlfile = '/ltraid4/ashao/HIM/hyak_store/HINDCAST/Rml-hind.nc';
himhind.wd = nc_varget(wdfile,'wd',[0 2 0 0],[inf 1 inf inf]);
himhind.rml = nc_varget(rmlfile,'Rml',[0 0 0 0],[inf 1 inf inf]);

himhind.time = nc_varget(wdfile,'Time');
himhind.time = (himhind.time - himhind.time(1))/365+1948;
%% truncate time index
iodtruncidx = iod.time < datenum(2008,1,1);
himtruncidx = himhind.time >= 1958;
load metrics
latrange = [-60 25];
lonrange = [20 -220];
latidx = metrics.geolat.data > min(latrange) & metrics.geolat.data < max(latrange);
lonidx = metrics.geolon.data > max(lonrange) | metrics.geolon.data < min(lonrange);
spatialidx = latidx & lonidx;

indian.wd = himhind.wd(himtruncidx,spatialidx);
indian.rml = himhind.rml(himtruncidx,spatialidx);
indian.time = himhind.time(himtruncidx);
iod.time = iod.time(iodtruncidx);
iod.idx = iod.idx(iodtruncidx);
%% Calculate mass flux at base of mixed layer
[ntime npts] = size(indian.subduction);
indian.mlflux = zeros(ntime/12,npts);
timelen = repmat(eomday(2014,1:12)',[1 npts])*86400;

sidx = 1;
eidx = 12;
counter = 0;
while eidx <= ntime;
    counter = counter + 1;
    idx = sidx:eidx;
    iod.annual(counter) = mean(iod.idx(idx));
    indian.mlflux(counter,:) = sum(indian.wd(idx,:).*timelen.*indian.rml(idx,:));
    sidx = sidx+12;
    eidx = sidx + 11;
    
end
%% Plot annual mass flux 
plotmat = zeros(size(metrics.geolat.data));
for t = 1:ntime/12
    clf
    
    worldmap(latrange,lonrange)
    plotmat(spatialidx) = indian.mlflux(t,:);
    pcolorm(metrics.geolat.data,metrics.geolon.data,plotmat)
    title(sprintf('Year %d',1948+t-1));    
    geoshow('landareas.shp')
    colorbar
    caxis([-1 1]*1e5)
    colormap(othercolor('RdYlBu5'))
    drawnow
    pause
    
end

%% Regress onto IUOD and plot R^2
for ipt = 1:npts
   
    [C lags] = xcov(indian.mlflux(:,ipt),iod.annual,5,'coeff');
    [null maxidx] = max(abs(C));
    indian.R(ipt) = C(maxidx);
    indian.lag(ipt) = lags(maxidx);
   
end
plotmat = zeros(size(metrics.geolat.data));
clf
subplot(2,1,1)
worldmap(latrange,lonrange)
plotmat(spatialidx) = indian.R.^2;
pcolorm(metrics.geolat.data,metrics.geolon.data,plotmat)
title('Correlation')
caxis([0 0.3])
geoshow('landareas.shp','FaceColor',[0.5 0.5 0.5])
colorbar

subplot(2,1,2)
worldmap(latrange,lonrange)
plotmat(spatialidx) = indian.lag;
pcolorm(metrics.geolat.data,metrics.geolon.data,plotmat)
title('Lags')
caxis([-5 5])
geoshow('landareas.shp','FaceColor',[0.5 0.5 0.5])
colorbar
