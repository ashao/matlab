iodtab = urlread('http://www.jamstec.go.jp/frcgc/research/d1/iod/DATA/dmi_HadISST.txt');
iodtab = str2num(iodtab);
iod.idx = iodtab(:,3);
iod.time = datenum(iodtab(:,1),iodtab(:,2),15);


%% Hindcasts
wdfile = '/ltraid4/ashao/HIM/hyak_store/HINDCAST/WD-hind.nc';
himhind.wd = nc_varget(wdfile,'wd',[0 1 0 0],[inf 1 inf inf]);

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
indian.time = himhind.time(himtruncidx);
iod.time = iod.time(iodtruncidx);
iod.idx = iod.idx(iodtruncidx);
ntime = length(indian.time);
% truncate spatially

%% Calculate monthly volume flux
indian.subduction = repmat(metrics.Ah.data(spatialidx)',[ntime 1]);
indian.subduction = indian.subduction.*indian.wd;

%% Annual sum of volume flux

[ntime npts] = size(indian.subduction);
indian.annualsubduction = zeros(ntime/12,npts);
timelen = repmat(eomday(2014,1:12)',[1 npts])*86400;

sidx = 1;
eidx = 12;
counter = 0;
while eidx <= ntime;
    counter = counter + 1;
    idx = sidx:eidx;
    iod.annual(counter) = mean(iod.idx(idx));
    indian.annualsubduction(counter,:) = sum(indian.subduction(idx,:).*timelen);
    sidx = sidx+12;
    eidx = sidx + 11;
    
end

%% Plot annual subduction 
plotmat = zeros(size(metrics.geolat.data));
for t = 1:ntime/12
    clf
    
    worldmap(latrange,lonrange)
    plotmat(spatialidx) = indian.annualsubduction(t,:);
    pcolorm(metrics.geolat.data,metrics.geolon.data,plotmat)
    title(sprintf('Year %d',1948+t-1));    
    geoshow('landareas.shp')
    colorbar
    caxis([-1 1]*1e12)
    colormap(othercolor('RdYlBu5'))
    drawnow
    
end
%% Plot variability
clf
    
    worldmap(latrange,lonrange)
    plotmat(spatialidx) = std(indian.annualsubduction)./mean(std(indian.annualsubduction));
    pcolorm(metrics.geolat.data,metrics.geolon.data,plotmat)
    title(sprintf('Year %d',1948+t-1));    
    geoshow('landareas.shp')
    colorbar
    caxis([-1 1]*1e12)
    colormap(othercolor('RdYlBu5'))
    drawnow
%% Regress onto IUOD and plot R^2
for ipt = 1:npts
   
    [C lags] = xcov(indian.annualsubduction(:,ipt),iod.annual,5,'coeff');
    [null maxidx] = max(abs(C));
    indian.R(ipt) = C(maxidx);
    indian.lag(ipt) = lags(maxidx);
    
end
%% Plot regressions and lags
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
title('Correlation')
caxis([-5 5])
geoshow('landareas.shp','FaceColor',[0.5 0.5 0.5])
colorbar
