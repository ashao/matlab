inpath = 'C:\Users\ashao\Data\altimetry\timevary\annual\processed\';
files = dir([inpath 't*.mat']);
nfiles = length(files);
trackfile = 'C:\Users\ashao\Data\vxxc_matlab\groundpath.mat';
load(trackfile)
nwindows = 19;
%%
timevary.north.lat= zeros(254,nwindows);
timevary.mean.lat = zeros(254,nwindows);
timevary.south.lat = zeros(254,nwindows);
timevary.north.lon= zeros(254,nwindows);
timevary.mean.lon = zeros(254,nwindows);
timevary.south.lon = zeros(254,nwindows);
numnull = 0;
nullidx = [];
for tidx = 1:nfiles
    
    load([inpath files(tidx).name]);
    for winidx = 1:nwindows
        
        
        optpar = cell2mat(opt_track(winidx).optpar');
        if ~isempty(optpar)
            lats = optpar(:,2);
            
            timevary.north.lat(tidx,winidx) = max(lats);
            timevary.north.lon(tidx,winidx) = interp1( ...
                groundpath(tidx).lat,groundpath(tidx).lon,timevary.north.lat(tidx,winidx));
            
            timevary.mean.lat(tidx,winidx) = nanmean(lats);
            timevary.mean.lon(tidx,winidx) = interp1( ...
                groundpath(tidx).lat,groundpath(tidx).lon,timevary.mean.lat(tidx,winidx));
            
            timevary.median.lat(tidx,winidx) = nanmedian(lats);
            timevary.median.lon(tidx,winidx) = interp1( ...
                groundpath(tidx).lat,groundpath(tidx).lon,timevary.median.lat(tidx,winidx));
            
            timevary.south.lat(tidx,winidx) = min(lats);
            timevary.south.lon(tidx,winidx) = interp1( ...
                groundpath(tidx).lat,groundpath(tidx).lon,timevary.south.lat(tidx,winidx));
            
        else
            numnull = numnull+1;
            nullidx(numnull) = tidx;
        end
    end
end
%%
nulls = unique(nullidx);
nnulls = length(nulls);
extents = {'north','mean','median','south'};
for i = 1:length(extents)
    
    for t = 1:nnulls

        timevary.(extents{i}).lat(nulls(t),:) = [];
        timevary.(extents{i}).lon(nulls(t),:) = [];
    end
end

%%
for winidx = 1:nwindows
    timevary.time(winidx) = opt_track(winidx).time;
end
%% Detrend the data

timevary.north.detrended = detrend(timevary.north.lat')';
timevary.mean.detrended = detrend(timevary.mean.lat')';
timevary.median.detrended = detrend(timevary.median.lat')';
timevary.south.detrended = detrend(timevary.south.lat')';

%% Calculate the trend
for tidx = 1:(nfiles-nnulls)
    
    [coeffs s] = polyfit(timevary.time,timevary.north.lat(tidx,:),1);
    Rinv = inv(s.R); c = (Rinv*Rinv')*s.normr^2/s.df;
    se = sqrt(diag(c))';
    t = coeffs ./ se;
    tnorth(tidx) = t(1);
    timevary.north.trend(tidx) = coeffs(1);
    
    [coeffs s] = polyfit(timevary.time,timevary.mean.lat(tidx,:),1);
    Rinv = inv(s.R); c = (Rinv*Rinv')*s.normr^2/s.df;
    se = sqrt(diag(c))';
    t = coeffs ./ se;
    tmean(tidx) =t(1);
    timevary.mean.trend(tidx) = coeffs(1);
    
    %     [coeffs s]= polyfit(timevary.time,timevary.median.lat(tidx,:),1);
    timevary.median.trend(tidx) = coeffs(1);
    
    [coeffs s]= polyfit(timevary.time,timevary.south.lat(tidx,:),1);
    Rinv = inv(s.R); c = (Rinv*Rinv')*s.normr^2/s.df;
    se = sqrt(diag(c))';
    t = coeffs ./ se;
    tsouth(tidx) =t(1);
    timevary.south.trend(tidx) = coeffs(1);
end

%%
clf
extents = {'north','mean','median','south'};
titles = {'Northern','Mean','Median','Southern'};
for i = 1:length(extents)
    subplot(5,1,i)
    
    plotlon = nanmean(timevary.(extents{i}).lon,2);
    [plotlon sortidx] = sort(plotlon);
    area(plotlon,smooth(timevary.(extents{i}).trend(sortidx),5,'mean')*110*365);
    set(gca,'xtick',0:60:360)
    grid on
    title([titles{i} ' Extent'])
    dar = daspect;
    ylim([-15 15])
    ylabel('Trend (km / yr)')
    timevary.(extents{i}).plotlon = plotlon;
    timevary.(extents{i}).plotidx = sortidx;
end
subplot(5,1,5)
axesm('MapProjection','Mercator','MapLatLimit',[-70 -20],'MapLonLimit',[0 360],'Frame','On','Grid','On');
geoshow('landareas.shp', 'FaceColor', [1 1 1]-0.1)
tightmap
daspect([2 4.5 1])
plabel;
mlabel([-70])
%% Correlation with ECMWF wind stress curl
load C:\Users\ashao\Data\ecmwf.annual.taucurl.mat
%%
latrange = [-60 -40];
lonrange = [0 360];
clf

timerange = [datenum(1993,1,1) datenum(2012,1,1)];
latidx = ecmwf.latitude >= min(latrange) & ecmwf.latitude <= max(latrange);
lonidx = ecmwf.longitude >= min(lonrange) & ecmwf.longitude <= max(lonrange);
tidx = ecmwf.time >= min(timerange) & ecmwf.time < max(timerange);

acc.taucurl_ann = ecmwf.taucurl_ann(tidx,latidx,lonidx);
acc.time = ecmwf.time(tidx);
acc.latitude = ecmwf.latitude(latidx);
acc.longitude = ecmwf.longitude(lonidx);
[longrid latgrid] = meshgrid(acc.longitude,acc.latitude);


notnan = ~isnan(acc.taucurl_ann(1,:,:));
wts = [];
wts(1,:,:) = cosd(latgrid);
wts = wts(:,notnan);
wts = repmat(wts,[length(acc.time),1]);

norm = sum(wts(1,:));
wts = wts./norm;

data = detrend(acc.taucurl_ann(1:end,notnan));

[U S V] = svds(data);

acc.pcs = (U*S);
acc.pervar = diag(S).^2/sum(diag(S).^2);

modepattern = zeros(size(latgrid));
modepattern(notnan) = V(:,1);
worldmap([-90 -30],[0 360])
pcolorm(latgrid,longrid,modepattern)
geoshow('landareas.shp', 'FaceColor', [1 1 1]-0.1)
curlidx = acc.pcs(:,1);
[C lags]= xcorr(detrend(curlidx)./std(detrend(curlidx)),detrend(sam.idx)./std(detrend(sam.idx)),'coeff');
[maxc I] = max(C)
lags(I)
%%
extents = {'north','mean','median','south'};
nmodes = 3;
nlags = 2;
fprintf('\n')
for i = 1:length(extents)
    lonidx = timevary.(extents{i}).plotlon > min(lonrange) & ...
        timevary.(extents{i}).plotlon < max(lonrange);
    [U S V] = svd(timevary.(extents{i}).detrended(lonidx,:)');
    timevary.(extents{i}).modevar = diag(S).^2/sum(diag(S).^2);
    timevary.(extents{i}).pcs = U*S;
%     subplot(5,1,i+1); hold on;
%     plot(cumsum(timevary.(extents{i}).modevar))
%     plot((timevary.(extents{i}).modevar))
%     
%     grid on;
%     ylim([0 1])
    if nmodes > 1
        varidx = sum(timevary.(extents{i}).pcs(:,1:nmodes),2);
    else
        varidx = timevary.(extents{i}).pcs(:,1);
    end
    
    [C lags] = xcorr(varidx,curlidx,nlags,'coeff');    
    [timevary.(extents{i}).curlcorr maxidx] = max(abs(C));    
    timevary.(extents{i}).curllag = lags(maxidx);
    
%     C = corrcoef(varidx,curlidx);
%     timevary.(extents{i}).curlcorr = C(1,2);
    
    fprintf('%s Curl \tCorrelation: %f Lag: %d\n',extents{i},timevary.(extents{i}).curlcorr,timevary.(extents{i}).curllag);
end
