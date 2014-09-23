inpath = '/ltraid3/ashao/uw-apl/projects/saf_altimetry/timevary/annual/processed/curvefit/';
files = dir([inpath 't*.mat']);
nfiles = length(files);
trackfile = '/ltraid4/aviso/alongtrack/sla/vxxc_matlab/groundpath.mat';
load(trackfile)
nwindows = 20;
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
    
    nwindows = length(opt_track);
    for winidx = 1:nwindows
        
        
        optpar = cell2mat(opt_track(winidx).optpar');
        if ~isempty(optpar) & nwindows == 20
%             ndata = cellfun(@length,opt_track(winidx).skewness);
            R2 = cell2mat(opt_track(winidx).R2);
            R2idx = R2 < 0.3;
%             any(R2idx)
            if any(R2idx)
                
                lats = optpar(R2idx,2);
                
                timevary.north.lat(tidx,winidx) = max(lats);
                timevary.north.lon(tidx,winidx) = interp1( ...
                    groundpath(tidx).lat,groundpath(tidx).lon,timevary.north.lat(tidx,winidx),'linear','extrap');
                
                timevary.mean.lat(tidx,winidx) = nanmean(lats);
                timevary.mean.lon(tidx,winidx) = interp1( ...
                    groundpath(tidx).lat,groundpath(tidx).lon,timevary.mean.lat(tidx,winidx),'linear','extrap');
                
                timevary.median.lat(tidx,winidx) = nanmedian(lats);
                timevary.median.lon(tidx,winidx) = interp1( ...
                    groundpath(tidx).lat,groundpath(tidx).lon,timevary.median.lat(tidx,winidx),'linear','extrap');
                
                timevary.south.lat(tidx,winidx) = min(lats);
                timevary.south.lon(tidx,winidx) = interp1( ...
                    groundpath(tidx).lat,groundpath(tidx).lon,timevary.south.lat(tidx,winidx),'linear','extrap');
            else
                timevary.north.lat(tidx,winidx) = NaN;
                timevary.south.lat(tidx,winidx) = NaN;
                timevary.mean.lat(tidx,winidx) = NaN;
                timevary.median.lat(tidx,winidx) = NaN;
                timevary.north.lon(tidx,winidx) = NaN;
                timevary.south.lon(tidx,winidx) = NaN;
                timevary.mean.lon(tidx,winidx) = NaN;
                timevary.median.lon(tidx,winidx) = NaN;
            end
        else
            numnull = numnull+1;
            nullidx(numnull) = tidx;
        end
    end
    
end
%%
nulls = unique(nullidx);
% nulls = 64:254;
nnulls = length(nulls)
extents = {'north','mean','median','south'};


for i = 1:length(extents)
        

        timevary.(extents{i}).lat(nulls,:) = [];
        timevary.(extents{i}).lon(nulls,:) = [];
    
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
extents = {'north','mean','median','south'};
for tidx = 1:(nfiles-nnulls)
    
    for i = 1:length(extents)
        extent = extents{i};
    
    mdl = LinearModel.fit(timevary.time,timevary.(extent).lat(tidx,:));
    timevary.(extent).trend(tidx)=mdl.Coefficients.Estimate(2);
    timevary.(extent).trend_pvalue(tidx)=mdl.Coefficients.pValue(2);
    
    end
end
%%
clf
hold on
extents = {'north','mean','south'};
titles = {'Northern','Mean','Southern'};
colors = {'r','k','b'};
for i = 1:length(extents)
    notnan = ~isnan(timevary.(extents{i}).lat(:,1)) & ~isnan(timevary.(extents{i}).lon(:,1));
    [meanlat plotlon] = meanm(timevary.(extents{i}).lat',timevary.(extents{i}).lon');
    
    
    [plotlon sortidx] = sort(wrapTo360(plotlon));
    
        timevary.(extents{i}).plotlon = plotlon;
%     timevary.(extents{i}).sortidx = sortidx;
    
%     timevary.(extents{i}).std = std(timevary.(extents{i}).lat,0,2);
%     timevary.(extents{i}).std = timevary.(extents{i}).std(timevary.(extents{i}).sortidx);

    plottrend = timevary.(extents{i}).trend;
    trend_pvalue = timevary.(extents{i}).trend_pvalue;;
    sigidx = trend_pvalue<0.05 & notnan';
    
    stem(plotlon(sigidx),plottrend(sigidx)*365*110,colors{i},'filled')    
    
%     set(gca,'xtick',0:60:360)
    grid on
    title([titles{i} ' Extent'])
    dar = daspect;
    ylim([-1 1]*30)
    ylabel('Trend (km / yr)')
    xlim([0 360])
%     
%     subplot(2,1,2); hold on
%     lats = timevary.(extents{i}).lat(sortidx,:);
%     xlim([ 1 20])
%     plot(lats(sigidx,:),colors{i})
end
% xlim([0 360])
% subplot(2,1,2)
% axesm('MapProjection','Mercator','MapLatLimit',[-70 -20],'MapLonLimit',[0 360],'Frame','On','Grid','On');
% geoshow('landareas.shp', 'FaceColor', [1 1 1]-0.1)
% tightmap
% daspect([2 4.5 1])
% plabel;
% mlabel([-70])
%% Correlation with ECMWF wind stress curl
load  /ltraid3/ashao/uw-apl/projects/saf_altimetry/ecmwf.annual.taucurl.mat
load  ~/uw-apl/projects/saf_altimetry/ecmwf.sam.mat
%%
latrange = [-70 -30];
lonrange = [30 120];
figure(1);clf

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
curlidx = acc.pcs(:,1)./std(acc.pcs(:,1));
samidx = interp1(sam.annual.time,sam.annual.idx,acc.time);
[C lags]= xcov(curlidx,samidx,'coeff');
[maxc I] = max(abs(C))
lags(I)

extents = {'north','mean','median','south'};
for i = 1:length(extents)
    scatterm(nanmean(timevary.(extents{i}).lat,2),nanmean(timevary.(extents{i}).lon,2))
end
    

figure(2);clf
for i=1:length(extents)
    for nmodes =  1:20 ;
        nlags = 3;
        fprintf('\nCorrelation to SAM and Curl %d Modes\n',nmodes)
        
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
        
        [C lags] = xcov(curlidx,varidx,nlags,'coeff');
        [timevary.(extents{i}).curlcorr(nmodes) maxidx] = max(abs(C));
        timevary.(extents{i}).curllag(nmodes) = lags(maxidx);
        [C lags] = xcov(samidx,varidx,nlags,'coeff');
        [timevary.(extents{i}).samcorr(nmodes) maxidx] = max(abs(C));
        timevary.(extents{i}).samlag(nmodes) = lags(maxidx);
        %     C = corrcoef(varidx,curlidx);
        %     timevary.(extents{i}).curlcorr = C(1,2);
        
        fprintf('%s Curl \tCorrelation: %f Lag: %d\n',extents{i},timevary.(extents{i}).curlcorr(nmodes),timevary.(extents{i}).curllag(nmodes));
        fprintf('%s SAM \tCorrelation: %f Lag: %d\n',extents{i},timevary.(extents{i}).samcorr(nmodes),timevary.(extents{i}).samlag(nmodes));
    end
    subplot(4,1,i); hold on;
    plot(timevary.(extents{i}).samcorr,'b');
    plot(timevary.(extents{i}).curlcorr,'k');
    xlim([1 nmodes]);ylim([0.3 1]);grid on;
end

%%

% %% Correlations with climate indices
%
% sam = sam_index;
% butter.b = [0.0033 0.0164 0.0328 0.0328 0.0164 0.0033]; % 4-month butterworth
% butter.a = [1 -2.4744 2.8110 -1.7038 0.5444 -0.0723];
%
% % butter.b = [0.0006 0.0029 0.0058 0.0058 0.0029 0.0006];
% % butter.a = [1 -3.3108 4.5846 -3.2750 1.1992 -0.1793];
% sam.filtidx = filtfilt(butter.b,butter.a, sam.idx)
% plot(sam.time,sam.idx,'k-')
% plot(sam.time,sam.filtidx,'LineWidth',2)
% xlim([datenum(1992,1,1) datenum(2012,1,1)])
%
% %%
% enso = read_psd('http://www.esrl.noaa.gov/psd/gcos_wgsp/Timeseries/Data/nino34.long.data');
% enso.filtidx = filtfilt(butter.b,butter.a,enso.idx);
% %%
% extents = {'north','mean','median','south'};
% clf
% hold on
% nlags = 5;
% nmodes = 3;
% minlon = 0;
% maxlon = 360;
% maxtime = 2013;
%     mintime = 1991;
%
%     fprintf('EOF Results: %d lags, %d modes, Lonrange: %d-%d\n',nlags,nmodes,minlon,maxlon)
%     subcounter =0 ;
%     timeidx = timevary.time > datenum(mintime,1,1) & timevary.time < datenum(maxtime,1,1);
%     for i = 1:length(extents)
% %         for tidx = 1:254
% %             timevary.(extents{i}).filt_season(tidx,:) = ...
% %                 annual_harms(timevary.(extents{i}).detrended(tidx,:)',timevary.time,10,5,365.25);
% %         end
% %
%         lonidx = timevary.(extents{i}).plotlon > minlon & timevary.(extents{i}).plotlon < maxlon;
%
%         [U S V] = svd(timevary.(extents{i}).detrended(lonidx,timeidx));
%         timevary.(extents{i}).modevar = diag(S).^2/sum(diag(S).^2);
%         timevary.(extents{i}).pcs = U*S;
%
%         samidx = interp1(sam.time,sam.filtidx,timevary.time(timeidx));
%         ensoidx = interp1(enso.time,enso.filtidx,timevary.time(timeidx));
%
%         [C lags] = xcorr(sum(timevary.(extents{i}).pcs(1:nmodes,:)),ensoidx,nlags,'coeff');
%
%         [timevary.(extents{i}).ensocorr maxidx] = max(abs(C));
%         timevary.(extents{i}).ensolag = lags(maxidx);
%
%         subcounter = subcounter+1;
%         subplot(4,2,subcounter)
%         plot(lags,C,'k-x')
%         ylim([-1 1]*0.35)
%         grid on;
%
%         [C lags] = xcorr(sum(timevary.(extents{i}).pcs(1:nmodes,:)),samidx,nlags,'coeff');
%         [timevary.(extents{i}).samcorr maxidx] = max(abs(C));
%         timevary.(extents{i}).samlag = lags(maxidx);
%         timevary.(extents{i}).ensomodevar = diag(S).^2/sum(diag(S).^2);
%
%         subcounter = subcounter+1;
%         subplot(4,2,subcounter)
%         plot(lags,C,'k-x')
%         ylim([-1 1]*0.35)
%         grid on;
%
%         fprintf('%s SAM\tCorrelation: %f Lag: %d\n',extents{i},timevary.(extents{i}).samcorr,timevary.(extents{i}).samlag);
%         fprintf('%s ENSO\tCorrelation: %f Lag: %d\n',extents{i},timevary.(extents{i}).ensocorr,timevary.(extents{i}).ensolag);
%
%     end
%
% end
% %%
% clf
% subplot(2,1,1); hold on;
% varpcs = sum(timevary.south.pcs(1:nmodes,:));
% varpcs = varpcs./std(varpcs);
% % plot(timevary.time,varpcs,'k-','LineWidth',2)
% % plot(timevary.time,samidx,'LineWidth',2)
% area(timevary.time,varpcs)
% area(timevary.time,samidx)
% grid on; box on;
% title('Correlation with SAM: R=0.29 Lag = 240 days')
% xlim([datenum(1994,1,1) datenum(2012,1,1)])
% datetick('KeepLimits');
% legend('Center SAF variability','SAM Index','Location','SouthEast')
%
% subplot(2,1,2); hold on;
% varpcs = sum(timevary.south.pcs(1:nmodes,:));
% varpcs = varpcs./std(varpcs);
% plot(timevary.time,varpcs,'k-','LineWidth',2)
% plotenso = (ensoidx-mean(ensoidx))/std(ensoidx);
% plot(timevary.time,plotenso,'LineWidth',2)
% grid on; box on;
% title('Correlation with ENSO: R=0.27 Lag = 480 days')
% xlim([datenum(1994,1,1) datenum(2012,1,1)])
% datetick('KeepLimits');
% legend('Southern ACC variability','ENSO Index','Location','SouthEast')