infiles.temp = '/ltraid3/ashao/uw-apl/data/offtrac/input/normalyear/temp.hind.nc';
him.global.sst = nc_varget(infiles.temp,'temp',[0 0 0 0],[inf 1 inf inf]);
[ntime nlat nlon] = size(him.global.sst);
him.time = nc_varget(infiles.temp,'Time');
%% 
load metrics
him.npac.idx = metrics.geolon.data >= -260 & ...
    metrics.geolon.data <= -100 & ...
    metrics.geolat.data >= 20 & ...
    metrics.geolat.data <= 65 & ...
    logical(metrics.wet.data);
him.npac.Ah = metrics.Ah.data(him.npac.idx);
him.npac.wts = repmat(him.npac.Ah./sum(him.npac.Ah),[1 ntime]);
him.npac.sst = him.global.sst(:,him.npac.idx)';
him.npac.lon = metrics.geolon.data(him.npac.idx);
him.npac.lat = metrics.geolat.data(him.npac.idx);
%%
% data = detrend(him.npac.sst,'constant');
data = him.npac.sst;
data = data.*him.npac.wts;
npts = sum(him.npac.idx(:));

data = detrend(him.npac.sst')';
for i=1:npts
    if mod(i,100)==0
        fprintf('Point %d/%d\n',i,npts)
    end
    data(i,:) = annual_harms(data(i,:),him.time,10,10,365);
end
data = data-repmat(mean(data),[4374 1]);;
him.npac.sst_detrend =  data;
% 
% for i = 1:ntime
%     temp(him.npac.idx) = data(:,i);
%     m_pcolor(metrics.geolon.data,metrics.geolat.data,temp)
%     shading flat
%     m_coast('patch',[0 0 0]);m_grid;
%     colorbar
%     caxis([-15 15])
%     pause
%     
% end
%%
[U S V] = svd(data');

him.npac.pcs = U*S;

%% Plot the time series maps in first mode
Ak = U(:,1) * S(1,1) * V(:,1)';


for i=1:720
temp = nan(210,360);
    temp(him.npac.idx) = Ak(i,:);
clf
m_pcolor(metrics.geolon.data,metrics.geolat.data,temp)
shading flat
caxis([-0.05 0.05])
colorbar
m_coast('patch',[0 0 0]);m_grid;
title(sprintf('Month %d',i))
pause

end


%%
m_proj('Mercator','lon',[-260 -100],'lat',[20 65])
temp = nan(210,360);
temp(him.npac.idx) = V(:,1);
subplot(2,1,1)
m_pcolor(metrics.geolon.data,metrics.geolat.data,temp)
shading flat
caxis([-0.05 0.05])
colorbar
m_coast('patch',[0 0 0]);m_grid;
subplot(2,1,2)
plot(him.time/365+1948,him.npac.pcs(:,1))
xlim([1948 2008])
grid on;


%% Load the mixed layer depth data from hindcasts
infiles.h = dir('/ltraid4/ashao/HIM/hyak_store/HINDCAST/ocean_month*.nc');
infiles.path = '/ltraid4/ashao/HIM/hyak_store/HINDCAST/';
him.global.mldepth = zeros(720,210,360);
sidx = 1;
for i=1:length(infiles.h)
    fprintf('%d/%d ',i,length(infiles.h));
    eidx = sidx+11;
    him.global.h(sidx:eidx,:,:) = sum(nc_varget( ...
        [infiles.path infiles.h(i).name],'h', ...
        [0 0 0 0],[inf 2 inf inf]),2);
    sidx = eidx+1;
end

%%
him.npac.mldepth = him.global.mldepth(:,him.npac.idx)';
data = detrend(him.npac.mldepth')';
for i=1:npts
    if mod(i,100)==0
        fprintf('Point %d/%d\n',i,npts)
    end
    data(i,:) = annual_harms(data(i,:),him.time,10,10,365);
end
data = data-repmat(mean(data),[4374 1]);;
him.npac.mldepth_detrend =  data;