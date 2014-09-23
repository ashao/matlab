infiles.z700mb = '/ltraid4/ecmwf/era-interim/monthly/700mb.geopotential.1979.2013.nc';
ecmwf.z700mb = nc_varget(infiles.z700mb,'z');
ecmwf.time = nc_varget(infiles.z700mb,'time')/24+datenum(1900,1,1);
ecmwf.lon = nc_varget(infiles.z700mb,'longitude');
ecmwf.lat = nc_varget(infiles.z700mb,'latitude');
%%
% Calculate the Southern Annular Mode from 700mb surface as the first EOF
% between 1979-2000

[ecmwf.longrid ecmwf.latgrid] = meshgrid(ecmwf.lon,ecmwf.lat);
sidx = ecmwf.latgrid>20 & ecmwf.latgrid < 90;
samtidx = ecmwf.time > datenum(1979,1,1) & ecmwf.time < datenum(2001,1,1);

sam.zhght = ecmwf.z700mb(samtidx,sidx);
[ntime npts] = size(sam.zhght);

%% Filter out the seasonal cycle
sam.filt.zhght = detrend(sam.zhght')';

fprintf('Filtering out seasonal cycle\n')
for i=1:npts
    if mod(i,10000)==0
    fprintf('Point: %d/%d\n',i,npts);
    end
    for mon=1:12
       tidx = mon:12:ntime;
       sam.filt.zhght(tidx,i)=sam.filt.zhght(tidx,i)- ...
           mean(sam.filt.zhght(tidx,i));       
    end    
end

wts = repmat(ecmwf.latgrid(sidx)',[ntime 1]);
sam.filt.zhght=sam.filt.zhght.*wts;
%%
[U S V] = svds(double(sam.filt.zhght),1);
sam.pc = U*S;;
sam.eof = nan(size(ecmwf.longrid));
sam.eof(sidx) = V;

%% Filter all the data and project onto the mode
sozhght=ecmwf.z700mb(:,sidx);
[ntime npts] = size(sozhght);
for i=1:npts
    for mon=1:12
        tidx = mon:12:ntime;
        sozhght(tidx,i)=sozhght(tidx,i)-mean(sozhght(tidx,i));
    end
end
wts = repmat(ecmwf.latgrid(sidx)',[ntime 1]);
sozhght=sozhght.*wts;
%%
sam.pc = -sozhght*V;
colormap(othercolor('BuDRd_12'))
m_proj('stereographic','lon',0,'lat',90,'radius',20)
subplot(3,1,[1 2])
m_contourf(ecmwf.lon,ecmwf.lat,sam.eof)
m_coast('line','Color','black','linewidth',2);
m_grid;
caxis([-1 1]*2.5e-2)
colorbar

subplot(3,1,3)
area(ecmwf.time,smooth(sam.pc./std(sam.pc(samtidx)),3,'mean'))
xlim([min(ecmwf.time),max(ecmwf.time)])
datetick('KeepLimits')
ylim([-2 2])
grid on
